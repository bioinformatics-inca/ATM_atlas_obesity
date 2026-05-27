from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from surface_biomarkers.io import load_cluster_artifacts


@dataclass
class DiscoveryConfig:
    model_name: str = "LGBM"
    positive_genes: int = 1
    negative_genes: int = 1
    min_expression_prop: float = 0.10
    shap_corr_positive: float = 0.15
    shap_corr_negative: float = -0.10
    max_technical_corr: Optional[float] = None
    technical_reference_genes: Tuple[str, ...] = ()
    top_n_heatmap: int = 10


def scaled_frame(scaler: object, X: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(scaler.transform(X), columns=X.columns, index=X.index)


def compute_shap_values(model: object, X: pd.DataFrame, background: Optional[pd.DataFrame] = None) -> np.ndarray:
    import shap

    explainer = shap.TreeExplainer(model, background) if background is not None else shap.TreeExplainer(model)
    values = explainer.shap_values(X, check_additivity=False)
    return values[1] if isinstance(values, list) else values


def rank_cluster_genes(
    X_raw: pd.DataFrame,
    X_scaled: pd.DataFrame,
    y_test: pd.Series,
    shap_values: np.ndarray,
    config: Optional[DiscoveryConfig] = None,
) -> pd.DataFrame:
    config = config or DiscoveryConfig()
    y = np.asarray(y_test).astype(int)
    target_mask = y == 1
    target = X_raw.loc[target_mask]
    rows = []

    for i, gene in enumerate(X_raw.columns):
        raw = X_raw[gene].to_numpy()
        scaled = X_scaled[gene].to_numpy()
        sv = shap_values[:, i]
        sv_target = sv[target_mask]
        scaled_target = scaled[target_mask]

        expression_prop = float((target[gene] > 0).mean()) if len(target) else 0.0
        if expression_prop < config.min_expression_prop:
            continue

        shap_corr = 0.0
        if np.std(scaled_target) > 0 and np.std(sv_target) > 0:
            shap_corr = float(np.corrcoef(scaled_target, sv_target)[0, 1])
            shap_corr = 0.0 if np.isnan(shap_corr) else shap_corr

        technical_corr = 0.0
        refs = [g for g in config.technical_reference_genes if g in X_raw.columns and g != gene]
        if refs:
            corrs = []
            for ref in refs:
                if np.std(raw) > 0 and np.std(X_raw[ref].to_numpy()) > 0:
                    corr = np.corrcoef(raw, X_raw[ref].to_numpy())[0, 1]
                    if not np.isnan(corr):
                        corrs.append(abs(float(corr)))
            technical_corr = max(corrs) if corrs else 0.0

        if config.max_technical_corr is not None and technical_corr > config.max_technical_corr:
            continue

        rows.append({
            "gene": gene,
            "mean_abs_shap": float(np.abs(sv_target).mean()) if len(sv_target) else 0.0,
            "mean_shap": float(sv_target.mean()) if len(sv_target) else 0.0,
            "mean_expression_target": float(target[gene].mean()) if len(target) else 0.0,
            "expression_prop_target": expression_prop,
            "expression_shap_corr": shap_corr,
            "technical_corr": technical_corr,
        })

    if not rows:
        return pd.DataFrame(columns=[
            "gene",
            "mean_abs_shap",
            "mean_shap",
            "mean_expression_target",
            "expression_prop_target",
            "expression_shap_corr",
            "technical_corr",
        ])
    return pd.DataFrame(rows).sort_values("mean_abs_shap", ascending=False).reset_index(drop=True)


def select_signature_genes(ranking: pd.DataFrame, config: Optional[DiscoveryConfig] = None) -> List[str]:
    config = config or DiscoveryConfig()
    if ranking.empty:
        return []

    pos = ranking[ranking["expression_shap_corr"] >= config.shap_corr_positive]
    neg = ranking[ranking["expression_shap_corr"] <= config.shap_corr_negative]

    genes: List[str] = []
    genes.extend(pos.head(config.positive_genes)["gene"].tolist())
    genes.extend(neg.head(config.negative_genes)["gene"].tolist())

    if len(genes) < config.positive_genes + config.negative_genes:
        for gene in ranking["gene"]:
            if gene not in genes:
                genes.append(gene)
            if len(genes) >= config.positive_genes + config.negative_genes:
                break
    return genes


def mean_abs_shap_matrix(
    base_dir: Union[str, Path],
    clusters: List[str],
    config: Optional[DiscoveryConfig] = None,
) -> pd.DataFrame:
    config = config or DiscoveryConfig()
    values = {}
    top_features = set()

    for cluster in clusters:
        artifacts = load_cluster_artifacts(Path(base_dir) / str(cluster), str(cluster), config.model_name)
        X_train_sc = scaled_frame(artifacts.scaler, artifacts.X_train)
        X_test_sc = scaled_frame(artifacts.scaler, artifacts.X_test)
        sv = compute_shap_values(artifacts.model, X_test_sc, X_train_sc)
        series = pd.Series(np.abs(sv).mean(axis=0), index=X_test_sc.columns)
        values[f"Cluster {cluster}"] = series
        top_features.update(series.nlargest(config.top_n_heatmap).index)

    heatmap = pd.DataFrame(values).loc[sorted(top_features)]
    return heatmap


def combined_shap_expression_matrix(
    base_dir: Union[str, Path],
    clusters: List[str],
    config: Optional[DiscoveryConfig] = None,
) -> pd.DataFrame:
    config = config or DiscoveryConfig()
    shap_df = mean_abs_shap_matrix(base_dir, clusters, config)
    expr = {}
    for cluster in clusters:
        artifacts = load_cluster_artifacts(Path(base_dir) / str(cluster), str(cluster), config.model_name)
        expr[f"Cluster {cluster}"] = artifacts.X_test.mean(axis=0)
    expr_df = pd.DataFrame(expr).reindex(shap_df.index)
    expr_norm = expr_df.div(expr_df.max(axis=1).replace(0, np.nan), axis=0).fillna(0)
    return shap_df * expr_norm


def calculate_exclusivity(score_matrix: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene, row in score_matrix.iterrows():
        ordered = row.sort_values(ascending=False)
        best = ordered.index[0]
        second = ordered.iloc[1] if len(ordered) > 1 else 0.0
        rows.append({
            "gene": gene,
            "best_cluster": best,
            "best_value": float(ordered.iloc[0]),
            "second_best": float(second),
            "exclusivity": float(ordered.iloc[0] - second),
        })
    return pd.DataFrame(rows).sort_values("exclusivity", ascending=False).reset_index(drop=True)


def discover_signatures(
    base_dir: Union[str, Path],
    clusters: List[str],
    config: Optional[DiscoveryConfig] = None,
) -> Tuple[Dict[str, List[str]], pd.DataFrame, pd.DataFrame]:
    config = config or DiscoveryConfig()
    signatures: Dict[str, List[str]] = {}
    rank_tables = []

    for cluster in clusters:
        artifacts = load_cluster_artifacts(Path(base_dir) / str(cluster), str(cluster), config.model_name)
        X_train_sc = scaled_frame(artifacts.scaler, artifacts.X_train)
        X_test_sc = scaled_frame(artifacts.scaler, artifacts.X_test)
        sv = compute_shap_values(artifacts.model, X_test_sc, X_train_sc)
        ranking = rank_cluster_genes(artifacts.X_test, X_test_sc, artifacts.y_test, sv, config)
        ranking.insert(0, "cluster", f"Cluster {cluster}")
        rank_tables.append(ranking)
        signatures[f"Cluster {cluster}"] = select_signature_genes(ranking, config)

    combined = combined_shap_expression_matrix(base_dir, clusters, config)
    exclusivity = calculate_exclusivity(combined)
    return signatures, pd.concat(rank_tables, ignore_index=True), exclusivity


def load_label_mapping(base_dir: Union[str, Path]) -> Dict[str, str]:
    mapping_path = Path(base_dir) / "label_mapping.csv"
    if not mapping_path.exists():
        return {}
    mapping_df = pd.read_csv(mapping_path)
    if "encoded_label" not in mapping_df.columns or "original_label" not in mapping_df.columns:
        return {}
    mapping = {}
    for _, row in mapping_df.iterrows():
        encoded = str(row["encoded_label"])
        original = str(row["original_label"])
        mapping[encoded] = original
        mapping[f"Cluster {encoded}"] = original
    return mapping


def _resolve_cluster_label(
    cluster: str,
    signatures: Dict[str, List[str]],
    label_mapping: Dict[str, str],
) -> str:
    encoded = str(cluster)
    cluster_name = f"Cluster {encoded}"
    real_name = label_mapping.get(encoded, label_mapping.get(cluster_name))
    if real_name in signatures:
        return real_name
    if cluster_name in signatures:
        return cluster_name
    if encoded in signatures:
        return encoded
    return real_name or cluster_name


def _normalize_cluster_names(cluster_names: Optional[Dict[str, str]]) -> Dict[str, str]:
    if not cluster_names:
        return {}
    normalized = {}
    for key, value in cluster_names.items():
        key_text = str(key)
        value_text = str(value)
        normalized[key_text] = value_text
        if key_text.startswith("Cluster "):
            normalized[key_text.replace("Cluster ", "", 1)] = value_text
        else:
            normalized[f"Cluster {key_text}"] = value_text
    return normalized


def evaluate_signatures(
    base_dir: Union[str, Path],
    clusters: List[str],
    signatures: Dict[str, List[str]],
    threshold: float = 0.5,
    cluster_map: Optional[Dict[str, str]] = None,
    cluster_names: Optional[Dict[str, str]] = None,
) -> pd.DataFrame:
    """Evaluate small gene signatures using the original X_test block logic.

    This intentionally mirrors the notebook workflow:
    1. read each ``cluster/X_test.csv``;
    2. assign every row in that file to the label from ``cluster_map``;
    3. score each signature by summing valid genes in ``X_test``;
    4. compute independent binary metrics for each signature.

    It does not read ``y_test.csv``.
    """
    from sklearn.metrics import accuracy_score, recall_score, roc_auc_score

    cluster_names_map = _normalize_cluster_names(cluster_names)
    all_scores = []
    all_y_true = []

    for cluster in clusters:
        folder = Path(base_dir) / str(cluster)
        path_test = folder / "X_test.csv"
        if not path_test.exists():
            continue

        X_test = pd.read_csv(path_test, index_col=0)
        if cluster_map is not None:
            mac_label = cluster_map[f"Cluster {cluster}"]
        elif str(cluster) in cluster_names_map and cluster_names_map[str(cluster)] in signatures:
            mac_label = cluster_names_map[str(cluster)]
        elif f"Cluster {cluster}" in cluster_names_map and cluster_names_map[f"Cluster {cluster}"] in signatures:
            mac_label = cluster_names_map[f"Cluster {cluster}"]
        elif f"Cluster {cluster}" in signatures:
            mac_label = f"Cluster {cluster}"
        elif str(cluster) in signatures:
            mac_label = str(cluster)
        else:
            mac_label = f"Cluster {cluster}"
        all_y_true.extend([mac_label] * len(X_test))

        batch_scores = pd.DataFrame(index=X_test.index)
        for signature_label, genes in signatures.items():
            valid_genes = [gene for gene in genes if gene in X_test.columns]
            batch_scores[signature_label] = X_test[valid_genes].sum(axis=1) if valid_genes else 0.0
        all_scores.append(batch_scores)

    if not all_scores:
        return pd.DataFrame()

    df_scores = pd.concat(all_scores)
    y_true_array = np.asarray(all_y_true)
    rows = []

    for target_label in signatures.keys():
        display_label = target_label
        if cluster_names_map:
            display_label = cluster_names_map.get(str(target_label), target_label)
        y_binary = (y_true_array == target_label).astype(int)
        y_score = df_scores[target_label]
        y_pred = (y_score > threshold).astype(int)
        auc = roc_auc_score(y_binary, y_score) if len(np.unique(y_binary)) > 1 else np.nan
        acc = round(accuracy_score(y_binary, y_pred), 4)
        rows.append({
            "Mac": target_label,
            "cluster": target_label,
            "cluster_label": display_label,
            "Acuracia (Indep)": acc,
            "accuracy": acc,
            "N_Positivos": int(y_binary.sum()),
            "auc": auc,
            "sensitivity": recall_score(y_binary, y_pred, zero_division=0),
            "threshold": threshold,
        })
    return pd.DataFrame(rows)
