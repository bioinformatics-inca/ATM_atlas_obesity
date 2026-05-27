from typing import Dict, List, Optional

import pandas as pd


def design_antibody_panel(
    markers: pd.DataFrame,
    application: str = "FACS",
    cluster_column: str = "cluster",
    gene_column: str = "gene",
    max_positive_markers: int = 2,
    max_negative_markers: int = 1,
    max_backup_markers: int = 2,
    min_positive_markers: int = 1,
    redundancy_matrix: Optional[pd.DataFrame] = None,
    redundancy_threshold: float = 0.85,
) -> pd.DataFrame:
    """Recommend compact antibody-marker panels from prioritized marker rows.

    This function designs marker combinations, not only individual-gene ranks.
    It expects the output of ``prioritize_signature_validation`` but also works
    with any DataFrame containing compatible score columns. Application-specific
    constraints are intentionally conservative: snRNA-seq evidence can nominate
    candidates, while protein localization and antibody evidence decide whether
    a marker is suitable for antibody-based assays.
    """
    required = {cluster_column, gene_column}
    missing = [column for column in required if column not in markers.columns]
    if missing:
        raise ValueError("markers is missing required columns: " + ", ".join(missing))

    def numeric(df: pd.DataFrame, column: str, default: float = 0.0) -> pd.Series:
        if column not in df.columns:
            return pd.Series(default, index=df.index, dtype=float)
        return pd.to_numeric(df[column], errors="coerce").fillna(default)

    def text_list(values: List[str]) -> str:
        return "; ".join(values)

    def redundancy(gene_a: str, gene_b: str) -> float:
        if redundancy_matrix is None:
            return 0.0
        if gene_a in redundancy_matrix.index and gene_b in redundancy_matrix.columns:
            value = redundancy_matrix.loc[gene_a, gene_b]
        elif gene_b in redundancy_matrix.index and gene_a in redundancy_matrix.columns:
            value = redundancy_matrix.loc[gene_b, gene_a]
        else:
            value = 0.0
        return float(pd.to_numeric(pd.Series([value]), errors="coerce").fillna(0).iloc[0])

    def mean_redundancy(genes: List[str]) -> float:
        if len(genes) < 2:
            return 0.0
        values = []
        for i, gene_a in enumerate(genes):
            for gene_b in genes[i + 1:]:
                values.append(abs(redundancy(gene_a, gene_b)))
        return float(sum(values) / len(values)) if values else 0.0

    app = application.strip().lower()
    constraints: Dict[str, Dict[str, float]] = {
        "facs": {"localization": 5.0, "antibody": 3.0, "application": 5.0},
        "flow": {"localization": 5.0, "antibody": 3.0, "application": 5.0},
        "flow cytometry": {"localization": 5.0, "antibody": 3.0, "application": 5.0},
        "cytof": {"localization": 5.0, "antibody": 4.0, "application": 5.0},
        "cite-seq": {"localization": 5.0, "antibody": 4.0, "application": 5.0},
        "citeseq": {"localization": 5.0, "antibody": 4.0, "application": 5.0},
        "if": {"localization": 2.0, "antibody": 4.0, "application": 4.0},
        "spatial": {"localization": 1.0, "antibody": 4.0, "application": 4.0},
        "ihc": {"localization": 1.0, "antibody": 4.0, "application": 4.0},
    }
    if app not in constraints:
        raise ValueError(
            "application must be one of: FACS, CyTOF, CITE-seq, IF, spatial, IHC"
        )
    limits = constraints[app]

    data = markers.copy()
    data[cluster_column] = data[cluster_column].astype(str)
    data[gene_column] = data[gene_column].astype(str)
    data["_final_marker_score"] = numeric(data, "final_marker_score")
    data["_transcript_score"] = numeric(data, "transcript_score")
    data["_specificity_score"] = numeric(data, "specificity_score")
    data["_prevalence_score"] = numeric(data, "prevalence_score")
    data["_localization_score"] = numeric(data, "localization_score")
    data["_antibody_availability_score"] = numeric(data, "antibody_availability_score")
    data["_application_compatibility_score"] = numeric(data, "application_compatibility_score")
    data["_penalty_score"] = numeric(data, "penalty_score")
    data["_mean_shap"] = numeric(data, "mean_shap")
    data["_expression_shap_corr"] = numeric(data, "expression_shap_corr")
    data["_is_signature"] = data["is_signature"].fillna(False).astype(bool) if "is_signature" in data.columns else True

    data["_passes_application"] = (
        (data["_localization_score"] >= limits["localization"])
        & (data["_antibody_availability_score"] >= limits["antibody"])
        & (data["_application_compatibility_score"] >= limits["application"])
    )
    data["_panel_rank_score"] = (
        data["_final_marker_score"]
        + 0.50 * data["_specificity_score"]
        + 0.25 * data["_prevalence_score"]
        - data["_penalty_score"]
    )

    rows = []
    for cluster, cluster_df in data.groupby(cluster_column, sort=False):
        candidates = cluster_df[cluster_df["_is_signature"]].copy()
        if candidates.empty:
            candidates = cluster_df.copy()

        positive_pool = candidates[
            candidates["_passes_application"]
            & (candidates["_mean_shap"] >= 0)
            & (candidates["_expression_shap_corr"] >= 0)
        ].sort_values("_panel_rank_score", ascending=False)
        negative_pool = candidates[
            candidates["_passes_application"]
            & (
                (candidates["_mean_shap"] < 0)
                | (candidates["_expression_shap_corr"] < 0)
            )
        ].sort_values("_panel_rank_score", ascending=False)

        selected_positive: List[str] = []
        redundancy_warnings: List[str] = []
        for _, row in positive_pool.iterrows():
            gene = row[gene_column]
            if any(abs(redundancy(gene, chosen)) >= redundancy_threshold for chosen in selected_positive):
                redundancy_warnings.append(gene)
                continue
            selected_positive.append(gene)
            if len(selected_positive) >= max_positive_markers:
                break

        selected_negative: List[str] = []
        for _, row in negative_pool.iterrows():
            gene = row[gene_column]
            if gene in selected_positive:
                continue
            if any(abs(redundancy(gene, chosen)) >= redundancy_threshold for chosen in selected_positive + selected_negative):
                redundancy_warnings.append(gene)
                continue
            selected_negative.append(gene)
            if len(selected_negative) >= max_negative_markers:
                break

        selected = selected_positive + selected_negative
        backup_pool = candidates[
            candidates["_passes_application"]
            & (~candidates[gene_column].isin(selected))
        ].sort_values("_panel_rank_score", ascending=False)
        backup_markers = backup_pool[gene_column].head(max_backup_markers).tolist()

        panel_genes = selected + backup_markers
        panel_rows = candidates[candidates[gene_column].isin(selected)]
        redundancy_penalty = mean_redundancy(selected)
        panel_score = (
            float(panel_rows["_panel_rank_score"].mean()) if not panel_rows.empty else 0.0
        ) - redundancy_penalty

        warnings = []
        if len(selected_positive) < min_positive_markers:
            warnings.append("insufficient_positive_markers")
        if not selected_negative:
            warnings.append("no_negative_exclusion_marker")
        if redundancy_warnings:
            warnings.append("redundant_markers_skipped:" + ",".join(sorted(set(redundancy_warnings))))
        if not backup_markers:
            warnings.append("no_backup_markers")
        if candidates["_passes_application"].sum() == 0:
            warnings.append("no_marker_passed_application_constraints")

        rows.append({
            "cluster": cluster,
            "application": application,
            "positive_markers": text_list(selected_positive),
            "negative_exclusion_markers": text_list(selected_negative),
            "backup_markers": text_list(backup_markers),
            "minimal_panel": text_list(selected),
            "panel_size": len(selected),
            "backup_count": len(backup_markers),
            "redundancy_penalty": round(float(redundancy_penalty), 4),
            "panel_score": round(float(panel_score), 4),
            "constraints_passed": len(selected_positive) >= min_positive_markers and candidates["_passes_application"].sum() > 0,
            "warnings": text_list(warnings),
            "all_candidate_markers": text_list(panel_genes),
        })

    return pd.DataFrame(rows).sort_values(
        ["constraints_passed", "panel_score"],
        ascending=[False, False],
    ).reset_index(drop=True)
