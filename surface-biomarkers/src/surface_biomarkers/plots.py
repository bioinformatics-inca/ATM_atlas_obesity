from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)

from surface_biomarkers.io import load_cluster_artifacts
from surface_biomarkers.signatures import compute_shap_values, scaled_frame


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


def display_cluster_name(cluster: str, label_mapping: Optional[Dict[str, str]] = None) -> str:
    if not label_mapping:
        return cluster
    return label_mapping.get(str(cluster), label_mapping.get(f"Cluster {cluster}", str(cluster)))


def plot_shap_heatmap(
    matrix: pd.DataFrame,
    output_path: Union[str, Path],
    normalize_by_gene: bool = True,
    title: Optional[str] = None,
    label_mapping: Optional[Dict[str, str]] = None,
    cluster_names: Optional[Dict[str, str]] = None,
    show: bool = True,
) -> Path:
    plot_df = matrix.T
    annot_df = plot_df.copy()
    if label_mapping is None:
        label_mapping = load_label_mapping(Path(output_path).parent)
    if cluster_names:
        label_mapping.update({str(key): str(value) for key, value in cluster_names.items()})
        label_mapping.update({f"Cluster {key}": str(value) for key, value in cluster_names.items()})
    plot_df.index = [display_cluster_name(idx, label_mapping) for idx in plot_df.index]
    annot_df.index = plot_df.index
    if normalize_by_gene:
        plot_df = plot_df.div(plot_df.max(axis=0).replace(0, pd.NA), axis=1).fillna(0)

    fig, ax = plt.subplots(figsize=(max(10, matrix.shape[0] * 0.7), max(4, matrix.shape[1] * 0.7)))
    sns.heatmap(
        plot_df,
        ax=ax,
        cmap="vlag",
        linewidths=0.4,
        linecolor="#eeeeee",
        annot=annot_df.round(3),
        annot_kws={"size": 9},
        fmt=".3f",
        cbar_kws={"label": "Normalized score" if normalize_by_gene else "Score", "shrink": 0.7},
        vmin=0 if normalize_by_gene else None,
        vmax=1 if normalize_by_gene else None,
    )
    if title:
        ax.set_title(title, fontsize=12, fontweight="bold", pad=12)
    ax.set_xlabel("Gene")
    ax.set_ylabel("Cluster")
    ax.tick_params(axis="x", labelsize=10, rotation=45)
    ax.tick_params(axis="y", labelsize=10, rotation=0)
    plt.tight_layout()
    output = Path(output_path)
    fig.savefig(output, bbox_inches="tight", dpi=300)
    if show:
        plt.show()
    plt.close(fig)
    return output


def plot_signature_accuracy(
    results: pd.DataFrame,
    output_path: Union[str, Path],
    cluster_names: Optional[Dict[str, str]] = None,
    show: bool = True,
) -> Path:
    fig, ax = plt.subplots(figsize=(max(8, len(results) * 1.8), 3))
    if "cluster_label" in results.columns:
        labels = results["cluster_label"].astype(str)
    elif "Mac" in results.columns:
        labels = results["Mac"].astype(str)
    else:
        labels = results["cluster"].astype(str)
    if cluster_names:
        mapping = {str(key): str(value) for key, value in cluster_names.items()}
        mapping.update({f"Cluster {key}": str(value) for key, value in cluster_names.items()})
        labels = labels.map(lambda label: mapping.get(str(label), str(label)))

    value_col = "accuracy" if "accuracy" in results.columns else "Acuracia (Indep)"
    values = results[value_col]
    ax.bar(labels, values, color="#6A5D7B", edgecolor="black", linewidth=0.5)
    for i, value in enumerate(values):
        ax.text(i, value + 0.015, f"{value:.2f}", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.set_ylabel("Accuracy")
    ax.set_ylim(0, 1.05)
    ax.grid(axis="y", alpha=0.15)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="x", rotation=30)
    plt.tight_layout()
    output = Path(output_path)
    fig.savefig(output, bbox_inches="tight", dpi=300)
    if show:
        plt.show()
    plt.close(fig)
    return output


def plot_panel_benchmark(
    panel_results: pd.DataFrame,
    output_path: Union[str, Path],
    cluster_names: Optional[Dict[str, str]] = None,
    title: str = "Antibody Panel Benchmark",
    show: bool = True,
) -> Path:
    required = {
        "cluster",
        "application",
        "panel_score",
        "panel_size",
        "backup_count",
        "redundancy_penalty",
        "constraints_passed",
    }
    missing = sorted(required - set(panel_results.columns))
    if missing:
        raise ValueError("panel_results is missing required columns: " + ", ".join(missing))

    plot_df = panel_results.copy()
    plot_df["cluster"] = plot_df["cluster"].astype(str)
    plot_df["application"] = plot_df["application"].astype(str)
    if cluster_names:
        mapping = {str(key): str(value) for key, value in cluster_names.items()}
        mapping.update({f"Cluster {key}": str(value) for key, value in cluster_names.items()})
        plot_df["cluster_label"] = plot_df["cluster"].map(lambda value: mapping.get(value, value))
    else:
        plot_df["cluster_label"] = plot_df["cluster"]
    plot_df["label"] = plot_df["cluster_label"] + "\n" + plot_df["application"]
    plot_df["panel_score"] = pd.to_numeric(plot_df["panel_score"], errors="coerce").fillna(0.0)
    plot_df["panel_size"] = pd.to_numeric(plot_df["panel_size"], errors="coerce").fillna(0).astype(int)
    plot_df["backup_count"] = pd.to_numeric(plot_df["backup_count"], errors="coerce").fillna(0).astype(int)
    plot_df["redundancy_penalty"] = pd.to_numeric(plot_df["redundancy_penalty"], errors="coerce").fillna(0.0)
    plot_df["constraints_passed"] = plot_df["constraints_passed"].astype(bool)
    plot_df = plot_df.sort_values(["constraints_passed", "panel_score"], ascending=[False, False])

    x = list(range(len(plot_df)))
    labels = plot_df["label"].tolist()
    colors = ["#3F7D5C" if passed else "#B75555" for passed in plot_df["constraints_passed"]]
    width = max(10, len(plot_df) * 1.2)
    fig, axes = plt.subplots(2, 2, figsize=(width, 8))
    axes = axes.ravel()

    axes[0].bar(x, plot_df["panel_score"], color=colors, edgecolor="black", linewidth=0.4)
    axes[0].set_title("Panel score", fontsize=11, fontweight="bold")
    axes[0].set_ylabel("Score")
    axes[0].grid(axis="y", alpha=0.18)

    axes[1].bar(x, plot_df["panel_size"], color="#5D7599", edgecolor="black", linewidth=0.4, label="Minimal panel")
    axes[1].bar(
        x,
        plot_df["backup_count"],
        bottom=plot_df["panel_size"],
        color="#C9A646",
        edgecolor="black",
        linewidth=0.4,
        label="Backup markers",
    )
    axes[1].set_title("Panel size and backups", fontsize=11, fontweight="bold")
    axes[1].set_ylabel("Marker count")
    axes[1].legend(frameon=False, fontsize=9)
    axes[1].grid(axis="y", alpha=0.18)

    axes[2].bar(x, plot_df["redundancy_penalty"], color="#7D688E", edgecolor="black", linewidth=0.4)
    axes[2].set_title("Redundancy penalty", fontsize=11, fontweight="bold")
    axes[2].set_ylabel("Mean marker redundancy")
    axes[2].grid(axis="y", alpha=0.18)

    pass_values = plot_df["constraints_passed"].astype(int)
    axes[3].bar(x, pass_values, color=colors, edgecolor="black", linewidth=0.4)
    axes[3].set_title("Application constraints", fontsize=11, fontweight="bold")
    axes[3].set_ylabel("Passed")
    axes[3].set_yticks([0, 1])
    axes[3].set_yticklabels(["No", "Yes"])
    axes[3].set_ylim(0, 1.2)
    axes[3].grid(axis="y", alpha=0.18)

    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=9)
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

    fig.suptitle(title, fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    output = Path(output_path)
    fig.savefig(output, bbox_inches="tight", dpi=300)
    if show:
        plt.show()
    plt.close(fig)
    return output


def cluster_metrics(y_true, y_pred, y_prob) -> Dict[str, float]:
    cm = confusion_matrix(y_true, y_pred)
    tn, fp, fn, tp = cm.ravel()
    specificity = tn / (tn + fp) if (tn + fp) else 0.0
    return {
        "auc": roc_auc_score(y_true, y_prob),
        "accuracy": accuracy_score(y_true, y_pred),
        "sensitivity": recall_score(y_true, y_pred, zero_division=0),
        "specificity": specificity,
        "precision": precision_score(y_true, y_pred, zero_division=0),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "tp": int(tp),
    }


def plot_confusion_matrix(
    y_true,
    y_pred,
    metrics: Dict[str, float],
    cluster: str,
    model_name: str,
    output_dir: Union[str, Path],
    display_cluster: Optional[str] = None,
    show: bool = True,
) -> Path:
    label = display_cluster or cluster
    cm = confusion_matrix(y_true, y_pred)
    fig, ax = plt.subplots(figsize=(5, 4.5))
    sns.heatmap(
        cm,
        annot=True,
        fmt="d",
        cmap="binary",
        xticklabels=["Others", label],
        yticklabels=["Others", label],
        annot_kws={"size": 12, "weight": "bold"},
        ax=ax,
    )
    ax.set_title(f"Confusion Matrix - {model_name}\n({label} vs All)", fontsize=11)
    ax.set_xlabel("Predicted label")
    ax.set_ylabel("True label")
    fig.text(
        0.5,
        0.01,
        (
            f"AUC: {metrics['auc']:.3f} | Accuracy: {metrics['accuracy']:.3f} | "
            f"Sensitivity: {metrics['sensitivity']:.3f} | Specificity: {metrics['specificity']:.3f}"
        ),
        ha="center",
        va="bottom",
        fontsize=9,
    )
    plt.tight_layout(rect=[0, 0.07, 1, 1])
    path = Path(output_dir) / f"confusion_{model_name}.svg"
    fig.savefig(path, bbox_inches="tight", dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return path


def plot_roc(
    y_true,
    y_prob,
    cluster: str,
    model_name: str,
    output_dir: Union[str, Path],
    display_cluster: Optional[str] = None,
    show: bool = True,
) -> Path:
    label = display_cluster or cluster
    fpr, tpr, _ = roc_curve(y_true, y_prob)
    auc = roc_auc_score(y_true, y_prob)
    fig, ax = plt.subplots(figsize=(5, 4.5))
    ax.fill_between(fpr, tpr, alpha=0.12, color="#2c7bb6")
    ax.plot(fpr, tpr, color="#2c7bb6", lw=2.5, label=f"AUC = {auc:.3f}")
    ax.plot([0, 1], [0, 1], color="#aaaaaa", lw=1.5, linestyle="--")
    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.05])
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title(f"ROC Curve - {model_name}\n({label} vs All)", fontsize=11)
    ax.legend(loc="lower right")
    ax.grid(True, linestyle=":", linewidth=0.6, color="#cccccc")
    plt.tight_layout()
    path = Path(output_dir) / f"roc_{model_name}.svg"
    fig.savefig(path, bbox_inches="tight", dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return path


def plot_shap_summary(
    model,
    X_train_scaled: pd.DataFrame,
    X_test_scaled: pd.DataFrame,
    cluster: str,
    model_name: str,
    output_dir: Union[str, Path],
    max_display: int = 10,
    display_cluster: Optional[str] = None,
    show: bool = True,
) -> Path:
    import shap
    from matplotlib.colors import LinearSegmentedColormap

    label = display_cluster or cluster
    cmap = LinearSegmentedColormap.from_list("shap_custom", ["#357EBDFF", "#D43F3AFF"], N=15)
    shap_values = compute_shap_values(model, X_test_scaled, X_train_scaled)
    plt.figure(figsize=(5, 6))
    shap.summary_plot(
        shap_values,
        X_test_scaled,
        max_display=max_display,
        show=False,
        plot_size=None,
        cmap=cmap,
    )
    fig = plt.gcf()
    fig.suptitle(
        f"SHAP Summary Plot - {model_name}\n({label} vs All)",
        fontsize=10,
        fontweight="bold",
        y=1.02,
    )
    plt.tight_layout()
    path = Path(output_dir) / f"shap_summary_{model_name}.svg"
    fig.savefig(path, bbox_inches="tight", dpi=150, facecolor="white")
    if show:
        plt.show()
    plt.close(fig)
    return path


def plot_cluster_panel(
    model,
    X_train_scaled: pd.DataFrame,
    X_test_scaled: pd.DataFrame,
    y_true,
    y_pred,
    y_prob,
    metrics: Dict[str, float],
    cluster: str,
    model_name: str,
    output_dir: Union[str, Path],
    max_display: int = 10,
    display_cluster: Optional[str] = None,
    show: bool = True,
) -> Path:
    import shap
    from matplotlib.colors import LinearSegmentedColormap

    label = display_cluster or cluster
    cmap = LinearSegmentedColormap.from_list("shap_custom", ["#357EBDFF", "#D43F3AFF"], N=15)
    shap_values = compute_shap_values(model, X_test_scaled, X_train_scaled)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))

    cm = confusion_matrix(y_true, y_pred)
    sns.heatmap(
        cm,
        annot=True,
        fmt="d",
        cmap="binary",
        xticklabels=["Others", label],
        yticklabels=["Others", label],
        annot_kws={"size": 11, "weight": "bold"},
        ax=axes[0],
    )
    axes[0].set_title("Confusion matrix", fontsize=11, fontweight="bold")
    axes[0].set_xlabel("Predicted")
    axes[0].set_ylabel("True")

    fpr, tpr, _ = roc_curve(y_true, y_prob)
    axes[1].fill_between(fpr, tpr, alpha=0.12, color="#2c7bb6")
    axes[1].plot(fpr, tpr, color="#2c7bb6", lw=2.4, label=f"AUC = {metrics['auc']:.3f}")
    axes[1].plot([0, 1], [0, 1], color="#aaaaaa", lw=1.2, linestyle="--")
    axes[1].set_xlim([-0.02, 1.02])
    axes[1].set_ylim([-0.02, 1.05])
    axes[1].set_title("ROC curve", fontsize=11, fontweight="bold")
    axes[1].set_xlabel("False Positive Rate")
    axes[1].set_ylabel("True Positive Rate")
    axes[1].legend(loc="lower right", fontsize=9)
    axes[1].grid(True, linestyle=":", linewidth=0.6, color="#cccccc")

    plt.sca(axes[2])
    shap.summary_plot(
        shap_values,
        X_test_scaled,
        max_display=max_display,
        show=False,
        plot_size=None,
        cmap=cmap,
        color_bar=False,
    )
    axes[2].set_title("SHAP summary", fontsize=11, fontweight="bold")

    fig.suptitle(
        (
            f"{label} vs All - {model_name} | "
            f"AUC={metrics['auc']:.3f}, Acc={metrics['accuracy']:.3f}, "
            f"Sens={metrics['sensitivity']:.3f}, Spec={metrics['specificity']:.3f}"
        ),
        fontsize=12,
        fontweight="bold",
    )
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    path = Path(output_dir) / f"panel_{model_name}.svg"
    fig.savefig(path, bbox_inches="tight", dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return path


def analyze_cluster(
    cluster_dir: Union[str, Path],
    cluster: Optional[str] = None,
    model_name: str = "LGBM",
    max_display: int = 10,
    display_cluster: Optional[str] = None,
    panel: bool = True,
    show_panel: bool = True,
) -> Tuple[Dict[str, float], List[Path]]:
    artifacts = load_cluster_artifacts(cluster_dir, cluster, model_name)
    X_train_scaled = scaled_frame(artifacts.scaler, artifacts.X_train)
    X_test_scaled = scaled_frame(artifacts.scaler, artifacts.X_test)
    y_pred = artifacts.model.predict(X_test_scaled)
    y_prob = artifacts.model.predict_proba(X_test_scaled)[:, 1]
    metrics = cluster_metrics(artifacts.y_test, y_pred, y_prob)
    metrics["cluster"] = artifacts.cluster
    metrics["cluster_label"] = display_cluster or artifacts.cluster
    metrics["model"] = model_name

    paths = [
        plot_confusion_matrix(
            artifacts.y_test,
            y_pred,
            metrics,
            artifacts.cluster,
            model_name,
            artifacts.path,
            display_cluster=display_cluster,
            show=False,
        ),
        plot_roc(
            artifacts.y_test,
            y_prob,
            artifacts.cluster,
            model_name,
            artifacts.path,
            display_cluster=display_cluster,
            show=False,
        ),
        plot_shap_summary(
            artifacts.model,
            X_train_scaled,
            X_test_scaled,
            artifacts.cluster,
            model_name,
            artifacts.path,
            max_display=max_display,
            display_cluster=display_cluster,
            show=False,
        ),
    ]
    if panel:
        paths.append(
            plot_cluster_panel(
                artifacts.model,
                X_train_scaled,
                X_test_scaled,
                artifacts.y_test,
                y_pred,
                y_prob,
                metrics,
                artifacts.cluster,
                model_name,
                artifacts.path,
                max_display=max_display,
                display_cluster=display_cluster,
                show=show_panel,
            )
        )
    pd.DataFrame([metrics]).to_csv(artifacts.path / f"metrics_{model_name}.csv", index=False)
    return metrics, paths


def analyze_clusters(
    base_dir: Union[str, Path],
    clusters: List[str],
    model_name: Union[str, List[str], Tuple[str, ...]] = "LGBM",
    max_display: int = 10,
    cluster_names: Optional[Dict[str, str]] = None,
    panel: bool = True,
    show_panel: bool = True,
) -> pd.DataFrame:
    model_names = [model_name] if isinstance(model_name, str) else list(model_name)
    label_mapping = load_label_mapping(base_dir)
    if cluster_names:
        label_mapping.update({str(key): str(value) for key, value in cluster_names.items()})
        label_mapping.update({f"Cluster {key}": str(value) for key, value in cluster_names.items()})
    rows = []
    for current_model in model_names:
        for cluster in clusters:
            cluster_label = display_cluster_name(str(cluster), label_mapping)
            metrics, _ = analyze_cluster(
                Path(base_dir) / str(cluster),
                cluster=str(cluster),
                model_name=current_model,
                max_display=max_display,
                display_cluster=cluster_label,
                panel=panel,
                show_panel=show_panel,
            )
            rows.append(metrics)
    metrics_df = pd.DataFrame(rows)
    ordered = [
        "cluster",
        "cluster_label",
        "model",
        "auc",
        "accuracy",
        "sensitivity",
        "specificity",
        "precision",
        "tn",
        "fp",
        "fn",
        "tp",
    ]
    metrics_df = metrics_df[[col for col in ordered if col in metrics_df.columns]]
    suffix = model_names[0] if len(model_names) == 1 else "all_models"
    metrics_df.to_csv(Path(base_dir) / f"metrics_all_clusters_{suffix}.csv", index=False)
    return metrics_df


def plot_confusion_and_roc(
    cluster_dir: Union[str, Path],
    cluster: Optional[str] = None,
    model_name: str = "LGBM",
    show: bool = True,
) -> List[Path]:
    metrics, paths = analyze_cluster(cluster_dir, cluster=cluster, model_name=model_name, panel=False, show_panel=False)
    if show:
        artifacts = load_cluster_artifacts(cluster_dir, cluster, model_name)
        X_test_scaled = scaled_frame(artifacts.scaler, artifacts.X_test)
        y_pred = artifacts.model.predict(X_test_scaled)
        y_prob = artifacts.model.predict_proba(X_test_scaled)[:, 1]
        plot_confusion_matrix(artifacts.y_test, y_pred, metrics, artifacts.cluster, model_name, artifacts.path, show=True)
        plot_roc(artifacts.y_test, y_prob, artifacts.cluster, model_name, artifacts.path, show=True)
    return paths[:2]
