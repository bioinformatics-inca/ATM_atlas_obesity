# surface-biomarkers

`surface-biomarkers` is a Python package for discovering compact surface-marker signatures from single-nuclei expression data.


- one-vs-rest machine-learning models for each annotated cluster;
- per-cluster ROC curves, confusion matrices, SHAP plots and metric tables;
- compact marker signatures using SHAP and expression direction;
- heatmaps of cluster-specific SHAP;
- independent accuracy estimates for small gene panels.

The package is designed for interactive use in Jupyter notebooks, while still providing a command-line interface for reproducible runs.

## Workflow Overview

```text
counts.csv
  |
  | load_counts_csv()
  v
cells x surface genes + celltype
  |
  | train_one_vs_rest()
  v
cluster_outputs/
  |
  | analyze_clusters()
  | discover_signatures()
  | evaluate_signatures()
  v
plots, metrics, SHAP rankings, marker signatures
```

## Installation

### Jupyter, simple install

```python
%pip install -e "ML_surface/surface-biomarkers"
```

### Jupyter, reproducible install

Use this when you want to reproduce the exact environment used during development:

```python
%pip install -r "ML_surface/surface-biomarkers/requirements.txt"
%pip install -e "ML_surface/surface-biomarkers" --no-deps
```

Restart the kernel after installation.

Pinned dependency versions:

```text
pandas==1.3.5
numpy==1.18.2
scikit-learn==0.22.2.post1
imbalanced-learn==0.6.2
lightgbm==2.3.1
matplotlib==3.0.2
seaborn==0.9.0
shap==0.42.1
```

## Input Format

The input file should be a CSV with cells in rows, genes in columns, and one annotation column:

```text
cell_id,PTPRC,MRC1,ABCA1,PLXND1,celltype
AAAC,1,3,0,2,0
AAAG,0,1,4,0,1
```

The annotation column is usually named `celltype`, but this can be changed with `target_column`.

## Loading Data

```python
from surface_biomarkers import load_counts_csv

dataset = load_counts_csv(
    "counts.csv",
    target_column="celltype",
)
```

### Use the internal surface-gene database

The package includes an internal `db_surface.csv`. To keep only surface genes:

```python
dataset = load_counts_csv(
    "counts.csv",
    target_column="celltype",
    gene_filter="surface",
)
```

By default, this keeps only:

Removes proteins whose names contain `Collagen`.

To keep all categories from the internal database:

```python
dataset = load_counts_csv(
    "counts.csv",
    target_column="celltype",
    gene_filter="surface",
    allowed_surface_categories=None,
)
```

### Use a custom gene list

```python
genes_interesse = ["PTPRC", "MRC1", "ABCA1", "PLXND1"]

dataset = load_counts_csv(
    "counts.csv",
    target_column="celltype",
    gene_filter="custom",
    genes=genes_interesse,
)
```

### Keep all genes

```python
dataset = load_counts_csv(
    "counts.csv",
    target_column="celltype",
    gene_filter="all",
)
```

### Encode labels

If your annotation column has labels such as `Macrof 1`, `Macrof2`, `Macrof3`, you can encode them as `0`, `1`, `2`:

```python
dataset = load_counts_csv(
    "counts.csv",
    target_column="celltype",
    gene_filter="surface",
    encode_labels=True,
)
```

The mapping is stored in:

```python
dataset.attrs["label_mapping"]
```

and, after training, saved to:

```text
cluster_outputs/label_mapping.csv
```

## Training One-vs-Rest Models

```python
from surface_biomarkers import TrainingConfig, train_one_vs_rest

config = TrainingConfig(
    resampling="SMOTE",
    cv=5,
    n_jobs=-1,
    verbose=1,
    grid_verbose=0,
    train_random_forest=True,
    train_lightgbm=True,
)

model_results = train_one_vs_rest(
    dataset,
    target_column="celltype",
    output_dir="cluster_outputs",
    config=config,
)

model_results
```

Verbosity options:

```python
verbose=0       # silent
verbose=1       # progress by cluster/model
verbose=2       # extra details
grid_verbose=2  # internal GridSearchCV messages
```

For each cluster, the package saves:

```text
cluster_outputs/0/X_train.csv
cluster_outputs/0/X_test.csv
cluster_outputs/0/y_train.csv
cluster_outputs/0/y_test.csv
cluster_outputs/0/scaler.pkl
cluster_outputs/0/lgb_model_0.pkl
cluster_outputs/0/rf_model_0.pkl
cluster_outputs/0/metrics_LGBM.csv
cluster_outputs/0/metrics_RF.csv
```

The global model summary is saved as:

```text
cluster_outputs/model_metrics.csv
```

## Per-Cluster Diagnostics

Use `analyze_clusters()` to generate confusion matrices, ROC curves, SHAP summaries, metric tables and a combined panel.

```python
from surface_biomarkers.plots import analyze_clusters

clusters = ["0", "1", "2", "3", "4", "5"]

metrics = analyze_clusters(
    base_dir="cluster_outputs",
    clusters=clusters,
    model_name=["LGBM", "RF"],
    max_display=10,
    panel=True,
    show_panel=True,
)

metrics
```

For each cluster/model, this saves:

```text
cluster_outputs/0/confusion_LGBM.svg
cluster_outputs/0/roc_LGBM.svg
cluster_outputs/0/shap_summary_LGBM.svg
cluster_outputs/0/panel_LGBM.svg
cluster_outputs/0/metrics_LGBM.csv
```

The panel contains:

```text
confusion matrix + ROC curve + SHAP summary
```

In Jupyter, only the panel is displayed by default. The individual plots are saved but not printed inline.

Use:

```python
show_panel=False
```

to save everything without displaying figures.

## Custom Cluster Names In Plots

If the dataset already uses numeric clusters but you want biological names in plots:

```python
cluster_names = {
    "0": "Mac1",
    "1": "Mac2",
    "2": "Mac3",
    "3": "Mac4",
    "4": "Mac5",
    "5": "Mac6",
}

metrics = analyze_clusters(
    base_dir="cluster_outputs",
    clusters=clusters,
    model_name=["LGBM", "RF"],
    cluster_names=cluster_names,
    show_panel=True,
)
```

If `encode_labels=True` was used during loading, the package automatically uses `cluster_outputs/label_mapping.csv`.

## Discover Compact Marker Signatures

```python
from surface_biomarkers import DiscoveryConfig, discover_signatures

discovery_config = DiscoveryConfig(
    model_name="LGBM",
    positive_genes=1,
    negative_genes=1,
    min_expression_prop=0.10,
    shap_corr_positive=0.15,
    shap_corr_negative=-0.10,
    top_n_heatmap=10,
)

signatures, gene_ranking, gene_exclusivity = discover_signatures(
    base_dir="cluster_outputs",
    clusters=clusters,
    config=discovery_config,
)

signatures
```

Example output:

```python
{
    "Cluster 0": ["PLXND1", "TLR1"],
    "Cluster 1": ["ABCA1"],
    "Cluster 2": ["CXADR", "FN1"],
}
```

Inspect the full ranking:

```python
gene_ranking.head(20)
```

Inspect the most exclusive genes:

```python
gene_exclusivity.head(20)
```

## SHAP Heatmaps

```python
from surface_biomarkers.signatures import (
    mean_abs_shap_matrix,
    combined_shap_expression_matrix,
)
from surface_biomarkers.plots import plot_shap_heatmap

shap_matrix = mean_abs_shap_matrix("cluster_outputs", clusters, discovery_config)
combined_matrix = combined_shap_expression_matrix("cluster_outputs", clusters, discovery_config)
```

Mean absolute SHAP heatmap:

```python
plot_shap_heatmap(
    shap_matrix,
    "cluster_outputs/shap_top_genes_heatmap.svg",
    title="Top genes by mean absolute SHAP",
    cluster_names=cluster_names,
    show=True,
)
```

SHAP x expression heatmap:

```python
plot_shap_heatmap(
    combined_matrix,
    "cluster_outputs/shap_expression_heatmap.svg",
    title="SHAP x expression score",
    cluster_names=cluster_names,
    show=True,
)
```


## Independent Accuracy Of Small Panels

`evaluate_signatures()` follows the original notebook logic:

1. read each `X_test.csv`;
2. assign all rows in that file to the cluster represented by its folder;
3. sum the expression of signature genes;
4. compute independent binary accuracy using `score > threshold`.

```python
from surface_biomarkers import evaluate_signatures

accuracy = evaluate_signatures(
    base_dir="cluster_outputs",
    clusters=clusters,
    signatures=signatures,
    threshold=0.5,
    cluster_names=cluster_names,
)

accuracy
```

Manual signatures also work:

```python
manual_signatures = {
    "Cluster 0": ["PLXND1", "TLR1"],
    "Cluster 1": ["ABCA1"],
    "Cluster 2": ["CXADR", "FN1"],
    "Cluster 3": ["LPL", "SEL1L3"],
    "Cluster 4": ["EMP1", "CD86"],
    "Cluster 5": ["IFNGR1", "ITGA4"],
}

accuracy = evaluate_signatures(
    base_dir="cluster_outputs",
    clusters=clusters,
    signatures=manual_signatures,
    threshold=0.5,
    cluster_names=cluster_names,
)
```

Plot the accuracy:

```python
from surface_biomarkers.plots import plot_signature_accuracy

plot_signature_accuracy(
    accuracy,
    "cluster_outputs/signature_accuracy.svg",
    cluster_names=cluster_names,
    show=True,
)
```

## Command-Line Usage

Train:

```bash
surface-biomarkers train counts.csv --target celltype --output cluster_outputs --resampling SMOTE
```

Analyze models:

```bash
surface-biomarkers analyze cluster_outputs --clusters 0 1 2 3 4 5 --model LGBM RF
```

Discover signatures:

```bash
surface-biomarkers discover cluster_outputs --clusters 0 1 2 3 4 5 --positive-genes 1 --negative-genes 1
```

The `discover` command saves:

```text
gene_ranking_by_cluster.csv
gene_exclusivity.csv
signature_specificity.csv
signature_accuracy.csv
shap_top_genes_heatmap.svg
shap_expression_heatmap.svg
signature_accuracy.svg
```

## Project Structure

```text
src/surface_biomarkers/
  data.py          # loading, gene filtering, label encoding
  training.py      # one-vs-rest model training
  signatures.py    # SHAP ranking, signatures, exclusivity, accuracy
  plots.py         # ROC, confusion matrix, SHAP, heatmaps, panels
  resources/
    db_surface.csv # internal surface-gene database from https://wlab.ethz.ch/surfaceome/
```

## Notes

- `gene_filter="surface"` uses the internal `db_surface.csv`.
- `evaluate_signatures()` intentionally reproduces the original notebook logic and does not read `y_test.csv`.
- `analyze_clusters()` saves individual plots and displays only the combined panel in Jupyter.
- Use all clusters together in `evaluate_signatures()` when reproducing the original independent-accuracy analysis.
