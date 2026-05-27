import pandas as pd
import matplotlib
from sklearn.ensemble import RandomForestClassifier

matplotlib.use("Agg")

from surface_biomarkers.data import filter_gene_list, filter_surface_genes, make_one_vs_rest, sanitize_name
from surface_biomarkers.panel_design import design_antibody_panel
from surface_biomarkers.plots import plot_panel_benchmark
from surface_biomarkers.signatures import calculate_exclusivity, rank_cluster_genes
from surface_biomarkers.training import _build_training_pipeline, _prefix_param_grid, get_resampler
from surface_biomarkers.validation import prioritize_signature_validation


def test_make_one_vs_rest():
    df = pd.DataFrame({"GeneA": [1, 2, 3], "celltype": ["A", "B", "A"]})
    out = make_one_vs_rest(df, "A")
    assert list(out["Class"]) == [1, 0, 1]


def test_sanitize_name():
    assert sanitize_name("Mac/1 alpha") == "Mac_1_alpha"


def test_filter_surface_genes_keeps_target_column():
    counts = pd.DataFrame({
        "PTPRC": [1, 0],
        "NOT_SURFACE": [5, 6],
        "MRC1": [2, 3],
        "celltype": ["Mac1", "Mac2"],
    })
    db_surface = pd.DataFrame({
        "ENTREZ gene symbol": ["PTPRC", "MRC1", "COL1A1"],
        "CSPA category": ["1 - high confidence", "2 - putative", "1 - high confidence"],
        "UP_Protein_name": ["CD45", "CD206", "Collagen alpha"],
    })
    out = filter_surface_genes(counts, db_surface)
    assert list(out.columns) == ["PTPRC", "MRC1", "celltype"]


def test_filter_gene_list_keeps_requested_genes():
    counts = pd.DataFrame({
        "PTPRC": [1, 0],
        "NOT_SURFACE": [5, 6],
        "MRC1": [2, 3],
        "celltype": ["Mac1", "Mac2"],
    })
    out = filter_gene_list(counts, ["MRC1"], keep_columns=("celltype",))
    assert list(out.columns) == ["MRC1", "celltype"]


def test_exclusivity():
    matrix = pd.DataFrame({"Cluster 0": [1.0, 0.2], "Cluster 1": [0.3, 0.8]}, index=["A", "B"])
    out = calculate_exclusivity(matrix)
    assert out.iloc[0]["gene"] == "A"
    assert round(out.iloc[0]["exclusivity"], 2) == 0.7


def test_training_pipeline_keeps_resampling_inside_cv():
    pipeline = _build_training_pipeline(
        RandomForestClassifier(random_state=1),
        get_resampler("RandomOverSampler", random_state=1),
    )

    assert list(pipeline.named_steps) == ["resampler", "scaler", "model"]
    assert _prefix_param_grid({"n_estimators": [10]}) == {"model__n_estimators": [10]}


def test_rank_cluster_genes_uses_target_cells_for_primary_shap_ranking():
    X_raw = pd.DataFrame({
        "TARGET_MARKER": [5.0, 6.0, 0.0, 0.0],
        "NEGATIVE_MARKER": [1.0, 1.0, 9.0, 10.0],
    })
    X_scaled = X_raw.copy()
    y_test = pd.Series([1, 1, 0, 0])
    shap_values = pd.DataFrame({
        "TARGET_MARKER": [3.0, 3.0, 0.1, 0.1],
        "NEGATIVE_MARKER": [0.2, 0.2, 9.0, 9.0],
    }).to_numpy()

    ranking = rank_cluster_genes(X_raw, X_scaled, y_test, shap_values)

    assert ranking.iloc[0]["gene"] == "TARGET_MARKER"
    assert ranking.iloc[0]["mean_abs_shap"] == 3.0


def test_prioritize_signature_validation_scores_signature_gene():
    ranking = pd.DataFrame({
        "cluster": ["Cluster 0", "Cluster 0"],
        "gene": ["PTPRC", "OTHER"],
        "mean_abs_shap": [2.0, 1.0],
        "mean_shap": [0.5, -0.1],
        "mean_expression_target": [10.0, 2.0],
        "expression_prop_target": [0.9, 0.2],
        "expression_shap_corr": [0.8, -0.3],
        "technical_corr": [0.0, 0.0],
    })
    surface_db = pd.DataFrame({
        "ENTREZ gene symbol": ["PTPRC"],
        "CSPA category": ["1 - high confidence"],
        "UniProt Cell surface": ["yes"],
        "# predicted TM domains": [1],
        "GPI": [0],
        "Phobius TM predicted yes/no": ["yes"],
    })
    protein_atlas = pd.DataFrame({
        "Gene": ["PTPRC"],
        "Antibody": ["HPA000001"],
        "Antibody RRID": ["AB_123"],
        "Reliability (IF)": ["Supported"],
        "Reliability (IH)": ["Approved"],
        "Subcellular location": ["Plasma membrane"],
        "Subcellular main location": ["Plasma membrane"],
    })

    out = prioritize_signature_validation(
        ranking,
        {"Cluster 0": ["PTPRC"]},
        surface_db,
        protein_atlas=protein_atlas,
    )

    assert list(out["gene"]) == ["PTPRC"]
    assert bool(out.iloc[0]["is_signature"])
    assert out.iloc[0]["surface_score"] == 8
    assert out.iloc[0]["snRNAseq_validation_class"] == "antibody_supported_snRNAseq_candidate"
    assert out.iloc[0]["FACS_compatibility"] == "high_priority_candidate"
    assert out.iloc[0]["Flow_cytometry_compatibility"] == "high_priority_candidate"
    assert out.iloc[0]["IF_compatibility"] == "compatible_candidate"


def test_prioritize_signature_validation_exposes_new_score_components():
    ranking = pd.DataFrame({
        "cluster": ["Cluster 0"],
        "gene": ["PTPRC"],
        "mean_abs_shap": [2.0],
        "mean_shap": [0.5],
        "mean_expression_target": [10.0],
        "expression_prop_target": [0.9],
        "expression_shap_corr": [0.8],
        "technical_corr": [0.0],
    })
    surface_db = pd.DataFrame({
        "ENTREZ gene symbol": ["PTPRC"],
        "CSPA category": ["1 - high confidence"],
        "UniProt Cell surface": ["yes"],
        "# predicted TM domains": [1],
        "GPI": [0],
        "Phobius TM predicted yes/no": ["yes"],
    })

    out = prioritize_signature_validation(
        ranking,
        {"Cluster 0": ["PTPRC"]},
        surface_db,
        protein_atlas=None,
    )

    for column in [
        "transcript_score",
        "specificity_score",
        "prevalence_score",
        "effect_size_score",
        "donor_robustness_score",
        "localization_score",
        "antibody_availability_score",
        "application_compatibility_score",
        "penalty_score",
        "final_marker_score",
    ]:
        assert column in out.columns


def test_final_marker_score_penalizes_rna_only_marker_without_protein_support():
    ranking = pd.DataFrame({
        "cluster": ["Cluster 0", "Cluster 0"],
        "gene": ["GOOD_SURFACE", "RNA_ONLY"],
        "mean_abs_shap": [8.0, 10.0],
        "mean_shap": [1.0, 1.5],
        "mean_expression_target": [5.0, 6.0],
        "mean_expression_rest": [0.5, 5.5],
        "expression_prop_target": [0.8, 0.9],
        "expression_prop_rest": [0.1, 0.8],
        "expression_shap_corr": [0.7, 0.9],
        "technical_corr": [0.0, 0.0],
        "donor_min_expression_prop": [0.7, 0.7],
        "donor_expression_prop_cv": [0.2, 0.2],
    })
    surface_db = pd.DataFrame({
        "ENTREZ gene symbol": ["GOOD_SURFACE"],
        "CSPA category": ["1 - high confidence"],
        "UniProt Cell surface": ["yes"],
        "# predicted TM domains": [1],
        "GPI": [0],
        "Phobius TM predicted yes/no": ["yes"],
    })
    protein_atlas = pd.DataFrame({
        "Gene": ["GOOD_SURFACE"],
        "Antibody": ["HPA000001"],
        "Antibody RRID": ["AB_123"],
        "Reliability (IF)": ["Supported"],
        "Reliability (IH)": ["Supported"],
        "Subcellular location": ["Plasma membrane"],
        "Subcellular main location": ["Plasma membrane"],
    })

    out = prioritize_signature_validation(
        ranking,
        {"Cluster 0": ["GOOD_SURFACE", "RNA_ONLY"]},
        surface_db,
        protein_atlas=protein_atlas,
    ).set_index("gene")

    assert out.loc["RNA_ONLY", "transcript_score"] > out.loc["GOOD_SURFACE", "transcript_score"]
    assert out.loc["RNA_ONLY", "localization_score"] == 0
    assert out.loc["RNA_ONLY", "antibody_availability_score"] == 0
    assert out.loc["RNA_ONLY", "penalty_score"] > out.loc["GOOD_SURFACE", "penalty_score"]
    assert out.loc["GOOD_SURFACE", "final_marker_score"] > out.loc["RNA_ONLY", "final_marker_score"]


def test_design_antibody_panel_selects_positive_negative_and_backups():
    markers = pd.DataFrame({
        "cluster": ["Cluster 0"] * 5,
        "gene": ["CD_A", "CD_B", "EXCLUDE_C", "BACKUP_D", "RNA_ONLY"],
        "is_signature": [True, True, True, True, True],
        "mean_shap": [1.2, 1.0, -0.8, 0.7, 2.0],
        "expression_shap_corr": [0.8, 0.7, -0.6, 0.6, 0.9],
        "final_marker_score": [9.0, 8.5, 7.5, 7.0, 9.5],
        "specificity_score": [8.0, 7.0, 8.0, 6.0, 9.0],
        "prevalence_score": [8.0, 8.0, 7.0, 6.0, 9.0],
        "localization_score": [9.0, 8.0, 8.0, 7.0, 0.0],
        "antibody_availability_score": [8.0, 8.0, 7.0, 6.0, 0.0],
        "application_compatibility_score": [8.5, 8.0, 7.5, 6.5, 2.0],
        "penalty_score": [0.0, 0.0, 0.0, 0.0, 5.0],
    })

    panel = design_antibody_panel(markers, application="FACS")
    row = panel.iloc[0]

    assert row["positive_markers"] == "CD_A; CD_B"
    assert row["negative_exclusion_markers"] == "EXCLUDE_C"
    assert row["backup_markers"] == "BACKUP_D"
    assert "RNA_ONLY" not in row["all_candidate_markers"]
    assert bool(row["constraints_passed"])


def test_design_antibody_panel_applies_redundancy_penalty():
    markers = pd.DataFrame({
        "cluster": ["Cluster 0"] * 3,
        "gene": ["CD_A", "CD_B", "CD_C"],
        "is_signature": [True, True, True],
        "mean_shap": [1.0, 0.9, 0.8],
        "expression_shap_corr": [0.8, 0.8, 0.8],
        "final_marker_score": [9.0, 8.8, 8.0],
        "specificity_score": [8.0, 8.0, 8.0],
        "prevalence_score": [8.0, 8.0, 8.0],
        "localization_score": [9.0, 9.0, 9.0],
        "antibody_availability_score": [8.0, 8.0, 8.0],
        "application_compatibility_score": [8.0, 8.0, 8.0],
        "penalty_score": [0.0, 0.0, 0.0],
    })
    redundancy = pd.DataFrame(
        [
            [1.0, 0.95, 0.10],
            [0.95, 1.0, 0.10],
            [0.10, 0.10, 1.0],
        ],
        index=["CD_A", "CD_B", "CD_C"],
        columns=["CD_A", "CD_B", "CD_C"],
    )

    panel = design_antibody_panel(markers, application="CyTOF", redundancy_matrix=redundancy)
    row = panel.iloc[0]

    assert row["positive_markers"] == "CD_A; CD_C"
    assert "redundant_markers_skipped:CD_B" in row["warnings"]


def test_plot_panel_benchmark_writes_svg(tmp_path):
    panel = pd.DataFrame({
        "cluster": ["Cluster 0", "Cluster 1"],
        "application": ["FACS", "FACS"],
        "positive_markers": ["CD_A; CD_B", "CD_C"],
        "negative_exclusion_markers": ["EXCLUDE_A", ""],
        "backup_markers": ["BACKUP_A", "BACKUP_B"],
        "minimal_panel": ["CD_A; CD_B; EXCLUDE_A", "CD_C"],
        "panel_size": [3, 1],
        "backup_count": [1, 1],
        "redundancy_penalty": [0.1, 0.0],
        "panel_score": [12.5, 8.0],
        "constraints_passed": [True, False],
        "warnings": ["", "insufficient_positive_markers"],
        "all_candidate_markers": ["CD_A; CD_B; EXCLUDE_A; BACKUP_A", "CD_C; BACKUP_B"],
    })

    output = plot_panel_benchmark(panel, tmp_path / "panel_benchmark.svg", show=False)

    assert output.exists()
    assert output.stat().st_size > 0
