from pathlib import Path
from typing import Dict, List, Optional, Union

import pandas as pd


SURFACE_GENE_COLUMN = "ENTREZ gene symbol"
HPA_GENE_COLUMN = "Gene"


def internal_resource_path(filename: str) -> Path:
    return Path(__file__).resolve().parent / "resources" / filename


def load_internal_surface_db() -> pd.DataFrame:
    return pd.read_csv(internal_resource_path("db_surface.csv"))


def load_internal_protein_atlas() -> pd.DataFrame:
    return pd.read_csv(internal_resource_path("proteinatlas.tsv"), sep="\t")


def _read_table(table: Union[str, Path, pd.DataFrame], sep: Optional[str] = None) -> pd.DataFrame:
    if isinstance(table, pd.DataFrame):
        return table.copy()
    path = Path(table)
    if sep is None:
        sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep)


def _minmax(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").fillna(0.0)
    min_value = values.min()
    max_value = values.max()
    if max_value == min_value:
        return pd.Series(0.0, index=series.index)
    return (values - min_value) / (max_value - min_value)


def _score_0_10(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").fillna(0.0).clip(lower=0.0, upper=10.0)


def _optional_numeric(df: pd.DataFrame, column: str, default: float = 0.0) -> pd.Series:
    if column not in df.columns:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[column], errors="coerce").fillna(default)


def _truthy(value) -> int:
    if pd.isna(value):
        return 0
    if isinstance(value, str):
        text = value.strip().lower()
        if text in {"yes", "true", "t", "1", "y", "positive"}:
            return 1
        if text in {"no", "false", "f", "0", "n", "negative", ""}:
            return 0
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(numeric):
        return 0
    return int(numeric > 0)


def _reliability_score(value) -> int:
    if pd.isna(value):
        return 0
    text = str(value).strip().lower()
    scores = {
        "enhanced": 3,
        "supported": 2,
        "approved": 1,
        "uncertain": -1,
    }
    return scores.get(text, 0)


def _join_unique(values: pd.Series) -> str:
    clean = [str(v) for v in values.dropna().astype(str) if str(v).strip()]
    if not clean:
        return ""
    return "; ".join(sorted(set(clean)))


def _max_reliability(values: pd.Series) -> str:
    if values.dropna().empty:
        return ""
    order = {"Enhanced": 3, "Supported": 2, "Approved": 1, "Uncertain": -1}
    best_value = ""
    best_score = -999
    for value in values.dropna().astype(str):
        score = order.get(value.strip(), _reliability_score(value))
        if score > best_score:
            best_score = score
            best_value = value
    return best_value


def _signature_pairs(signatures: Dict[str, List[str]]) -> pd.DataFrame:
    rows = []
    for cluster, genes in signatures.items():
        for gene in genes:
            rows.append({"cluster": cluster, "gene": gene, "is_signature": True})
    if not rows:
        return pd.DataFrame(columns=["cluster", "gene", "is_signature"])
    return pd.DataFrame(rows).drop_duplicates(["cluster", "gene"])


def _prepare_surface_db(surface_db: Union[str, Path, pd.DataFrame]) -> pd.DataFrame:
    surface = (
        load_internal_surface_db()
        if isinstance(surface_db, str) and surface_db == "internal"
        else _read_table(surface_db)
    )
    if SURFACE_GENE_COLUMN not in surface.columns:
        raise ValueError(f"surface_db must contain '{SURFACE_GENE_COLUMN}'.")

    for column in [
        "CSPA category",
        "UniProt Cell surface",
        "# predicted TM domains",
        "GPI",
        "Phobius TM predicted yes/no",
    ]:
        if column not in surface.columns:
            surface[column] = pd.NA

    surface["gene"] = surface[SURFACE_GENE_COLUMN].astype(str)
    surface["_uniprot_cell_surface"] = surface["UniProt Cell surface"].map(_truthy)
    surface["_tm_domain"] = pd.to_numeric(surface["# predicted TM domains"], errors="coerce").fillna(0).gt(0).astype(int)
    surface["_gpi_anchor"] = pd.to_numeric(surface["GPI"], errors="coerce").fillna(0).gt(0).astype(int)
    surface["_phobius_tm"] = surface["Phobius TM predicted yes/no"].map(_truthy)

    return surface.groupby("gene", as_index=False).agg({
        "CSPA category": _join_unique,
        "UniProt Cell surface": _join_unique,
        "# predicted TM domains": "max",
        "GPI": "max",
        "_uniprot_cell_surface": "max",
        "_tm_domain": "max",
        "_gpi_anchor": "max",
        "_phobius_tm": "max",
    })


def _prepare_protein_atlas(protein_atlas: Optional[Union[str, Path, pd.DataFrame]]) -> pd.DataFrame:
    columns = [
        "gene",
        "possible_antibodies",
        "possible_antibody_rrids",
        "Reliability (IF)",
        "Reliability (IH)",
        "Subcellular location",
        "Subcellular main location",
        "_has_hpa_antibody",
        "_IF_reliability_score",
        "_IHC_reliability_score",
    ]
    if protein_atlas is None:
        return pd.DataFrame(columns=columns)

    atlas = (
        load_internal_protein_atlas()
        if isinstance(protein_atlas, str) and protein_atlas == "internal"
        else _read_table(protein_atlas)
    )
    if HPA_GENE_COLUMN not in atlas.columns:
        raise ValueError(f"protein_atlas must contain '{HPA_GENE_COLUMN}'.")

    for column in [
        "Antibody",
        "Antibody RRID",
        "Reliability (IF)",
        "Reliability (IH)",
        "Subcellular location",
        "Subcellular main location",
    ]:
        if column not in atlas.columns:
            atlas[column] = pd.NA

    atlas["gene"] = atlas[HPA_GENE_COLUMN].astype(str)
    atlas["_has_hpa_antibody"] = atlas["Antibody"].notna().astype(int)
    atlas["_IF_reliability_score"] = atlas["Reliability (IF)"].map(_reliability_score)
    atlas["_IHC_reliability_score"] = atlas["Reliability (IH)"].map(_reliability_score)

    prepared = atlas.groupby("gene", as_index=False).agg({
        "Antibody": _join_unique,
        "Antibody RRID": _join_unique,
        "Reliability (IF)": _max_reliability,
        "Reliability (IH)": _max_reliability,
        "Subcellular location": _join_unique,
        "Subcellular main location": _join_unique,
        "_has_hpa_antibody": "max",
        "_IF_reliability_score": "max",
        "_IHC_reliability_score": "max",
    })
    return prepared.rename(columns={
        "Antibody": "possible_antibodies",
        "Antibody RRID": "possible_antibody_rrids",
    })


def classify_snRNAseq_candidate(row):
    if row["surface_score"] >= 5 and (
        row["_IF_reliability_score"] >= 1 or row["_IHC_reliability_score"] >= 1
    ):
        return "antibody_supported_snRNAseq_candidate"

    if row["surface_score"] >= 5:
        return "surface_supported_snRNAseq_candidate"

    if row["model_marker_score"] > 0:
        return "candidate_only_snRNAseq"

    return "not_recommended"


def classify_facs(row):
    if row["surface_score"] >= 5:
        return "high_priority_candidate"

    if row["surface_score"] >= 3:
        return "candidate_only_snRNAseq"

    return "not_recommended"


def classify_if(row):
    if row["_IF_reliability_score"] >= 1:
        return "compatible_candidate"

    return "needs_validation"


def classify_ihc(row):
    if row["_IHC_reliability_score"] >= 1:
        return "compatible_candidate"

    return "needs_validation"


def add_marker_score_components(output: pd.DataFrame) -> pd.DataFrame:
    """Add separated marker-prioritization scores.

    The components keep transcriptomic evidence separate from protein
    localization, antibody evidence and panel/application suitability. This
    prevents strong RNA markers from being ranked as strong antibody markers
    when the protein-level evidence is weak.
    """
    output = output.copy()
    expression_prop = _optional_numeric(output, "expression_prop_target")
    mean_expression = _optional_numeric(output, "mean_expression_target")
    mean_abs_shap_norm = _optional_numeric(output, "mean_abs_shap_norm")
    mean_expression_norm = _optional_numeric(output, "mean_expression_norm")
    expression_corr = _optional_numeric(output, "expression_shap_corr")

    output["transcript_score"] = _score_0_10(
        10 * (0.70 * mean_abs_shap_norm + 0.30 * mean_expression_norm)
    )
    output["prevalence_score"] = _score_0_10(10 * expression_prop.clip(lower=0.0, upper=1.0))

    mean_expression_rest = _optional_numeric(output, "mean_expression_rest", default=-1.0)
    has_expression_rest = "mean_expression_rest" in output.columns
    if has_expression_rest:
        expression_ratio = mean_expression / (mean_expression + mean_expression_rest.clip(lower=0.0) + 1e-9)
        expression_delta_norm = _minmax(mean_expression - mean_expression_rest)
        output["effect_size_score"] = _score_0_10(
            10 * (0.50 * expression_ratio + 0.50 * expression_delta_norm)
        )
    else:
        output["effect_size_score"] = _score_0_10(10 * mean_expression_norm)

    expression_prop_rest = _optional_numeric(output, "expression_prop_rest", default=-1.0)
    has_prop_rest = "expression_prop_rest" in output.columns
    if has_prop_rest:
        prop_delta = (expression_prop - expression_prop_rest).clip(lower=0.0, upper=1.0)
        if has_expression_rest:
            expression_ratio = mean_expression / (mean_expression + mean_expression_rest.clip(lower=0.0) + 1e-9)
        else:
            expression_ratio = mean_expression_norm
        output["specificity_score"] = _score_0_10(10 * (0.60 * prop_delta + 0.40 * expression_ratio))
    else:
        output["specificity_score"] = _score_0_10(10 * expression_corr.clip(lower=0.0, upper=1.0))

    if "donor_robustness_score" in output.columns:
        output["donor_robustness_score"] = _score_0_10(output["donor_robustness_score"])
    elif {"donor_min_expression_prop", "donor_expression_prop_cv"}.issubset(output.columns):
        donor_min = _optional_numeric(output, "donor_min_expression_prop").clip(lower=0.0, upper=1.0)
        donor_cv = _optional_numeric(output, "donor_expression_prop_cv", default=1.0).clip(lower=0.0, upper=1.0)
        output["donor_robustness_score"] = _score_0_10(10 * (0.70 * donor_min + 0.30 * (1 - donor_cv)))
    else:
        output["donor_robustness_score"] = 5.0

    output["localization_score"] = _score_0_10(10 * _optional_numeric(output, "surface_score") / 9.0)

    if_score = _optional_numeric(output, "_IF_reliability_score").clip(lower=0.0, upper=3.0) / 3.0
    ihc_score = _optional_numeric(output, "_IHC_reliability_score").clip(lower=0.0, upper=3.0) / 3.0
    has_antibody = _optional_numeric(output, "_has_hpa_antibody").clip(lower=0.0, upper=1.0)
    output["antibody_availability_score"] = _score_0_10(
        10 * (0.40 * has_antibody + 0.30 * if_score + 0.30 * ihc_score)
    )

    output["application_compatibility_score"] = _score_0_10(
        0.45 * output["localization_score"]
        + 0.30 * output["antibody_availability_score"]
        + 0.25 * output["prevalence_score"]
    )

    technical_corr = _optional_numeric(output, "technical_corr")
    penalty = pd.Series(0.0, index=output.index)
    penalty += ((output["transcript_score"] >= 7.0) & (output["localization_score"] < 4.0)).astype(float) * 3.0
    penalty += (output["antibody_availability_score"] == 0.0).astype(float) * 2.0
    penalty += (output["prevalence_score"] < 2.0).astype(float) * 2.0
    penalty += (technical_corr > 0.7).astype(float) * 2.0
    penalty += (output["donor_robustness_score"] < 4.0).astype(float) * 1.0
    output["penalty_score"] = _score_0_10(penalty)

    raw_final = (
        0.25 * output["transcript_score"]
        + 0.15 * output["specificity_score"]
        + 0.10 * output["prevalence_score"]
        + 0.10 * output["effect_size_score"]
        + 0.10 * output["donor_robustness_score"]
        + 0.15 * output["localization_score"]
        + 0.10 * output["antibody_availability_score"]
        + 0.05 * output["application_compatibility_score"]
        - output["penalty_score"]
    )
    output["final_marker_score"] = _score_0_10(raw_final)
    return output


def prioritize_signature_validation(
    ranking: pd.DataFrame,
    signatures: Dict[str, List[str]],
    surface_db: Union[str, Path, pd.DataFrame] = "internal",
    protein_atlas: Optional[Union[str, Path, pd.DataFrame]] = "internal",
    signature_only: bool = True,
) -> pd.DataFrame:
    """Prioritize SHAP-derived signatures for translational validation.

    The model was trained on single-nuclei RNA-seq, so surface assays such as
    FACS, CyTOF and CITE-seq are treated as candidates that still need wet-lab
    validation. The scores below combine model evidence, expression evidence,
    surface-protein evidence and antibody availability without claiming that
    snRNA-seq alone validates an assay.
    """
    required = [
        "cluster",
        "gene",
        "mean_abs_shap",
        "mean_shap",
        "mean_expression_target",
        "expression_prop_target",
        "expression_shap_corr",
    ]
    missing = [column for column in required if column not in ranking.columns]
    if missing:
        raise ValueError("ranking is missing required columns: " + ", ".join(missing))

    ranked = ranking.copy()
    ranked["gene"] = ranked["gene"].astype(str)
    ranked["cluster"] = ranked["cluster"].astype(str)

    ranked["mean_abs_shap_norm"] = _minmax(ranked["mean_abs_shap"])
    ranked["expression_prop_norm"] = _minmax(ranked["expression_prop_target"])
    ranked["mean_expression_norm"] = _minmax(ranked["mean_expression_target"])
    ranked["model_marker_score"] = (
        4 * ranked["mean_abs_shap_norm"]
        + 2 * ranked["expression_prop_norm"]
        + ranked["mean_expression_norm"]
    )

    signature_df = _signature_pairs(signatures)
    output = ranked.merge(signature_df, on=["cluster", "gene"], how="left")
    output["is_signature"] = output["is_signature"].fillna(False).astype(bool)
    if signature_only:
        output = output[output["is_signature"]].copy()

    surface = _prepare_surface_db(surface_db)
    output = output.merge(surface, on="gene", how="left")
    output["is_in_surface_db"] = output["CSPA category"].notna().astype(int)

    for column in ["_uniprot_cell_surface", "_tm_domain", "_gpi_anchor", "_phobius_tm"]:
        output[column] = output[column].fillna(0).astype(int)

    output["surface_score"] = (
        2 * output["is_in_surface_db"]
        + 3 * output["_uniprot_cell_surface"]
        + 2 * output["_tm_domain"]
        + output["_gpi_anchor"]
        + output["_phobius_tm"]
    )

    atlas = _prepare_protein_atlas(protein_atlas)
    output = output.merge(atlas, on="gene", how="left")
    output["_has_hpa_antibody"] = output["_has_hpa_antibody"].fillna(0).astype(int)
    output["_IF_reliability_score"] = output["_IF_reliability_score"].fillna(0).astype(int)
    output["_IHC_reliability_score"] = output["_IHC_reliability_score"].fillna(0).astype(int)

    output["facs_score"] = output["model_marker_score"] + output["surface_score"]
    output["cite_score"] = output["facs_score"]
    output["cytof_score"] = output["facs_score"]
    output["if_score"] = output["model_marker_score"] + 2 * output["_has_hpa_antibody"] + 2 * output["_IF_reliability_score"]
    output["ihc_score"] = output["model_marker_score"] + 2 * output["_has_hpa_antibody"] + 2 * output["_IHC_reliability_score"]
    output = add_marker_score_components(output)

    # snRNA-seq supports candidate prioritization, but protein-level assays
    # still require experimental validation before markers are considered validated.
    output["snRNAseq_validation_class"] = output.apply(
        classify_snRNAseq_candidate,
        axis=1
    )

    output["FACS_compatibility"] = output.apply(classify_facs, axis=1)
    output["Flow_cytometry_compatibility"] = output["FACS_compatibility"]
    output["CITEseq_compatibility"] = output["FACS_compatibility"]
    output["CyTOF_compatibility"] = output["FACS_compatibility"]

    output["IF_compatibility"] = output.apply(classify_if, axis=1)
    output["IHC_compatibility"] = output.apply(classify_ihc, axis=1)

    fill_columns = [
        "possible_antibodies",
        "possible_antibody_rrids",
        "Reliability (IF)",
        "Reliability (IH)",
        "Subcellular location",
        "Subcellular main location",
        "CSPA category",
        "UniProt Cell surface",
    ]
    for column in fill_columns:
        if column not in output.columns:
            output[column] = ""
        output[column] = output[column].fillna("")
    for column in ["# predicted TM domains", "GPI"]:
        if column not in output.columns:
            output[column] = 0
        output[column] = pd.to_numeric(output[column], errors="coerce").fillna(0)

    columns = [
        "cluster",
        "gene",
        "is_signature",
        "mean_shap",
        "mean_expression_target",
        "expression_prop_target",
        "expression_shap_corr",
        "surface_score",
        "model_marker_score",
        "facs_score",
        "cite_score",
        "cytof_score",
        "if_score",
        "ihc_score",
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
        "snRNAseq_validation_class",
        "FACS_compatibility",
        "Flow_cytometry_compatibility",
        "CITEseq_compatibility",
        "CyTOF_compatibility",
        "IF_compatibility",
        "IHC_compatibility",
        "possible_antibodies",
        "possible_antibody_rrids",
        "Reliability (IF)",
        "Reliability (IH)",
        "CSPA category",
        "UniProt Cell surface",
        "# predicted TM domains",
        "GPI",
    ]
    return output[columns].sort_values(
        ["is_signature", "final_marker_score", "facs_score", "if_score", "ihc_score"],
        ascending=[False, False, False, False, False],
    ).reset_index(drop=True)


def save_signature_validation_report(
    ranking: pd.DataFrame,
    signatures: Dict[str, List[str]],
    surface_db: Union[str, Path, pd.DataFrame] = "internal",
    protein_atlas: Optional[Union[str, Path, pd.DataFrame]] = "internal",
    output_path: Union[str, Path] = "signature_validation_priority.csv",
    signature_only: bool = True,
) -> pd.DataFrame:
    report = prioritize_signature_validation(
        ranking=ranking,
        signatures=signatures,
        surface_db=surface_db,
        protein_atlas=protein_atlas,
        signature_only=signature_only,
    )
    report.to_csv(output_path, index=False)
    return report
