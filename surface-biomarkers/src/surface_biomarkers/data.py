import re
from pathlib import Path
from typing import List, Optional, Tuple, Union

import pandas as pd


def sanitize_name(value: object) -> str:
    text = str(value).replace("/", "_").replace("\\", "_")
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    return text.strip("_") or "cluster"


def load_counts_csv(
    path: Union[str, Path],
    target_column: str = "celltype",
    gene_filter: str = "all",
    genes: Optional[Union[List[str], Tuple[str, ...], str, Path]] = None,
    encode_labels: bool = False,
    surface_gene_column: str = "ENTREZ gene symbol",
    surface_category_column: str = "CSPA category",
    allowed_surface_categories: Optional[Tuple[str, ...]] = ("1 - high confidence", "2 - putative"),
    protein_name_column: str = "UP_Protein_name",
    exclude_protein_pattern: Optional[str] = "Collagen",
) -> pd.DataFrame:
    """Load counts and optionally filter genes.

    The input count matrix is expected to have cells in rows, genes in columns,
    and one annotation column such as ``celltype``.

    ``gene_filter="all"`` keeps every gene.
    ``gene_filter="surface"`` uses the internal ``db_surface.csv`` packaged with
    this library.
    ``gene_filter="custom"`` keeps the genes supplied in ``genes``.
    If ``encode_labels=True``, labels in ``target_column`` are converted to
    ``"0"``, ``"1"``, ... and the original labels are stored in
    ``df.attrs["label_mapping"]``.
    """
    df = pd.read_csv(path, index_col=0)
    if target_column not in df.columns:
        raise ValueError(f"Column '{target_column}' was not found in {path}.")
    if df[target_column].isna().any():
        raise ValueError(f"Column '{target_column}' contains missing labels.")
    if encode_labels:
        df = encode_target_labels(df, target_column=target_column)

    if gene_filter == "all":
        return df
    if gene_filter == "surface":
        surface_genes = load_internal_surface_db()
        df = filter_surface_genes(
            df,
            surface_genes,
            gene_column=surface_gene_column,
            category_column=surface_category_column,
            allowed_categories=allowed_surface_categories,
            protein_name_column=protein_name_column,
            exclude_protein_pattern=exclude_protein_pattern,
            keep_columns=(target_column,),
        )
    elif gene_filter == "custom":
        if genes is None:
            raise ValueError("gene_filter='custom' requires genes.")
        df = filter_gene_list(df, genes, keep_columns=(target_column,))
    else:
        raise ValueError("gene_filter must be one of: 'all', 'surface', 'custom'.")
    return df


def encode_target_labels(df: pd.DataFrame, target_column: str = "celltype") -> pd.DataFrame:
    """Encode target labels as strings "0", "1", ... and store the mapping."""
    out = df.copy()
    labels = sorted(out[target_column].dropna().unique())
    mapping = {original: str(i) for i, original in enumerate(labels)}
    out[target_column] = out[target_column].map(mapping)
    out.attrs["label_mapping"] = pd.DataFrame({
        "encoded_label": [mapping[label] for label in labels],
        "original_label": labels,
    })
    return out


def load_internal_surface_db() -> pd.DataFrame:
    """Load the packaged surface-gene database."""
    resource = Path(__file__).resolve().parent / "resources" / "db_surface.csv"
    return pd.read_csv(resource)


def _load_gene_list(genes: Union[List[str], Tuple[str, ...], str, Path]) -> List[str]:
    if isinstance(genes, (list, tuple)):
        return [str(g) for g in genes]
    path = Path(genes)
    if path.suffix.lower() == ".csv":
        df = pd.read_csv(path)
        if df.shape[1] == 0:
            return []
        return df.iloc[:, 0].dropna().astype(str).tolist()
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def make_one_vs_rest(
    df: pd.DataFrame,
    cluster: object,
    target_column: str = "celltype",
    class_column: str = "Class",
) -> pd.DataFrame:
    out = df.copy()
    out[class_column] = (out[target_column] == cluster).astype(int)
    return out.drop(columns=[target_column])


def filter_surface_genes(
    counts: pd.DataFrame,
    surface_genes: pd.DataFrame,
    gene_column: str = "ENTREZ gene symbol",
    category_column: str = "CSPA category",
    allowed_categories: Optional[Tuple[str, ...]] = ("1 - high confidence", "2 - putative"),
    protein_name_column: str = "UP_Protein_name",
    exclude_protein_pattern: Optional[str] = "Collagen",
    keep_columns: Tuple[str, ...] = ("celltype",),
) -> pd.DataFrame:
    genes = surface_genes.copy()
    if gene_column not in genes.columns:
        raise ValueError(f"Surface gene column '{gene_column}' was not found in db_surface.")
    if allowed_categories is not None and category_column in genes:
        genes = genes[genes[category_column].isin(allowed_categories)]
    if exclude_protein_pattern and protein_name_column in genes:
        mask = ~genes[protein_name_column].fillna("").str.contains(
            exclude_protein_pattern, case=False, regex=True
        )
        genes = genes[mask]

    keep_genes = set(genes[gene_column].dropna().astype(str))
    selected = [col for col in counts.columns if col in keep_genes]
    selected.extend([col for col in keep_columns if col in counts.columns])
    if not selected:
        raise ValueError("No surface genes from db_surface were found in the count matrix.")
    return counts.loc[:, selected].copy()


def filter_gene_list(
    counts: pd.DataFrame,
    genes: Union[List[str], Tuple[str, ...], str, Path],
    keep_columns: Tuple[str, ...] = ("celltype",),
) -> pd.DataFrame:
    keep_genes = set(_load_gene_list(genes))
    selected = [col for col in counts.columns if col in keep_genes]
    selected.extend([col for col in keep_columns if col in counts.columns])
    if not selected:
        raise ValueError("None of the requested genes were found in the count matrix.")
    return counts.loc[:, selected].copy()
