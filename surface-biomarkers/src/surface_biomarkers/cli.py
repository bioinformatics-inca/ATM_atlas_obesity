import argparse
import json
from pathlib import Path
from typing import List, Optional

import pandas as pd

from surface_biomarkers.data import load_counts_csv
from surface_biomarkers.plots import analyze_clusters, plot_signature_accuracy, plot_shap_heatmap
from surface_biomarkers.signatures import (
    DiscoveryConfig,
    combined_shap_expression_matrix,
    discover_signatures,
    evaluate_signatures,
    mean_abs_shap_matrix,
)
from surface_biomarkers.training import TrainingConfig, train_one_vs_rest
from surface_biomarkers.validation import save_signature_validation_report


def _load_signatures(path: str):
    signature_path = Path(path)
    if signature_path.suffix.lower() == ".json":
        with signature_path.open("r", encoding="utf-8") as handle:
            return json.load(handle)

    df = pd.read_csv(signature_path)
    if {"cluster", "gene"}.issubset(df.columns):
        return df.groupby("cluster")["gene"].apply(lambda values: [str(v) for v in values]).to_dict()
    if {"cluster", "genes"}.issubset(df.columns):
        signatures = {}
        for _, row in df.iterrows():
            genes = [gene.strip() for gene in str(row["genes"]).replace(";", ",").split(",") if gene.strip()]
            signatures[str(row["cluster"])] = genes
        return signatures
    raise ValueError("signatures file must be JSON or CSV with columns cluster,gene or cluster,genes.")


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="surface-biomarkers")
    sub = p.add_subparsers(dest="command", required=True)

    train = sub.add_parser("train")
    train.add_argument("counts_csv")
    train.add_argument("--target", default="celltype")
    train.add_argument("--output", default="cluster_outputs")
    train.add_argument("--resampling", default="SMOTE")
    train.add_argument("--cv", type=int, default=5)
    train.add_argument("--lgbm-only", action="store_true")
    train.add_argument("--verbose", type=int, default=1)
    train.add_argument("--grid-verbose", type=int, default=0)

    discover = sub.add_parser("discover")
    discover.add_argument("base_dir")
    discover.add_argument("--clusters", nargs="+", required=True)
    discover.add_argument("--positive-genes", type=int, default=1)
    discover.add_argument("--negative-genes", type=int, default=1)
    discover.add_argument("--threshold", type=float, default=0.5)

    analyze = sub.add_parser("analyze")
    analyze.add_argument("base_dir")
    analyze.add_argument("--clusters", nargs="+", required=True)
    analyze.add_argument("--model", nargs="+", default=["LGBM"], choices=["LGBM", "RF"])
    analyze.add_argument("--max-display", type=int, default=10)

    validate = sub.add_parser("validate-signatures")
    validate.add_argument("base_dir")
    validate.add_argument("--ranking", default="gene_ranking_by_cluster.csv")
    validate.add_argument("--signatures", required=True)
    validate.add_argument("--surface-db", default="internal")
    validate.add_argument("--protein-atlas", default="internal")
    validate.add_argument("--output", default="signature_validation_priority.csv")
    validate.add_argument("--all-ranked", action="store_true")

    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = parser().parse_args(argv)

    if args.command == "train":
        df = load_counts_csv(args.counts_csv, args.target)
        config = TrainingConfig(
            resampling=args.resampling,
            cv=args.cv,
            verbose=args.verbose,
            grid_verbose=args.grid_verbose,
            train_random_forest=not args.lgbm_only,
            train_lightgbm=True,
        )
        print(train_one_vs_rest(df, args.target, args.output, config))
        return 0

    if args.command == "discover":
        base = Path(args.base_dir)
        config = DiscoveryConfig(
            positive_genes=args.positive_genes,
            negative_genes=args.negative_genes,
        )
        signatures, ranking, exclusivity = discover_signatures(base, args.clusters, config)
        shap_matrix = mean_abs_shap_matrix(base, args.clusters, config)
        combined_matrix = combined_shap_expression_matrix(base, args.clusters, config)
        accuracy = evaluate_signatures(base, args.clusters, signatures, args.threshold)

        with (base / "signatures.json").open("w", encoding="utf-8") as handle:
            json.dump(signatures, handle, indent=2)
        pd.DataFrame(
            [
                {"cluster": cluster, "gene": gene}
                for cluster, genes in signatures.items()
                for gene in genes
            ]
        ).to_csv(base / "signatures.csv", index=False)
        ranking.to_csv(base / "gene_ranking_by_cluster.csv", index=False)
        exclusivity.to_csv(base / "gene_exclusivity.csv", index=False)
        accuracy.to_csv(base / "signature_accuracy.csv", index=False)
        plot_shap_heatmap(
            shap_matrix,
            base / "shap_top_genes_heatmap.svg",
            title="Top genes by mean absolute SHAP",
            show=False,
        )
        plot_shap_heatmap(
            combined_matrix,
            base / "shap_expression_heatmap.svg",
            title="SHAP x expression score",
            show=False,
        )
        plot_signature_accuracy(accuracy, base / "signature_accuracy.svg", show=False)
        print("Signatures:")
        for cluster, genes in signatures.items():
            print(f"{cluster}: {', '.join(genes)}")
        return 0

    if args.command == "analyze":
        metrics = analyze_clusters(
            args.base_dir,
            args.clusters,
            model_name=args.model,
            max_display=args.max_display,
            show_panel=False,
        )
        print(metrics)
        return 0

    if args.command == "validate-signatures":
        base = Path(args.base_dir)
        ranking_path = Path(args.ranking)
        if not ranking_path.is_absolute():
            ranking_path = base / ranking_path
        output_path = Path(args.output)
        if not output_path.is_absolute():
            output_path = base / output_path
        ranking = pd.read_csv(ranking_path)
        signatures = _load_signatures(args.signatures)
        report = save_signature_validation_report(
            ranking=ranking,
            signatures=signatures,
            surface_db=args.surface_db,
            protein_atlas=None if args.protein_atlas.lower() == "none" else args.protein_atlas,
            output_path=output_path,
            signature_only=not args.all_ranked,
        )
        print(report)
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
