"""Surface biomarker discovery for single-cell expression matrices."""

from surface_biomarkers.data import load_counts_csv, make_one_vs_rest, sanitize_name
from surface_biomarkers.signatures import DiscoveryConfig

__all__ = [
    "DiscoveryConfig",
    "TrainingConfig",
    "design_antibody_panel",
    "discover_signatures",
    "evaluate_signatures",
    "load_counts_csv",
    "make_one_vs_rest",
    "prioritize_signature_validation",
    "sanitize_name",
    "save_signature_validation_report",
    "train_one_vs_rest",
]


def __getattr__(name: str):
    if name in {"TrainingConfig", "train_one_vs_rest"}:
        from surface_biomarkers.training import TrainingConfig, train_one_vs_rest

        return {"TrainingConfig": TrainingConfig, "train_one_vs_rest": train_one_vs_rest}[name]
    if name in {"discover_signatures", "evaluate_signatures"}:
        from surface_biomarkers.signatures import discover_signatures, evaluate_signatures

        return {
            "discover_signatures": discover_signatures,
            "evaluate_signatures": evaluate_signatures,
        }[name]
    if name in {"prioritize_signature_validation", "save_signature_validation_report"}:
        from surface_biomarkers.validation import (
            prioritize_signature_validation,
            save_signature_validation_report,
        )

        return {
            "prioritize_signature_validation": prioritize_signature_validation,
            "save_signature_validation_report": save_signature_validation_report,
        }[name]
    if name == "design_antibody_panel":
        from surface_biomarkers.panel_design import design_antibody_panel

        return design_antibody_panel
    raise AttributeError(name)
