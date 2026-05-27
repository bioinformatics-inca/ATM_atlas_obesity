import json
import pickle
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import pandas as pd
from imblearn.combine import SMOTEENN
from imblearn.over_sampling import ADASYN, SMOTE, RandomOverSampler
from imblearn.pipeline import Pipeline
from lightgbm import LGBMClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.preprocessing import StandardScaler

from surface_biomarkers.data import make_one_vs_rest, sanitize_name


@dataclass
class TrainingConfig:
    resampling: str = "SMOTE"
    test_size: float = 0.2
    random_state: int = 42
    cv: int = 5
    scoring: str = "roc_auc"
    n_jobs: int = -1
    verbose: int = 1
    grid_verbose: int = 0
    train_random_forest: bool = True
    train_lightgbm: bool = True
    rf_param_grid: Dict[str, List[Any]] = field(default_factory=lambda: {
        "n_estimators": [100, 500],
        "max_features": ["sqrt", "log2"],
        "max_depth": [5, 10],
        "min_samples_split": [2, 5],
        "min_samples_leaf": [1, 2],
        "bootstrap": [True, False],
    })
    lgbm_param_grid: Dict[str, List[Any]] = field(default_factory=lambda: {
        "learning_rate": [0.01, 0.05],
        "n_estimators": [60, 100],
        "num_leaves": [8, 12],
        "objective": ["binary"],
        "metric": ["binary_logloss"],
    })


def get_resampler(name: str, random_state: int):
    options = {
        "RandomOverSampler": RandomOverSampler,
        "SMOTE": SMOTE,
        "ADASYN": ADASYN,
        "SMOTEENN": SMOTEENN,
        "none": None,
        "None": None,
    }
    if name not in options:
        raise ValueError(f"Unknown resampling method: {name}")
    cls = options[name]
    return None if cls is None else cls(random_state=random_state)


def _save_pickle(obj: object, path: Path) -> None:
    with path.open("wb") as handle:
        pickle.dump(obj, handle)


def _log(config: TrainingConfig, message: str, level: int = 1) -> None:
    if config.verbose >= level:
        print(message, flush=True)


def _prefix_param_grid(param_grid: Dict[str, List[Any]], step_name: str = "model") -> Dict[str, List[Any]]:
    return {
        key if key.startswith(f"{step_name}__") else f"{step_name}__{key}": value
        for key, value in param_grid.items()
    }


def _build_training_pipeline(estimator: object, resampler: Optional[object]) -> Pipeline:
    # Resampling and scaling must live inside the CV pipeline to avoid leaking
    # synthetic samples or preprocessing statistics across validation folds.
    steps = []
    if resampler is not None:
        steps.append(("resampler", resampler))
    steps.extend([
        ("scaler", StandardScaler()),
        ("model", estimator),
    ])
    return Pipeline(steps)


def train_cluster_models(
    binary_df: pd.DataFrame,
    cluster_name: str,
    output_dir: Union[str, Path],
    class_column: str = "Class",
    config: Optional[TrainingConfig] = None,
) -> Dict[str, Any]:
    config = config or TrainingConfig()
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    _log(config, f"\n[cluster {cluster_name}] preparing data")

    X = binary_df.drop(columns=[class_column])
    y = binary_df[class_column]
    positives = int(y.sum())
    negatives = int((y == 0).sum())
    _log(config, f"[cluster {cluster_name}] cells: {len(y)} | positives: {positives} | negatives: {negatives}")

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=config.test_size,
        random_state=config.random_state,
        stratify=y,
    )

    X_train.to_csv(output / "X_train.csv")
    X_test.to_csv(output / "X_test.csv")
    y_train.to_csv(output / "y_train.csv")
    y_test.to_csv(output / "y_test.csv")

    _log(config, f"[cluster {cluster_name}] resampling: {config.resampling}")

    jobs: List[Tuple[str, object, Dict[str, List[Any]], str]] = []
    if config.train_random_forest:
        jobs.append((
            "RF",
            RandomForestClassifier(random_state=config.random_state),
            config.rf_param_grid,
            f"rf_model_{cluster_name}.pkl",
        ))
    if config.train_lightgbm:
        jobs.append((
            "LGBM",
            LGBMClassifier(random_state=config.random_state, verbosity=-1),
            config.lgbm_param_grid,
            f"lgb_model_{cluster_name}.pkl",
        ))

    row: Dict[str, Any] = {"Cluster": cluster_name}
    for model_label, estimator, grid, filename in jobs:
        _log(config, f"[cluster {cluster_name}] training {model_label} with GridSearchCV")
        pipeline = _build_training_pipeline(estimator, get_resampler(config.resampling, config.random_state))
        search = GridSearchCV(
            pipeline,
            _prefix_param_grid(grid),
            cv=config.cv,
            scoring=config.scoring,
            n_jobs=config.n_jobs,
            verbose=config.grid_verbose,
        )
        search.fit(X_train, y_train)
        y_pred = search.predict(X_test)
        y_prob = search.predict_proba(X_test)[:, 1]
        row[f"{model_label}_AUC"] = roc_auc_score(y_test, y_prob)
        row[f"{model_label}_Accuracy"] = accuracy_score(y_test, y_pred)
        row[f"{model_label}_BestParams"] = {
            key.replace("model__", "", 1): value
            for key, value in search.best_params_.items()
        }
        best_pipeline = search.best_estimator_
        _save_pickle(best_pipeline.named_steps["scaler"], output / "scaler.pkl")
        _save_pickle(best_pipeline.named_steps["model"], output / filename)
        _log(
            config,
            (
                f"[cluster {cluster_name}] {model_label} done | "
                f"AUC={row[f'{model_label}_AUC']:.3f} | "
                f"Accuracy={row[f'{model_label}_Accuracy']:.3f}"
            ),
        )

    with (output / "training_config.json").open("w", encoding="utf-8") as handle:
        json.dump(asdict(config), handle, indent=2)
    _log(config, f"[cluster {cluster_name}] saved outputs to {output}")
    return row


def train_one_vs_rest(
    dataset: pd.DataFrame,
    target_column: str = "celltype",
    output_dir: Union[str, Path] = "cluster_outputs",
    config: Optional[TrainingConfig] = None,
) -> pd.DataFrame:
    config = config or TrainingConfig()
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    label_mapping = dataset.attrs.get("label_mapping")
    if label_mapping is not None:
        label_mapping.to_csv(output / "label_mapping.csv", index=False)
        _log(config, f"Saved label mapping to {output / 'label_mapping.csv'}")
    rows = []
    labels = sorted(dataset[target_column].unique())
    _log(config, f"Training {len(labels)} one-vs-rest models into {output}")

    for label in labels:
        cluster_name = sanitize_name(label)
        _log(config, f"\n=== {cluster_name} vs all ===")
        binary = make_one_vs_rest(dataset, label, target_column=target_column, class_column="Class")
        row = train_cluster_models(binary, cluster_name, output / cluster_name, config=config)
        row["OriginalCluster"] = label
        rows.append(row)

    results = pd.DataFrame(rows)
    results.to_csv(output / "model_metrics.csv", index=False)
    _log(config, f"\nFinished. Summary saved to {output / 'model_metrics.csv'}")
    return results
