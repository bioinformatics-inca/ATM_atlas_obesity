import pickle
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union

import pandas as pd


@dataclass
class ClusterArtifacts:
    cluster: str
    path: Path
    X_train: pd.DataFrame
    X_test: pd.DataFrame
    y_test: pd.Series
    scaler: object
    model: object


def load_cluster_artifacts(
    cluster_dir: Union[str, Path],
    cluster: Optional[str] = None,
    model_name: str = "LGBM",
) -> ClusterArtifacts:
    path = Path(cluster_dir)
    cluster_id = str(cluster or path.name)
    X_train = pd.read_csv(path / "X_train.csv", index_col=0)
    X_test = pd.read_csv(path / "X_test.csv", index_col=0)
    y_test = pd.read_csv(path / "y_test.csv", index_col=0).squeeze("columns")

    with (path / "scaler.pkl").open("rb") as handle:
        scaler = pickle.load(handle)

    prefix = "lgb" if model_name.upper() == "LGBM" else "rf"
    with (path / f"{prefix}_model_{cluster_id}.pkl").open("rb") as handle:
        model = pickle.load(handle)

    return ClusterArtifacts(cluster_id, path, X_train, X_test, y_test, scaler, model)
