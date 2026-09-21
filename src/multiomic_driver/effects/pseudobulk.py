"""Aggregate cells into perturbation-treatment-replicate units."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd


def aggregate_pseudobulk(
    matrix,
    metadata: pd.DataFrame,
    groupby: Sequence[str],
    *,
    method: str = "mean",
) -> tuple[np.ndarray, pd.DataFrame]:
    """Aggregate rows of a matrix for each metadata group."""
    if matrix.shape[0] != len(metadata):
        raise ValueError("Matrix rows must match metadata rows")
    missing = set(groupby) - set(metadata.columns)
    if missing:
        raise ValueError(f"Missing grouping columns: {sorted(missing)}")
    if method not in {"mean", "sum"}:
        raise ValueError("method must be 'mean' or 'sum'")
    vectors: list[np.ndarray] = []
    labels: list[tuple] = []
    grouped = metadata.groupby(list(groupby), sort=True, observed=True).indices
    for label, indices in grouped.items():
        block = matrix[indices]
        aggregate = block.mean(axis=0) if method == "mean" else block.sum(axis=0)
        vectors.append(np.asarray(aggregate).ravel())
        labels.append(label if isinstance(label, tuple) else (label,))
    return np.vstack(vectors), pd.DataFrame(labels, columns=groupby)

