"""Simple, interpretable cross-modal score integration."""

from __future__ import annotations

import numpy as np
import pandas as pd


def zscore(values: pd.Series) -> pd.Series:
    """Standardise scores with population variance; constants map to zero."""
    numeric = values.astype(float)
    standard_deviation = numeric.std(ddof=0)
    if np.isclose(standard_deviation, 0.0):
        return pd.Series(0.0, index=values.index, dtype=float)
    return (numeric - numeric.mean()) / standard_deviation


def combine_multimodal_scores(scores: pd.DataFrame) -> pd.DataFrame:
    """Combine RNA and ATAC z-scores and assign deterministic descending ranks."""
    required = {"candidate", "rna_score", "atac_score"}
    missing = required - set(scores.columns)
    if missing:
        raise ValueError(f"Missing score columns: {sorted(missing)}")
    result = scores.copy()
    result["rna_z"] = zscore(result["rna_score"])
    result["atac_z"] = zscore(result["atac_score"])
    result["multiomic_score"] = (result["rna_z"] + result["atac_z"]) / 2.0
    result = result.sort_values(
        ["multiomic_score", "candidate"], ascending=[False, True], kind="stable"
    ).reset_index(drop=True)
    result["multiomic_rank"] = np.arange(1, len(result) + 1)
    return result

