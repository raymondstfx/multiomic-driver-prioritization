"""Simple, interpretable cross-modal score integration."""

from __future__ import annotations

import numpy as np
import pandas as pd


def zscore(values: pd.Series) -> pd.Series:
    """Standardise finite scores using population SD (``ddof=0``)."""
    numeric = pd.to_numeric(values, errors="raise").astype(float)
    if not np.isfinite(numeric).all():
        raise ValueError("Cannot standardize scores containing NaN or infinity")
    standard_deviation = numeric.std(ddof=0)
    if np.isclose(standard_deviation, 0.0):
        raise ValueError("Cannot standardize scores because score variance is zero")
    return (numeric - numeric.mean()) / standard_deviation


def _validate_candidates(table: pd.DataFrame, modality: str) -> set[str]:
    if "candidate" not in table:
        raise ValueError(f"{modality} table is missing candidate")
    if table["candidate"].isna().any():
        raise ValueError(f"{modality} candidates cannot be missing")
    if table["candidate"].duplicated().any():
        raise ValueError(f"{modality} table contains duplicate candidates")
    return set(table["candidate"].astype(str))


def combine_modality_rankings(
    rna_scores: pd.DataFrame, atac_scores: pd.DataFrame
) -> pd.DataFrame:
    """Validate, standardize, and equally combine modality rankings."""
    rna_candidates = _validate_candidates(rna_scores, "RNA")
    atac_candidates = _validate_candidates(atac_scores, "ATAC")
    if rna_candidates != atac_candidates:
        rna_only = sorted(rna_candidates - atac_candidates)
        atac_only = sorted(atac_candidates - rna_candidates)
        raise ValueError(
            f"RNA/ATAC candidate sets differ; RNA-only={rna_only}, "
            f"ATAC-only={atac_only}"
        )
    result = rna_scores.merge(
        atac_scores, on="candidate", how="inner", validate="one_to_one"
    )
    result["rna_z"] = zscore(result["rna_score"])
    result["atac_z"] = zscore(result["atac_score"])
    result["multiomic_score"] = (result["rna_z"] + result["atac_z"]) / 2.0
    result["multiomic_rank"] = result["multiomic_score"].rank(
        method="min", ascending=False
    ).astype(int)
    result = result.sort_values(
        ["multiomic_rank", "candidate"], ascending=[True, True], kind="stable"
    ).reset_index(drop=True)
    return result


def combine_multimodal_scores(scores: pd.DataFrame) -> pd.DataFrame:
    """Backward-compatible wrapper for a pre-merged score table."""
    required = {"candidate", "rna_score", "atac_score"}
    missing = required - set(scores.columns)
    if missing:
        raise ValueError(f"Missing score columns: {sorted(missing)}")
    return combine_modality_rankings(
        scores[["candidate", "rna_score"]],
        scores[["candidate", "atac_score"]],
    )
