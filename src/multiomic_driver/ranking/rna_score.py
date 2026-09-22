"""RNA effect scoring."""

from __future__ import annotations

import numpy as np
import pandas as pd


def score_rna_effect(effect_vector) -> float:
    """Return the Euclidean magnitude of a background-corrected RNA effect."""
    return float(np.linalg.norm(np.asarray(effect_vector, dtype=float)))


def rank_rna_effects(effect_summary: pd.DataFrame) -> pd.DataFrame:
    """Extract verified mean-effect norms and rank RNA candidates."""
    required = {
        "target_gene",
        "mean_effect_norm",
        "replicate_cosine_similarity",
        "r1_effect_norm",
        "r2_effect_norm",
    }
    missing = required - set(effect_summary.columns)
    if missing:
        raise ValueError(f"RNA effect summary is missing: {sorted(missing)}")
    if effect_summary["target_gene"].duplicated().any():
        raise ValueError("RNA effect summary contains duplicate candidates")
    scores = pd.to_numeric(effect_summary["mean_effect_norm"], errors="raise")
    if not np.isfinite(scores).all():
        raise ValueError("RNA scores must be finite")
    if scores.lt(0).any():
        raise ValueError("RNA effect norms cannot be negative")
    result = pd.DataFrame(
        {
            "candidate": effect_summary["target_gene"].astype(str),
            "rna_score": scores.astype(float),
            "rna_replicate_cosine": effect_summary[
                "replicate_cosine_similarity"
            ].astype(float),
            "rna_r1_effect_norm": effect_summary["r1_effect_norm"].astype(float),
            "rna_r2_effect_norm": effect_summary["r2_effect_norm"].astype(float),
        }
    )
    result["rna_rank"] = result["rna_score"].rank(
        method="min", ascending=False
    ).astype(int)
    return result.sort_values(
        ["rna_rank", "candidate"], kind="stable"
    ).reset_index(drop=True)
