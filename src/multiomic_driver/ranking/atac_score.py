"""ATAC effect scoring."""

from __future__ import annotations

import numpy as np


def score_atac_effect(effect_vector) -> float:
    """Return the Euclidean magnitude of a background-corrected ATAC effect."""
    return float(np.linalg.norm(np.asarray(effect_vector, dtype=float)))

