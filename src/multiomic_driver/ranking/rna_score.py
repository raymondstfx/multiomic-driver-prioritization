"""RNA effect scoring."""

from __future__ import annotations

import numpy as np


def score_rna_effect(effect_vector) -> float:
    """Return the Euclidean magnitude of a background-corrected RNA effect."""
    return float(np.linalg.norm(np.asarray(effect_vector, dtype=float)))

