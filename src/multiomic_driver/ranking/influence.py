"""Experimental influence analogue; not used by the primary ranking."""

from __future__ import annotations

import numpy as np


def influence_analogue(effect_vector) -> float:
    """Sum absolute downstream effects as a transparent influence analogue."""
    return float(np.abs(np.asarray(effect_vector, dtype=float)).sum())
