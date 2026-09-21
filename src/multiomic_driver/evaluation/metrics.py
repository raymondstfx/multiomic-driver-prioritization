"""Small-sample ranking diagnostics."""

from __future__ import annotations

from collections.abc import Iterable


def top_k_recovery(ranks: Iterable[int | None], k: int) -> float:
    """Return the fraction of supplied candidates observed within rank k."""
    values = list(ranks)
    if not values:
        raise ValueError("At least one candidate rank is required")
    return sum(rank is not None and rank <= k for rank in values) / len(values)

