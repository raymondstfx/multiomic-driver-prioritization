"""Difference-in-Differences background correction."""

from __future__ import annotations


def difference_in_differences(
    treated_perturbed,
    treated_control,
    untreated_perturbed,
    untreated_control,
):
    """Calculate a scalar or vector background-corrected treatment effect."""
    return (
        treated_perturbed
        - treated_control
        - untreated_perturbed
        + untreated_control
    )

