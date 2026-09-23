"""Compatibility helper for simple contrasts; not used by the primary workflow."""

from __future__ import annotations


def perturbation_effect(perturbed, non_targeting):
    """Return the perturbation-minus-control contrast."""
    return perturbed - non_targeting
