"""Tests for isolated Phase 7 robustness helpers."""

import numpy as np
import pandas as pd

from multiomic_driver.evaluation.sensitivity import (
    combine_sensitivity_scores,
    common_guide_sets,
    guide_composition,
    latent_centroid_effects,
)


def _metadata() -> pd.DataFrame:
    rows = []
    for condition in ("DMSO", "Dasatinib"):
        for replicate in (1, 2):
            for target, guide, ntc in (("A", "A_1", False), ("NTC1", "NTC1", True)):
                for _ in range(3):
                    rows.append(
                        {
                            "condition": condition,
                            "replicate": replicate,
                            "target_gene": target,
                            "guide": guide,
                            "is_non_targeting": ntc,
                            "guide_assignment_status": "singlet",
                        }
                    )
    return pd.DataFrame(rows)


def test_robust_effects_and_common_guide_sets() -> None:
    metadata = _metadata()
    matrix = np.arange(len(metadata) * 2, dtype=float).reshape(-1, 2)
    for method in ("mean", "winsorised_mean", "trimmed_mean", "median"):
        replicates, summary = latent_centroid_effects(
            matrix, metadata, ["A"], min_cells_per_group=3, method=method
        )
        assert len(replicates) == 2
        assert len(summary) == 1
        assert np.isfinite(summary.loc[0, "score"])
    composition = guide_composition(metadata)
    assert common_guide_sets(composition, min_cells_per_group=3) == {
        "A": {"A_1"},
        "NTC1": {"NTC1"},
    }


def test_weighted_sensitivity_combination_keeps_statistical_ties() -> None:
    result = combine_sensitivity_scores(
        pd.Series({"A": 2.0, "B": 1.0, "C": 1.0}),
        pd.Series({"A": 1.0, "B": 2.0, "C": 2.0}),
        rna_weight=0.5,
    )
    assert result["multiomic_rank"].tolist() == [1, 1, 1]
