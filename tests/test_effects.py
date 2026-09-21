import numpy as np
import pandas as pd
import pytest

from multiomic_driver.effects.background_corrected import difference_in_differences
from multiomic_driver.effects.replicate_aware import (
    InsufficientGroupSizeError,
    MissingExperimentalGroupError,
    candidate_replicate_effects,
    compute_modality_effects,
    validate_eligibility_for_effects,
)


def test_difference_in_differences() -> None:
    result = difference_in_differences(
        treated_perturbed=10,
        treated_control=4,
        untreated_perturbed=5,
        untreated_control=3,
    )
    assert result == 4


def test_vector_difference_in_differences() -> None:
    result = difference_in_differences(
        np.array([10, 2]), np.array([4, 1]), np.array([5, 4]), np.array([3, 2])
    )
    np.testing.assert_array_equal(result, np.array([4, -1]))


def _replicate_design() -> tuple[np.ndarray, pd.DataFrame]:
    rows = []
    vectors = []
    values = {
        (1, "DMSO", "GENE"): [5, 2],
        (1, "DMSO", "NTC"): [3, 1],
        (1, "Dasatinib", "GENE"): [10, 4],
        (1, "Dasatinib", "NTC"): [4, 1],
        (2, "DMSO", "GENE"): [4, 2],
        (2, "DMSO", "NTC"): [1, 1],
        (2, "Dasatinib", "GENE"): [8, 1],
        (2, "Dasatinib", "NTC"): [2, 1],
    }
    for (replicate, condition, target), vector in values.items():
        rows.append(
            {
                "condition": condition,
                "replicate": replicate,
                "guide_assignment_status": "singlet",
                "target_gene": target,
                "is_non_targeting": target == "NTC",
            }
        )
        vectors.append(vector)
    return np.asarray(vectors, dtype=float), pd.DataFrame(rows)


def test_replicate_effects_are_computed_before_summarizing() -> None:
    matrix, metadata = _replicate_design()
    replicate_effects, summary = compute_modality_effects(
        matrix,
        metadata,
        ["GENE"],
        min_cells_per_group=1,
        dimension_prefix="PC",
    )

    np.testing.assert_array_equal(
        replicate_effects[["PC1", "PC2"]].to_numpy(), [[4, 2], [3, -1]]
    )
    np.testing.assert_allclose(summary[["PC1", "PC2"]].iloc[0], [3.5, 0.5])
    np.testing.assert_allclose(
        summary.at[0, "replicate_cosine_similarity"], 1 / np.sqrt(2)
    )


def test_missing_and_below_minimum_groups_raise_distinct_errors() -> None:
    matrix, metadata = _replicate_design()
    missing = ~(
        metadata["target_gene"].eq("GENE")
        & metadata["condition"].eq("Dasatinib")
        & metadata["replicate"].eq(1)
    )
    with pytest.raises(MissingExperimentalGroupError):
        candidate_replicate_effects(
            matrix[missing],
            metadata.loc[missing].reset_index(drop=True),
            "GENE",
            min_cells_per_group=1,
        )
    with pytest.raises(InsufficientGroupSizeError):
        candidate_replicate_effects(
            matrix,
            metadata,
            "GENE",
            min_cells_per_group=2,
        )


def test_effect_gate_rejects_missing_group_converted_to_zero() -> None:
    eligibility = pd.DataFrame(
        {
            "target_gene": ["MISSING", "LOW"],
            "dmso_r1_cells": [25, 18],
            "dmso_r2_cells": [25, 25],
            "dasatinib_r1_cells": [pd.NA, 25],
            "dasatinib_r2_cells": [25, 25],
            "min_group_cells": [pd.NA, 18],
            "all_four_groups_present": [False, True],
            "phase4_eligible": [False, False],
            "exclusion_reason": [
                "missing_Dasatinib_R1",
                "below_min_DMSO_R1(18<20)",
            ],
        }
    )
    validate_eligibility_for_effects(eligibility, min_cells_per_group=20)
    eligibility.loc[0, "dasatinib_r1_cells"] = 0
    eligibility.loc[0, "min_group_cells"] = 0
    with pytest.raises(ValueError, match="Low-count exclusion lacks reason"):
        validate_eligibility_for_effects(eligibility, min_cells_per_group=20)
