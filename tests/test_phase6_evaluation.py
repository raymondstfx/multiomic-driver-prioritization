"""Synthetic tests for Phase 6 descriptive and external validation logic."""

import numpy as np
import pandas as pd
import pytest

from multiomic_driver.evaluation.ranking_evaluation import (
    assign_cross_modal_profile,
    hic2_validation_table,
    rank_comparison_table,
    validate_ranking,
)
from multiomic_driver.evaluation.zfpm2 import (
    hic2_ntc_contrasts,
    replicate_aware_did,
    summarize_gene_expression,
)


def _ranking() -> pd.DataFrame:
    candidates = ["A", "HIC2", "B"]
    return pd.DataFrame(
        {
            "candidate": candidates,
            "rna_score": [3.0, 2.0, 1.0],
            "atac_score": [2.0, 3.0, 1.0],
            "rna_z": [1.0, 0.5, -1.5],
            "atac_z": [0.5, 1.0, -1.5],
            "multiomic_score": [0.75, 0.75, -1.5],
            "rna_rank": [1, 2, 3],
            "atac_rank": [2, 1, 3],
            "multiomic_rank": [1, 2, 3],
            "rna_replicate_cosine": [0.8, 0.9, 0.1],
            "atac_replicate_cosine": [0.7, 0.6, -0.2],
        }
    )


def test_rank_shifts_and_profiles_are_derived() -> None:
    comparison = rank_comparison_table(_ranking())
    hic2 = comparison.set_index("candidate").loc["HIC2"]
    assert hic2["rna_to_multiomic_shift"] == 0
    assert hic2["atac_to_multiomic_shift"] == -1
    assert hic2["cross_modal_profile"] == "both_high"
    assert assign_cross_modal_profile(1.0, -0.1) == "rna_dominant"
    assert assign_cross_modal_profile(-0.1, 1.0) == "atac_dominant"
    assert assign_cross_modal_profile(0.0, 0.0) == "both_below_average"


def test_hic2_lookup_is_post_hoc_and_zfpm2_is_excluded() -> None:
    ranking = _ranking()
    validate_ranking(ranking, expected_count=3)
    hic2 = hic2_validation_table(ranking).iloc[0]
    assert hic2["multiomic_rank"] == 2
    assert hic2["validation_role"] == "experimentally supported pooled perturbation"

    invalid = pd.concat(
        [ranking.iloc[:2], ranking.iloc[[2]].assign(candidate="ZFPM2")],
        ignore_index=True,
    )
    with pytest.raises(ValueError, match="must not be ranked"):
        validate_ranking(invalid, expected_count=3)


def _expression_inputs() -> tuple[np.ndarray, pd.DataFrame]:
    rows = []
    values = []
    means = {
        ("DMSO", "1", "HIC2"): 2.0,
        ("DMSO", "1", "NTC"): 1.0,
        ("Dasatinib", "1", "HIC2"): 5.0,
        ("Dasatinib", "1", "NTC"): 2.0,
        ("DMSO", "2", "HIC2"): 3.0,
        ("DMSO", "2", "NTC"): 2.0,
        ("Dasatinib", "2", "HIC2"): 6.0,
        ("Dasatinib", "2", "NTC"): 3.0,
    }
    for (condition, replicate, group), value in means.items():
        for _ in range(2):
            rows.append(
                {
                    "condition": condition,
                    "replicate": replicate,
                    "guide_assignment_status": "singlet",
                    "target_gene": "HIC2" if group == "HIC2" else "NTC1",
                    "is_non_targeting": group == "NTC",
                }
            )
            values.append(value)
    return np.asarray(values), pd.DataFrame(rows)


def test_grouped_zfpm2_summary_and_replicate_aware_contrast() -> None:
    expression, metadata = _expression_inputs()
    summary = summarize_gene_expression(expression, metadata)
    assert len(summary) == 8
    assert set(summary["group"]) == {"HIC2", "NTC"}
    assert (summary["n_cells"] == 2).all()
    assert (summary["fraction_expressing"] == 1.0).all()

    contrasts = hic2_ntc_contrasts(summary)
    effects, did_summary = replicate_aware_did(contrasts)
    assert effects["replicate_specific_did"].tolist() == [2.0, 2.0]
    assert did_summary.iloc[0]["mean_replicate_did"] == 2.0
    assert bool(did_summary.iloc[0]["replicate_direction_agreement"])
