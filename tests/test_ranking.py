"""Tests for label-free modality and multi-omic ranking."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from multiomic_driver.ranking.atac_score import rank_atac_effects
from multiomic_driver.ranking.multimodal_score import (
    combine_modality_rankings,
    zscore,
)
from multiomic_driver.ranking.rna_score import rank_rna_effects


def _effect_summary(candidates=("A", "B", "C")) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "target_gene": candidates,
            "mean_effect_norm": [1.0, 3.0, 3.0],
            "replicate_cosine_similarity": [-0.2, 0.8, 0.4],
            "r1_effect_norm": [1.1, 3.2, 2.8],
            "r2_effect_norm": [0.9, 2.8, 3.2],
        }
    )


def test_modality_ranking_uses_mean_norm_and_deterministic_ties() -> None:
    rna = rank_rna_effects(_effect_summary())
    atac = rank_atac_effects(_effect_summary().iloc[::-1].reset_index(drop=True))

    assert rna["candidate"].tolist() == ["B", "C", "A"]
    assert rna["rna_rank"].tolist() == [1, 1, 3]
    assert atac["candidate"].tolist() == ["B", "C", "A"]
    assert atac["atac_rank"].tolist() == [1, 1, 3]
    assert rna.set_index("candidate").at["A", "rna_score"] == 1.0
    assert rna.set_index("candidate").at["A", "rna_replicate_cosine"] == -0.2


def test_population_zscore_and_zero_variance_error() -> None:
    result = zscore(pd.Series([1.0, 2.0, 3.0]))
    np.testing.assert_allclose(result.mean(), 0.0, atol=1e-12)
    np.testing.assert_allclose(result.std(ddof=0), 1.0, atol=1e-12)
    np.testing.assert_allclose(result, [-np.sqrt(1.5), 0, np.sqrt(1.5)])
    with pytest.raises(ValueError, match="variance is zero"):
        zscore(pd.Series([2.0, 2.0, 2.0]))


def test_multiomic_combination_is_order_independent_and_equal_weighted() -> None:
    rna = pd.DataFrame(
        {
            "candidate": ["C", "A", "B"],
            "rna_score": [2.0, 1.0, 4.0],
            "rna_rank": [2, 3, 1],
        }
    )
    atac = pd.DataFrame(
        {
            "candidate": ["B", "C", "A"],
            "atac_score": [2.0, 5.0, 1.0],
            "atac_rank": [2, 1, 3],
        }
    )
    result = combine_modality_rankings(rna, atac)

    np.testing.assert_allclose(
        result["multiomic_score"], (result["rna_z"] + result["atac_z"]) / 2
    )
    shuffled = combine_modality_rankings(
        rna.sample(frac=1, random_state=3), atac.sample(frac=1, random_state=4)
    )
    pd.testing.assert_frame_equal(result, shuffled)


def test_candidate_mismatch_and_duplicates_fail_loudly() -> None:
    rna = pd.DataFrame({"candidate": ["A", "B"], "rna_score": [1.0, 2.0]})
    atac = pd.DataFrame({"candidate": ["A", "C"], "atac_score": [1.0, 2.0]})
    with pytest.raises(ValueError, match="candidate sets differ"):
        combine_modality_rankings(rna, atac)
    duplicated = pd.DataFrame(
        {"candidate": ["A", "A"], "rna_score": [1.0, 2.0]}
    )
    with pytest.raises(ValueError, match="duplicate candidates"):
        combine_modality_rankings(
            duplicated,
            pd.DataFrame(
                {"candidate": ["A", "B"], "atac_score": [1.0, 2.0]}
            ),
        )


def test_known_labels_are_not_used_or_inserted() -> None:
    ordinary = _effect_summary(("A", "B", "C"))
    labelled = ordinary.copy()
    labelled.loc[0, "target_gene"] = "HIC2"
    ordinary_scores = rank_rna_effects(ordinary)
    labelled_scores = rank_rna_effects(labelled)
    assert sorted(ordinary_scores["rna_score"]) == sorted(labelled_scores["rna_score"])

    rna = labelled_scores[["candidate", "rna_score"]]
    atac = rank_atac_effects(labelled)[["candidate", "atac_score"]]
    result = combine_modality_rankings(rna, atac)
    assert "HIC2" in set(result["candidate"])
    assert "ZFPM2" not in set(result["candidate"])


def test_nonfinite_primary_score_is_rejected() -> None:
    effects = _effect_summary()
    effects.loc[0, "mean_effect_norm"] = np.nan
    with pytest.raises(ValueError, match="finite"):
        rank_rna_effects(effects)
