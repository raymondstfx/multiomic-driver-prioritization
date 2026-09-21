import pandas as pd

from multiomic_driver.ranking.multimodal_score import combine_multimodal_scores


def test_zscore_combination_and_deterministic_ranking() -> None:
    scores = pd.DataFrame(
        {
            "candidate": ["B", "A", "C"],
            "rna_score": [2.0, 2.0, 1.0],
            "atac_score": [2.0, 2.0, 1.0],
        }
    )
    result = combine_multimodal_scores(scores)
    assert result["candidate"].tolist() == ["A", "B", "C"]
    assert result["multiomic_rank"].tolist() == [1, 2, 3]

