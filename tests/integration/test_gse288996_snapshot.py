"""Optional regression snapshot for locally generated real-data artifacts."""

from pathlib import Path

import pandas as pd
import pytest

pytestmark = pytest.mark.real_data


def test_gse288996_primary_snapshot() -> None:
    root = Path(__file__).resolve().parents[2]
    ranking_path = root / "results/summary/final_candidate_ranking.csv"
    metadata_path = root / "data/interim/cell_metadata.csv"
    if not ranking_path.exists() or not metadata_path.exists():
        pytest.skip("Local real-data outputs are unavailable")
    ranking = pd.read_csv(ranking_path)
    metadata = pd.read_csv(metadata_path, usecols=["cell_id", "target_gene"])
    assert len(metadata) == 53_898
    assert metadata["cell_id"].is_unique
    assert len(ranking) == 12
    assert ranking.nsmallest(5, "display_order")["candidate"].tolist() == [
        "YEATS4",
        "GPBP1L1",
        "ZBED6",
        "HIC2",
        "KMT2B",
    ]
    assert "ZFPM2" not in set(ranking["candidate"])
