"""Combine modality-specific effect magnitudes into a candidate ranking."""

from pathlib import Path

import pandas as pd
from _common import run_stage

from multiomic_driver.ranking.multimodal_score import combine_multimodal_scores
from multiomic_driver.utils.io import require_path


def main(config: dict) -> None:
    tables = Path(config["paths"]["results"]) / "tables"
    rna = pd.read_csv(require_path(tables / "rna_scores.csv"))
    atac = pd.read_csv(require_path(tables / "atac_scores.csv"))
    merged = rna.merge(atac, on="candidate", validate="one_to_one")
    ranking = combine_multimodal_scores(merged)
    tables.mkdir(parents=True, exist_ok=True)
    ranking.to_csv(tables / "candidate_ranking.csv", index=False)


if __name__ == "__main__":
    run_stage("05_rank_candidates", main)

