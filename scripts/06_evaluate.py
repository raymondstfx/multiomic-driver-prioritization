"""Evaluate externally supported candidates after ranking."""

from pathlib import Path

import pandas as pd
from _common import run_stage

from multiomic_driver.evaluation.validation import evaluate_known_candidates
from multiomic_driver.utils.io import require_path


def main(config: dict) -> None:
    tables = Path(config["paths"]["results"]) / "tables"
    ranking = pd.read_csv(require_path(tables / "candidate_ranking.csv"))
    summary = evaluate_known_candidates(
        ranking, config["validation"]["rankable_candidates"]
    )
    summary.to_csv(tables / "validation_summary.csv", index=False)


if __name__ == "__main__":
    run_stage("06_evaluate", main)
