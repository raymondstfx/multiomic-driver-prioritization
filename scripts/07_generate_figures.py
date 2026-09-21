"""Generate figures from completed result tables."""

from pathlib import Path

import pandas as pd
from _common import run_stage

from multiomic_driver.utils.io import require_path
from multiomic_driver.visualization.plots import save_experimental_design
from multiomic_driver.visualization.ranking_plots import plot_rna_vs_atac


def main(config: dict) -> None:
    results = Path(config["paths"]["results"])
    ranking = pd.read_csv(require_path(results / "tables" / "candidate_ranking.csv"))
    figures = results / "figures"
    save_experimental_design(figures / "experimental_design.png")
    plot_rna_vs_atac(ranking, figures / "rna_vs_atac_scores.png")
    # TODO: add ranking and validated-rank panels after final table columns are fixed.


if __name__ == "__main__":
    run_stage("07_generate_figures", main)

