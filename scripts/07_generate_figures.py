"""Legacy convenience entry point; Phase 5 now generates these figures."""

import pandas as pd
from _common import run_stage

from multiomic_driver.utils.artifacts import figure_dir, table_dir
from multiomic_driver.utils.io import require_path
from multiomic_driver.visualization.plots import save_experimental_design
from multiomic_driver.visualization.ranking_plots import plot_rna_vs_atac


def main(config: dict) -> None:
    ranking = pd.read_csv(
        require_path(table_dir(config, "phase05") / "candidate_ranking.csv")
    )
    figures = figure_dir(config, "phase05")
    save_experimental_design(figures / "experimental_design.png")
    plot_rna_vs_atac(ranking, figures / "rna_vs_atac_scores.png")


if __name__ == "__main__":
    run_stage("07_generate_figures", main)
