"""Evaluate fixed Phase 5 rankings without tuning or recomputing their scores."""

import logging
from pathlib import Path

import pandas as pd
from _common import run_stage

from multiomic_driver.evaluation.ranking_evaluation import (
    final_candidate_evaluation,
    hic2_validation_table,
    literature_evidence_table,
    rank_comparison_table,
    rank_correlation_summary,
    replicate_reliability_table,
    top_candidate_comparison,
    validate_ranking,
)
from multiomic_driver.utils.artifacts import figure_dir, table_dir
from multiomic_driver.utils.io import load_yaml, require_path
from multiomic_driver.visualization.ranking_plots import plot_top5_candidate_evidence

LOGGER = logging.getLogger("06_evaluate")


def _write(table: pd.DataFrame, path: Path) -> None:
    table.to_csv(path, index=False)
    LOGGER.info("Generated table %s (%d rows)", path, len(table))


def main(config: dict) -> None:
    ranking_tables = table_dir(config, "phase05")
    tables = table_dir(config, "phase06")
    figures = figure_dir(config, "phase06")

    ranking_path = require_path(
        ranking_tables / "candidate_ranking.csv", "Phase 5 ranking"
    )
    ranking = pd.read_csv(ranking_path)
    settings = config["validation"]
    validate_ranking(
        ranking,
        expected_count=int(settings["expected_primary_candidates"]),
        required_candidates=settings["rankable_candidates"],
        excluded_candidates=settings["external_validation_genes"],
    )
    LOGGER.info(
        "Validated ranking input %s with %d candidates", ranking_path, len(ranking)
    )

    comparison = rank_comparison_table(ranking)
    correlations = rank_correlation_summary(ranking)
    reliability = replicate_reliability_table(
        ranking,
        moderate_threshold=float(settings["replicate_consistency_moderate"]),
        positive_threshold=float(settings["replicate_consistency_positive"]),
    )
    hic2 = hic2_validation_table(ranking)
    top5 = top_candidate_comparison(ranking, int(settings["top_n_candidates"]))

    evidence_path = Path(config["_project_root"]) / settings["literature_evidence_path"]
    evidence_config = load_yaml(evidence_path)
    top5_with_ranks = top5.merge(
        ranking[["candidate", "rna_rank", "atac_rank"]], on="candidate"
    )
    literature = literature_evidence_table(top5_with_ranks, evidence_config)
    final = final_candidate_evaluation(ranking, comparison, literature)

    for table, filename in (
        (comparison, "modality_rank_comparison.csv"),
        (correlations, "rank_correlation_summary.csv"),
        (reliability, "replicate_reliability_summary.csv"),
        (hic2, "hic2_validation_summary.csv"),
        (top5, "top5_candidate_comparison.csv"),
        (literature, "top_candidate_literature_evidence.csv"),
        (final, "final_candidate_evaluation.csv"),
    ):
        _write(table, tables / filename)

    plot_top5_candidate_evidence(top5, figures / "top5_candidate_evidence.png")
    LOGGER.info("Generated figure %s", figures / "top5_candidate_evidence.png")
    LOGGER.info("Top candidates: %s", ", ".join(top5["candidate"]))
    LOGGER.info(
        "HIC2 ranks: RNA=%d, ATAC=%d, multi-omic=%d",
        hic2.iloc[0]["rna_rank"],
        hic2.iloc[0]["atac_rank"],
        hic2.iloc[0]["multiomic_rank"],
    )
    LOGGER.info(
        "Rank correlations: %s",
        correlations.set_index("comparison")["spearman_rho"].round(3).to_dict(),
    )


if __name__ == "__main__":
    run_stage("06_evaluate", main)
