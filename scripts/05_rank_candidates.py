"""Rank primary candidates from verified Phase 4 modality effects."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
from _common import run_stage

from multiomic_driver.ranking.atac_score import rank_atac_effects
from multiomic_driver.ranking.multimodal_score import combine_modality_rankings
from multiomic_driver.ranking.rna_score import rank_rna_effects
from multiomic_driver.utils.io import require_path
from multiomic_driver.visualization.ranking_plots import (
    plot_candidate_ranking,
    plot_modality_rank_comparison,
    plot_replicate_consistency,
    plot_rna_vs_atac,
    plot_rna_vs_atac_zscores,
)

LOGGER = logging.getLogger("05_rank_candidates")


def _eligible_candidates(table: pd.DataFrame) -> set[str]:
    values = table["phase4_eligible"]
    if not pd.api.types.is_bool_dtype(values.dtype):
        values = values.astype(str).str.lower().eq("true")
    return set(table.loc[values, "target_gene"].astype(str))


def main(config: dict) -> None:
    tables = Path(config["paths"]["results"]) / "tables"
    figures = Path(config["paths"]["results"]) / "figures"
    rna_path = require_path(tables / "rna_effect_summary.csv")
    atac_path = require_path(tables / "atac_effect_summary.csv")
    eligibility_path = require_path(tables / "candidate_phase4_eligibility.csv")
    LOGGER.info("RNA input: %s", rna_path)
    LOGGER.info("ATAC input: %s", atac_path)
    if config["ranking"]["combine_method"] != "mean_zscore":
        raise ValueError("Phase 5 primary ranking requires mean_zscore")

    eligibility = pd.read_csv(eligibility_path)
    expected = _eligible_candidates(eligibility)
    rna = rank_rna_effects(pd.read_csv(rna_path))
    atac = rank_atac_effects(pd.read_csv(atac_path))
    if set(rna["candidate"]) != expected or set(atac["candidate"]) != expected:
        raise ValueError("Phase 4 summaries do not match the Phase 3.5 primary set")
    ranking = combine_modality_rankings(rna, atac)
    external = set(config["validation"]["external_validation_genes"])
    if external & set(ranking["candidate"]):
        raise ValueError("External validation genes cannot be ranked perturbations")

    columns = [
        "candidate",
        "rna_score",
        "rna_z",
        "rna_rank",
        "rna_replicate_cosine",
        "rna_r1_effect_norm",
        "rna_r2_effect_norm",
        "atac_score",
        "atac_z",
        "atac_rank",
        "atac_replicate_cosine",
        "atac_r1_effect_norm",
        "atac_r2_effect_norm",
        "multiomic_score",
        "multiomic_rank",
    ]
    ranking = ranking[columns]
    output = tables / "candidate_ranking.csv"
    ranking.to_csv(output, index=False)
    diagnostics = ranking[
        [
            "candidate",
            "rna_r1_effect_norm",
            "rna_r2_effect_norm",
            "rna_replicate_cosine",
            "atac_r1_effect_norm",
            "atac_r2_effect_norm",
            "atac_replicate_cosine",
        ]
    ]
    diagnostics.to_csv(tables / "candidate_ranking_diagnostics.csv", index=False)

    rna_sd = ranking["rna_score"].std(ddof=0)
    atac_sd = ranking["atac_score"].std(ddof=0)
    summary = pd.DataFrame(
        [
            ("candidate_count", len(ranking)),
            ("rna_score_mean", ranking["rna_score"].mean()),
            ("rna_score_population_sd", rna_sd),
            ("atac_score_mean", ranking["atac_score"].mean()),
            ("atac_score_population_sd", atac_sd),
            (
                "rna_atac_spearman_correlation",
                ranking["rna_z"].corr(ranking["atac_z"], method="spearman"),
            ),
            ("top_rna_candidate", ranking.nsmallest(1, "rna_rank").iloc[0]["candidate"]),
            ("top_atac_candidate", ranking.nsmallest(1, "atac_rank").iloc[0]["candidate"]),
            ("top_multiomic_candidate", ranking.iloc[0]["candidate"]),
        ],
        columns=["metric", "value"],
    )
    summary.to_csv(tables / "ranking_summary.csv", index=False)

    plot_rna_vs_atac(ranking, figures / "rna_vs_atac_scores.png")
    plot_rna_vs_atac_zscores(ranking, figures / "rna_vs_atac_zscores.png")
    plot_candidate_ranking(ranking, figures / "candidate_ranking.png")
    plot_modality_rank_comparison(
        ranking, figures / "modality_rank_comparison.png"
    )
    plot_replicate_consistency(ranking, figures / "replicate_consistency.png")

    LOGGER.info("Ranked %d primary candidates", len(ranking))
    LOGGER.info(
        "RNA scores: range %.6g-%.6g, mean %.6g, population SD %.6g",
        ranking["rna_score"].min(),
        ranking["rna_score"].max(),
        ranking["rna_score"].mean(),
        rna_sd,
    )
    LOGGER.info(
        "ATAC scores: range %.6g-%.6g, mean %.6g, population SD %.6g",
        ranking["atac_score"].min(),
        ranking["atac_score"].max(),
        ranking["atac_score"].mean(),
        atac_sd,
    )
    LOGGER.info(
        "Top candidates: RNA=%s, ATAC=%s, multi-omic=%s",
        summary.set_index("metric").at["top_rna_candidate", "value"],
        summary.set_index("metric").at["top_atac_candidate", "value"],
        summary.set_index("metric").at["top_multiomic_candidate", "value"],
    )
    hic2 = ranking.loc[ranking["candidate"].eq("HIC2")]
    if not hic2.empty:
        row = hic2.iloc[0]
        hic2[
            [
                "candidate",
                "rna_rank",
                "atac_rank",
                "multiomic_rank",
                "rna_replicate_cosine",
                "atac_replicate_cosine",
            ]
        ].to_csv(tables / "phase5_hic2_evaluation.csv", index=False)
        LOGGER.info(
            "Post-hoc HIC2: RNA rank=%d, ATAC rank=%d, multi-omic rank=%d, "
            "RNA cosine=%.3f, ATAC cosine=%.3f",
            row["rna_rank"],
            row["atac_rank"],
            row["multiomic_rank"],
            row["rna_replicate_cosine"],
            row["atac_replicate_cosine"],
        )
    LOGGER.info("Output: %s", output)


if __name__ == "__main__":
    run_stage("05_rank_candidates", main)
