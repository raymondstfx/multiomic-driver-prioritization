"""Evaluate downstream ZFPM2 expression without assigning it a perturbation rank."""

import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from _common import run_stage
from scipy import sparse

from multiomic_driver.evaluation.ranking_evaluation import validate_ranking
from multiomic_driver.evaluation.zfpm2 import (
    hic2_ntc_contrasts,
    replicate_aware_did,
    summarize_gene_expression,
)
from multiomic_driver.utils.artifacts import figure_dir, table_dir
from multiomic_driver.utils.io import require_path
from multiomic_driver.visualization.ranking_plots import plot_zfpm2_hic2_vs_ntc

LOGGER = logging.getLogger("06b_zfpm2_external_validation")


def _zfpm2_vectors(
    rna_path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, pd.DataFrame]:
    rna = ad.read_h5ad(rna_path, backed="r")
    try:
        if "feature_name" not in rna.var.columns:
            raise ValueError("RNA feature metadata does not contain feature_name")
        matches = np.flatnonzero(rna.var["feature_name"].astype(str).eq("ZFPM2"))
        if len(matches) != 1:
            raise ValueError(f"Expected one ZFPM2 RNA feature, found {len(matches)}")
        column = rna[:, matches].X
        values = (
            column.toarray().ravel()
            if sparse.issparse(column)
            else np.asarray(column).ravel()
        )
        counts_column = rna[:, matches].layers["counts"]
        counts = (
            counts_column.toarray().ravel()
            if sparse.issparse(counts_column)
            else np.asarray(counts_column).ravel()
        )
        metadata = rna.obs.copy()
        libraries = metadata["total_counts"].to_numpy(dtype=float)
        LOGGER.info("ZFPM2 is available as RNA feature %s", rna.var_names[matches[0]])
        return values, counts, libraries, metadata
    finally:
        rna.file.close()


def main(config: dict) -> None:
    ranking_tables = table_dir(config, "phase05")
    tables = table_dir(config, "phase06")
    figures = figure_dir(config, "phase06")
    processed = Path(config["paths"]["processed"])

    ranking = pd.read_csv(require_path(ranking_tables / "candidate_ranking.csv"))
    settings = config["validation"]
    validate_ranking(
        ranking,
        expected_count=int(settings["expected_primary_candidates"]),
        required_candidates=settings["rankable_candidates"],
        excluded_candidates=settings["external_validation_genes"],
    )
    expression, counts, libraries, metadata = _zfpm2_vectors(
        require_path(processed / "rna_processed.h5ad", "processed RNA object")
    )
    summary = summarize_gene_expression(
        expression, metadata, raw_counts=counts, library_sizes=libraries
    )
    metrics = [
        "mean_expression",
        "median_expression",
        "fraction_expressing",
        "mean_positive_expression",
        "pseudobulk_log1p_cpm",
    ]
    contrast_tables = []
    effect_tables = []
    did_tables = []
    for metric in metrics:
        metric_contrasts = hic2_ntc_contrasts(summary, metric=metric)
        metric_effects, metric_summary = replicate_aware_did(metric_contrasts)
        contrast_tables.append(metric_contrasts)
        effect_tables.append(metric_effects)
        did_tables.append(metric_summary)
    contrasts = pd.concat(contrast_tables, ignore_index=True)
    replicate_effects = pd.concat(effect_tables, ignore_index=True)
    did_summary = pd.concat(did_tables, ignore_index=True)

    candidate_context = []
    for candidate in ranking["candidate"].astype(str):
        candidate_summary = summarize_gene_expression(
            expression,
            metadata,
            target_candidate=candidate,
            raw_counts=counts,
            library_sizes=libraries,
        )
        candidate_contrasts = hic2_ntc_contrasts(
            candidate_summary,
            metric="mean_expression",
            target_candidate=candidate,
        )
        _, candidate_did = replicate_aware_did(candidate_contrasts)
        candidate_context.append(candidate_did)
    context = pd.concat(candidate_context, ignore_index=True)
    context["absolute_mean_replicate_did"] = context["mean_replicate_did"].abs()
    context["absolute_did_rank"] = context["absolute_mean_replicate_did"].rank(
        method="min", ascending=False
    ).astype(int)
    context = context.sort_values(["absolute_did_rank", "target_perturbation"])

    for table, filename in (
        (summary, "zfpm2_expression_by_group.csv"),
        (contrasts, "zfpm2_hic2_ntc_contrasts.csv"),
        (replicate_effects, "zfpm2_replicate_did.csv"),
        (did_summary, "zfpm2_validation_summary.csv"),
        (context, "zfpm2_across_candidate_context.csv"),
    ):
        path = tables / filename
        table.to_csv(path, index=False)
        LOGGER.info("Generated table %s (%d rows)", path, len(table))

    LOGGER.info(
        "ZFPM2 HIC2/NTC group sizes: %s",
        summary.set_index(["condition", "replicate", "group"])["n_cells"].to_dict(),
    )
    plot_zfpm2_hic2_vs_ntc(
        summary,
        contrasts.loc[contrasts["metric"].eq("mean_expression")],
        figures / "zfpm2_hic2_vs_ntc.png",
    )
    LOGGER.info("Generated figure %s", figures / "zfpm2_hic2_vs_ntc.png")


if __name__ == "__main__":
    run_stage("06b_zfpm2_external_validation", main)
