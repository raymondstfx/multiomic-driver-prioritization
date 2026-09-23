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
from multiomic_driver.utils.io import require_path
from multiomic_driver.visualization.ranking_plots import plot_zfpm2_hic2_vs_ntc

LOGGER = logging.getLogger("06b_zfpm2_external_validation")


def _zfpm2_vector(rna_path: Path) -> tuple[np.ndarray, pd.DataFrame]:
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
        metadata = rna.obs.copy()
        LOGGER.info("ZFPM2 is available as RNA feature %s", rna.var_names[matches[0]])
        return values, metadata
    finally:
        rna.file.close()


def main(config: dict) -> None:
    tables = Path(config["paths"]["results"]) / "tables"
    figures = Path(config["paths"]["results"]) / "figures"
    processed = Path(config["paths"]["processed"])

    ranking = pd.read_csv(require_path(tables / "candidate_ranking.csv"))
    settings = config["validation"]
    validate_ranking(
        ranking,
        expected_count=int(settings["expected_primary_candidates"]),
        required_candidates=settings["rankable_candidates"],
        excluded_candidates=settings["external_validation_genes"],
    )
    expression, metadata = _zfpm2_vector(
        require_path(processed / "rna_processed.h5ad", "processed RNA object")
    )
    summary = summarize_gene_expression(expression, metadata)
    contrasts = hic2_ntc_contrasts(summary)
    replicate_effects, did_summary = replicate_aware_did(contrasts)

    for table, filename in (
        (summary, "zfpm2_expression_by_group.csv"),
        (contrasts, "zfpm2_hic2_ntc_contrasts.csv"),
        (replicate_effects, "zfpm2_replicate_did.csv"),
        (did_summary, "zfpm2_validation_summary.csv"),
    ):
        path = tables / filename
        table.to_csv(path, index=False)
        LOGGER.info("Generated table %s (%d rows)", path, len(table))

    LOGGER.info(
        "ZFPM2 HIC2/NTC group sizes: %s",
        summary.set_index(["condition", "replicate", "group"])["n_cells"].to_dict(),
    )
    plot_zfpm2_hic2_vs_ntc(
        summary, contrasts, figures / "zfpm2_hic2_vs_ntc.png"
    )
    LOGGER.info("Generated figure %s", figures / "zfpm2_hic2_vs_ntc.png")


if __name__ == "__main__":
    run_stage("06b_zfpm2_external_validation", main)
