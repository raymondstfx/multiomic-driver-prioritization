"""Preprocess the merged ATAC AnnData object."""

import logging
from pathlib import Path

import anndata as ad
from _common import run_stage

from multiomic_driver.preprocessing.atac import preprocess_atac
from multiomic_driver.preprocessing.qc import (
    summary_table,
    update_cross_modality_tables,
    validate_metadata,
)
from multiomic_driver.utils.io import require_path
from multiomic_driver.visualization.plots import save_embedding

LOGGER = logging.getLogger("03_preprocess_atac")


def main(config: dict) -> None:
    source = require_path(Path(config["paths"]["interim"]) / "atac_merged.h5ad")
    LOGGER.info("Input: %s", source)
    result = ad.read_h5ad(source)
    validate_metadata(result, "ATAC")
    cells_before, features_before = result.shape
    LOGGER.info("Input shape: %s", result.shape)
    LOGGER.info("ATAC thresholds: %s", config["atac"])
    settings = {**config["atac"], "random_seed": config["random_seed"]}
    result = preprocess_atac(result, **settings, copy=False)
    validate_metadata(result, "ATAC")

    tables = Path(config["paths"]["results"]) / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    summary_table(
        {
            "cells_before": cells_before,
            "cells_after": result.n_obs,
            "cells_removed": cells_before - result.n_obs,
            "input_called_peaks": features_before,
            "consensus_peaks_before_filtering": result.uns["preprocessing"][
                "consensus_features"
            ],
            "features_after_filtering": result.n_vars,
            "input_peaks_collapsed_by_harmonization": features_before
            - result.uns["preprocessing"]["consensus_features"],
            "consensus_features_removed_by_filtering": result.uns[
                "preprocessing"
            ]["consensus_features"]
            - result.n_vars,
            "lsi_dimensions": result.obsm["X_lsi"].shape[1],
            "lsi1_retained": True,
        }
    ).to_csv(tables / "atac_qc_summary.csv", index=False)
    retention = update_cross_modality_tables(result, "ATAC", tables)
    if (retention["cell_count"] == 0).any():
        raise ValueError("ATAC QC lost at least one HIC2 or NTC condition/replicate group")

    rna_output = Path(config["paths"]["processed"]) / "rna_processed.h5ad"
    if rna_output.exists():
        rna = ad.read_h5ad(rna_output, backed="r")
        overlap = len(result.obs_names.intersection(rna.obs_names))
        summary_table(
            {
                "rna_cells": rna.n_obs,
                "atac_cells": result.n_obs,
                "shared_cell_ids": overlap,
                "rna_only_cell_ids": rna.n_obs - overlap,
                "atac_only_cell_ids": result.n_obs - overlap,
            }
        ).to_csv(tables / "post_qc_cell_overlap.csv", index=False)
        rna.file.close()

    figures = Path(config["paths"]["results"]) / "figures"
    for color in ("condition", "replicate"):
        save_embedding(
            result,
            basis="X_lsi",
            color=color,
            output_path=figures / f"atac_lsi_{color}.png",
            title=f"ATAC LSI by {color}",
            axis_prefix="LSI",
        )

    output = Path(config["paths"]["processed"]) / "atac_processed.h5ad"
    output.parent.mkdir(parents=True, exist_ok=True)
    result.write_h5ad(output, compression="lzf")
    LOGGER.info(
        "Removed %d cells; harmonized %d called peaks to %d consensus peaks; "
        "filtered %d consensus peaks",
        cells_before - result.n_obs,
        features_before,
        result.uns["preprocessing"]["consensus_features"],
        result.uns["preprocessing"]["consensus_features"] - result.n_vars,
    )
    LOGGER.info(
        "Output shape: %s; LSI dimensions: %d (LSI1 retained)",
        result.shape,
        result.obsm["X_lsi"].shape[1],
    )
    LOGGER.info("HIC2/NTC retention:\n%s", retention.to_string(index=False))
    LOGGER.info("Output: %s", output)


if __name__ == "__main__":
    run_stage("03_preprocess_atac", main)
