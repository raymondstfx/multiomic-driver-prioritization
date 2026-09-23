"""Preprocess the merged RNA AnnData object."""

import logging
from pathlib import Path

import anndata as ad
from _common import run_stage

from multiomic_driver.preprocessing.qc import (
    summary_table,
    update_cross_modality_tables,
    validate_metadata,
)
from multiomic_driver.preprocessing.rna import preprocess_rna
from multiomic_driver.utils.artifacts import figure_dir, table_dir
from multiomic_driver.utils.io import require_path
from multiomic_driver.visualization.plots import save_embedding

LOGGER = logging.getLogger("02_preprocess_rna")


def main(config: dict) -> None:
    source = require_path(Path(config["paths"]["interim"]) / "rna_merged.h5ad")
    LOGGER.info("Input: %s", source)
    result = ad.read_h5ad(source)
    validate_metadata(result, "RNA")
    cells_before, genes_before = result.shape
    LOGGER.info("Input shape: %s", result.shape)
    LOGGER.info("RNA thresholds: %s", config["rna"])
    settings = {**config["rna"], "random_seed": config["random_seed"]}
    result = preprocess_rna(result, **settings, copy=False, preserve_counts=True)
    validate_metadata(result, "RNA")

    tables = table_dir(config, "phase03")
    metrics = {
        "cells_before": cells_before,
        "cells_after": result.n_obs,
        "cells_removed": cells_before - result.n_obs,
        "genes_before": genes_before,
        "genes_after": result.n_vars,
        "genes_removed": genes_before - result.n_vars,
        "mitochondrial_genes_detected": int(result.var["mt"].sum()),
        "median_total_counts": float(result.obs["total_counts"].median()),
        "median_genes_by_counts": float(result.obs["n_genes_by_counts"].median()),
        "selected_hvg": int(result.var["highly_variable"].sum()),
        "pca_dimensions": result.obsm["X_pca"].shape[1],
    }
    if "pct_counts_mt" in result.obs:
        metrics["median_pct_counts_mt"] = float(result.obs["pct_counts_mt"].median())
    summary_table(metrics).to_csv(tables / "rna_qc_summary.csv", index=False)
    retention = update_cross_modality_tables(result, "RNA", tables)
    if (retention["cell_count"] == 0).any():
        raise ValueError(
            "RNA QC lost at least one HIC2 or NTC condition/replicate group"
        )

    figures = figure_dir(config, "phase03")
    for color in ("condition", "replicate"):
        save_embedding(
            result,
            basis="X_pca",
            color=color,
            output_path=figures / f"rna_pca_{color}.png",
            title=f"RNA PCA by {color}",
            axis_prefix="PC",
        )
    output = Path(config["paths"]["processed"]) / "rna_processed.h5ad"
    output.parent.mkdir(parents=True, exist_ok=True)
    result.write_h5ad(output, compression="lzf")
    LOGGER.info(
        "Removed %d cells and %d genes",
        cells_before - result.n_obs,
        genes_before - result.n_vars,
    )
    LOGGER.info(
        "Output shape: %s; PCA dimensions: %d",
        result.shape,
        result.obsm["X_pca"].shape[1],
    )
    LOGGER.info("HIC2/NTC retention:\n%s", retention.to_string(index=False))
    LOGGER.info("Output: %s", output)


if __name__ == "__main__":
    run_stage("02_preprocess_rna", main)
