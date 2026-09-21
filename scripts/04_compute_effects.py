"""Compute replicate-aware background-corrected modality effects."""

import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from _common import run_stage

from multiomic_driver.effects.replicate_aware import (
    compute_modality_effects,
    validate_eligibility_for_effects,
)
from multiomic_driver.preprocessing.qc import validate_metadata
from multiomic_driver.utils.io import require_path

LOGGER = logging.getLogger("04_compute_effects")


def main(config: dict) -> None:
    processed = Path(config["paths"]["processed"])
    tables = Path(config["paths"]["results"]) / "tables"
    rna_path = require_path(processed / "rna_processed.h5ad")
    atac_path = require_path(processed / "atac_processed.h5ad")
    eligibility_path = require_path(tables / "candidate_phase4_eligibility.csv")
    eligibility = pd.read_csv(eligibility_path)
    minimum = int(config["phase4"]["min_cells_per_group"])
    if not bool(config["phase4"]["require_all_groups"]):
        raise ValueError("Primary Phase 4 effects require all four groups")
    if config["pseudobulk"]["method"] != "mean":
        raise ValueError("Phase 4 latent-space pseudobulk currently requires mean")
    validate_eligibility_for_effects(
        eligibility, min_cells_per_group=minimum
    )
    candidates = eligibility.loc[
        eligibility["phase4_eligible"].astype(bool), "target_gene"
    ].tolist()
    if not candidates:
        raise ValueError("No Phase 4-eligible targeting perturbations")

    rna = ad.read_h5ad(rna_path, backed="r")
    atac = ad.read_h5ad(atac_path, backed="r")
    validate_metadata(rna, "RNA")
    validate_metadata(atac, "ATAC")
    if not rna.obs_names.equals(atac.obs_names):
        raise ValueError("Processed RNA and ATAC cell IDs or order do not match")
    LOGGER.info(
        "Computing replicate-specific effects for %d candidates at minimum %d cells",
        len(candidates),
        minimum,
    )
    modality_inputs = {
        "rna": (np.asarray(rna.obsm["X_pca"]), "PC"),
        "atac": (np.asarray(atac.obsm["X_lsi"]), "LSI"),
    }
    manifest = []
    for modality, (matrix, prefix) in modality_inputs.items():
        replicate_effects, summaries = compute_modality_effects(
            matrix,
            rna.obs,
            candidates,
            min_cells_per_group=minimum,
            dimension_prefix=prefix,
        )
        replicate_effects.insert(0, "modality", modality.upper())
        summaries.insert(0, "modality", modality.upper())
        replicate_effects.to_csv(
            tables / f"{modality}_replicate_effects.csv", index=False
        )
        summaries.to_csv(tables / f"{modality}_effect_summary.csv", index=False)
        manifest.append(
            {
                "modality": modality.upper(),
                "candidates": len(summaries),
                "replicate_effect_rows": len(replicate_effects),
                "dimensions": matrix.shape[1],
                "replicates_merged_before_effect_estimation": False,
                "consistency_metric": "cosine_similarity_R1_vs_R2",
            }
        )
        LOGGER.info(
            "%s: wrote %d replicate effects and %d mean summaries",
            modality.upper(),
            len(replicate_effects),
            len(summaries),
        )
    pd.DataFrame(manifest).to_csv(tables / "phase4_effect_manifest.csv", index=False)
    rna.file.close()
    atac.file.close()


if __name__ == "__main__":
    run_stage("04_compute_effects", main)
