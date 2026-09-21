"""Review Phase 3 QC and generate the Phase 4 eligibility gate."""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
from _common import run_stage

from multiomic_driver.evaluation.eligibility import (
    evaluate_candidate_eligibility,
    mitochondrial_qc_by_group,
    mitochondrial_qc_summary,
    ntc_group_coverage,
    readiness_summary,
)
from multiomic_driver.preprocessing.qc import validate_metadata
from multiomic_driver.utils.io import require_path
from multiomic_driver.visualization.plots import save_mitochondrial_qc_figures

LOGGER = logging.getLogger("03b_phase4_readiness")


def main(config: dict) -> None:
    processed = Path(config["paths"]["processed"])
    rna_path = require_path(processed / "rna_processed.h5ad")
    atac_path = require_path(processed / "atac_processed.h5ad")
    LOGGER.info("RNA input: %s", rna_path)
    LOGGER.info("ATAC input: %s", atac_path)
    rna = ad.read_h5ad(rna_path, backed="r")
    atac = ad.read_h5ad(atac_path, backed="r")
    validate_metadata(rna, "RNA")
    validate_metadata(atac, "ATAC")
    if not rna.obs_names.equals(atac.obs_names):
        raise ValueError("Processed RNA and ATAC cell IDs or order do not match")
    LOGGER.info("Validated %d paired processed cells", rna.n_obs)

    settings = config["phase4"]
    minimum = int(settings["min_cells_per_group"])
    require_all = bool(settings["require_all_groups"])
    LOGGER.info(
        "Eligibility settings: min_cells_per_group=%d, require_all_groups=%s",
        minimum,
        require_all,
    )

    tables = Path(config["paths"]["results"]) / "tables"
    figures = Path(config["paths"]["results"]) / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    mt_summary = mitochondrial_qc_summary(rna.obs)
    mt_by_group = mitochondrial_qc_by_group(rna.obs)
    mt_summary.to_csv(tables / "rna_mt_qc_summary.csv", index=False)
    mt_by_group.to_csv(tables / "rna_mt_qc_by_group.csv", index=False)
    save_mitochondrial_qc_figures(rna.obs, figures)

    eligibility = evaluate_candidate_eligibility(
        rna.obs,
        min_cells_per_group=minimum,
        require_all_groups=require_all,
    )
    ntc = ntc_group_coverage(rna.obs, min_cells_per_group=minimum)
    readiness = readiness_summary(eligibility, ntc)
    eligibility.to_csv(tables / "candidate_phase4_eligibility.csv", index=False)
    eligibility.loc[eligibility["phase4_eligible"]].to_csv(
        tables / "phase4_primary_candidates.csv", index=False
    )
    eligibility.loc[~eligibility["phase4_eligible"]].to_csv(
        tables / "phase4_secondary_candidates.csv", index=False
    )
    ntc.to_csv(tables / "ntc_phase4_coverage.csv", index=False)
    readiness.to_csv(tables / "phase4_readiness_summary.csv", index=False)

    if "ZFPM2" in set(eligibility["target_gene"]):
        raise ValueError("ZFPM2 incorrectly appeared as a pooled perturbation")
    hic2 = eligibility.loc[eligibility["target_gene"].eq("HIC2")]
    expected_hic2 = [183, 230, 325, 320]
    count_columns = [
        "dmso_r1_cells",
        "dmso_r2_cells",
        "dasatinib_r1_cells",
        "dasatinib_r2_cells",
    ]
    if len(hic2) != 1 or hic2[count_columns].iloc[0].tolist() != expected_hic2:
        raise ValueError("HIC2 counts changed unexpectedly from Phase 3")
    if not bool(hic2.iloc[0]["phase4_eligible"]):
        raise ValueError("HIC2 is not Phase 4 eligible")
    if not bool(readiness.set_index("metric").at["phase4_ready", "value"]):
        raise ValueError("Phase 4 readiness gate failed; see readiness summary")

    LOGGER.info("Mitochondrial QC by group:\n%s", mt_by_group.to_string(index=False))
    LOGGER.info(
        "Candidate eligibility:\n%s",
        eligibility.to_string(index=False),
    )
    LOGGER.info("NTC coverage:\n%s", ntc.to_string(index=False))
    LOGGER.info("Readiness summary:\n%s", readiness.to_string(index=False))
    rna.file.close()
    atac.file.close()


if __name__ == "__main__":
    run_stage("03b_phase4_readiness", main)
