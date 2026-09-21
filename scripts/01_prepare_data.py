"""Create aligned metadata and sparse RNA/ATAC interim objects."""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path

import anndata as ad
import pandas as pd
from _common import run_stage
from anndata.experimental import concat_on_disk

from multiomic_driver.data.align import (
    build_cell_metadata,
    classify_guide_calls,
    metadata_for_modality,
)
from multiomic_driver.data.discovery import discover_experiments
from multiomic_driver.data.loader import (
    load_guide_call_matrix,
    read_barcodes,
    read_features,
    read_matrix_market_header,
    stream_multiome_to_h5ad,
    validate_matrix_dimensions,
)
from multiomic_driver.utils.io import require_path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _validate_shared_arc_files(
    rna_matrix: Path,
    atac_matrix: Path,
    rna_barcodes: pd.Index,
    atac_barcodes: pd.Index,
) -> None:
    if not rna_barcodes.equals(atac_barcodes):
        raise ValueError(
            f"RNA and ATAC barcode order differs: {rna_matrix.name} / {atac_matrix.name}"
        )
    if read_matrix_market_header(rna_matrix) != read_matrix_market_header(atac_matrix):
        raise ValueError("RNA and ATAC Matrix Market headers differ")
    if rna_matrix.stat().st_size != atac_matrix.stat().st_size:
        raise ValueError("RNA and ATAC shared ARC matrix sizes differ")
    if _sha256(rna_matrix) != _sha256(atac_matrix):
        raise ValueError("RNA and ATAC shared ARC matrix contents differ")


def main(config: dict) -> None:
    logger = logging.getLogger("01_prepare_data")
    raw = require_path(config["paths"]["raw"], "raw GSE288996 directory")
    archive = require_path(config["paths"]["guide_calls"], "official guide calls")
    expected_archive_hash = config["guide_calls"]["sha256"].lower()
    observed_archive_hash = _sha256(archive)
    if observed_archive_hash != expected_archive_hash:
        raise ValueError(
            "Official guide-call archive checksum mismatch: "
            f"expected {expected_archive_hash}, observed {observed_archive_hash}"
        )
    experiments = discover_experiments(raw)
    interim = Path(config["paths"]["interim"])
    parts = interim / "_ingestion_parts"
    interim.mkdir(parents=True, exist_ok=True)
    parts.mkdir(parents=True, exist_ok=True)

    all_metadata: list[pd.DataFrame] = []
    rna_parts: list[Path] = []
    atac_parts: list[Path] = []
    for experiment in experiments:
        logger.info(
            "Preparing %s replicate %s", experiment.condition, experiment.replicate
        )
        rna_barcodes = read_barcodes(experiment.rna_barcodes)
        atac_barcodes = read_barcodes(experiment.atac_barcodes)
        features = read_features(experiment.features)
        validate_matrix_dimensions(experiment.rna_matrix, rna_barcodes, features)
        validate_matrix_dimensions(experiment.atac_matrix, atac_barcodes, features)
        _validate_shared_arc_files(
            experiment.rna_matrix,
            experiment.atac_matrix,
            rna_barcodes,
            atac_barcodes,
        )

        calls = load_guide_call_matrix(archive, experiment.group_key)
        assignments = classify_guide_calls(calls)
        metadata = build_cell_metadata(
            rna_barcodes,
            atac_barcodes,
            assignments,
            condition=experiment.condition,
            replicate=experiment.replicate,
            rna_sample_id=experiment.rna_sample_id,
            atac_sample_id=experiment.atac_sample_id,
            guide_sample_id=experiment.guide_sample_id,
        )
        all_metadata.append(metadata)
        obs = metadata_for_modality(metadata, rna_barcodes, "multiome")
        rna_part = parts / f"{experiment.group_key}_rna.h5ad"
        atac_part = parts / f"{experiment.group_key}_atac.h5ad"
        rna_nnz, atac_nnz = stream_multiome_to_h5ad(
            experiment.rna_matrix,
            experiment.rna_barcodes,
            experiment.features,
            obs,
            rna_part,
            atac_part,
        )
        logger.info(
            "Wrote sample parts: %s RNA entries, %s ATAC entries",
            rna_nnz,
            atac_nnz,
        )
        rna_parts.append(rna_part)
        atac_parts.append(atac_part)

    cell_metadata = pd.concat(all_metadata, ignore_index=True)
    if cell_metadata["cell_id"].duplicated().any():
        raise ValueError("Global cell identifiers are not unique")
    cell_metadata.to_csv(interim / "cell_metadata.csv", index=False)
    logger.info("Wrote %s cell metadata rows", len(cell_metadata))

    rna_output = interim / "rna_merged.h5ad"
    atac_output = interim / "atac_merged.h5ad"
    for output in (rna_output, atac_output):
        if output.exists():
            output.unlink()
    concat_on_disk(
        rna_parts,
        rna_output,
        join="inner",
        merge="same",
        max_loaded_elems=10_000_000,
    )
    concat_on_disk(
        atac_parts,
        atac_output,
        join="outer",
        merge="first",
        fill_value=0,
        max_loaded_elems=10_000_000,
    )
    logger.info("Merged RNA and ATAC H5AD outputs")

    rna = ad.read_h5ad(rna_output, backed="r")
    try:
        if rna.n_obs != len(cell_metadata):
            raise ValueError("Merged RNA cell count does not match cell metadata")
    finally:
        rna.file.close()
    atac = ad.read_h5ad(atac_output, backed="r")
    try:
        if atac.n_obs != len(cell_metadata):
            raise ValueError("Merged ATAC cell count does not match cell metadata")
    finally:
        atac.file.close()

    for path in [*rna_parts, *atac_parts]:
        path.unlink()
    parts.rmdir()


if __name__ == "__main__":
    run_stage("01_prepare_data", main)
