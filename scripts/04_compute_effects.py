"""Construct pseudobulk units and background-corrected modality effects."""

from pathlib import Path

from _common import run_stage

from multiomic_driver.utils.io import require_path


def main(config: dict) -> None:
    processed = Path(config["paths"]["processed"])
    require_path(processed / "rna_processed.h5ad")
    require_path(processed / "atac_lsi.npy")
    require_path(Path(config["paths"]["interim"]) / "cell_metadata.csv")
    raise NotImplementedError(
        "Effect table construction awaits confirmation of guide and non-targeting labels."
    )


if __name__ == "__main__":
    run_stage("04_compute_effects", main)

