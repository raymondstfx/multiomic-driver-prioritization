"""Preprocess the merged RNA AnnData object."""

from pathlib import Path

import anndata as ad
from _common import run_stage

from multiomic_driver.preprocessing.rna import preprocess_rna
from multiomic_driver.utils.io import require_path


def main(config: dict) -> None:
    source = require_path(Path(config["paths"]["interim"]) / "rna_merged.h5ad")
    result = preprocess_rna(ad.read_h5ad(source), **config["rna"])
    output = Path(config["paths"]["processed"]) / "rna_processed.h5ad"
    output.parent.mkdir(parents=True, exist_ok=True)
    result.write_h5ad(output)


if __name__ == "__main__":
    run_stage("02_preprocess_rna", main)

