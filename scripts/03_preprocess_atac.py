"""Preprocess the merged ATAC sparse matrix."""

from pathlib import Path

from _common import run_stage
from scipy import sparse

from multiomic_driver.preprocessing.atac import preprocess_atac
from multiomic_driver.utils.io import require_path


def main(config: dict) -> None:
    source = require_path(Path(config["paths"]["interim"]) / "atac_merged.npz")
    settings = {**config["atac"], "random_seed": config["random_seed"]}
    tfidf, lsi = preprocess_atac(sparse.load_npz(source), **settings)
    output = Path(config["paths"]["processed"])
    output.mkdir(parents=True, exist_ok=True)
    sparse.save_npz(output / "atac_tfidf.npz", tfidf)
    import numpy as np

    np.save(output / "atac_lsi.npy", lsi)


if __name__ == "__main__":
    run_stage("03_preprocess_atac", main)

