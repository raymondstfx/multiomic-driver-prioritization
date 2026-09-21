import gzip
from pathlib import Path

import anndata as ad
import pandas as pd
import pytest

from multiomic_driver.data.loader import (
    load_guide_assignments,
    read_barcodes,
    read_features,
    stream_multiome_to_h5ad,
    validate_matrix_dimensions,
)
from multiomic_driver.data.metadata import parse_sample_name


def test_metadata_parsing() -> None:
    result = parse_sample_name("GSM8780564_K562_Dasatinib_RNA_1")
    assert result == {
        "sample_id": "GSM8780564",
        "cell_line": "K562",
        "condition": "Dasatinib",
        "modality": "RNA",
        "replicate": 1,
    }


def test_missing_file_error() -> None:
    with pytest.raises(FileNotFoundError, match="Required input file"):
        load_guide_assignments(Path("definitely_missing.tsv"))


def _write_gzip(path: Path, text: str) -> None:
    with gzip.open(path, "wt") as stream:
        stream.write(text)


def test_matrix_dimension_validation_and_streaming(tmp_path: Path) -> None:
    matrix = tmp_path / "matrix.mtx.gz"
    barcodes = tmp_path / "barcodes.tsv.gz"
    features = tmp_path / "features.tsv.gz"
    _write_gzip(barcodes, "CELL-A-1\nCELL-B-1\n")
    _write_gzip(
        features,
        "gene-1\tGENE1\tGene Expression\tchr1\t1\t2\n"
        "gene-2\tGENE2\tGene Expression\tchr1\t3\t4\n"
        "chr1:10-20\tchr1:10-20\tPeaks\tchr1\t10\t20\n",
    )
    _write_gzip(
        matrix,
        "%%MatrixMarket matrix coordinate integer general\n"
        "% synthetic\n"
        "3 2 4\n"
        "1 1 2\n"
        "3 1 1\n"
        "2 2 3\n"
        "3 2 4\n",
    )
    barcode_index = read_barcodes(barcodes)
    feature_table = read_features(features)
    assert validate_matrix_dimensions(matrix, barcode_index, feature_table) == (3, 2, 4)
    obs = pd.DataFrame(index=pd.Index(["sample:CELL-A-1", "sample:CELL-B-1"]))
    rna_path = tmp_path / "rna.h5ad"
    atac_path = tmp_path / "atac.h5ad"
    assert stream_multiome_to_h5ad(
        matrix,
        barcodes,
        features,
        obs,
        rna_path,
        atac_path,
    ) == (2, 2)
    rna = ad.read_h5ad(rna_path)
    atac = ad.read_h5ad(atac_path)
    assert rna.shape == (2, 2)
    assert atac.shape == (2, 1)
    assert rna.X.toarray().tolist() == [[2, 0], [0, 3]]
    assert atac.X.toarray().tolist() == [[1], [4]]


def test_duplicate_barcode_handling(tmp_path: Path) -> None:
    barcodes = tmp_path / "barcodes.tsv.gz"
    _write_gzip(barcodes, "CELL-1\nCELL-1\n")
    with pytest.raises(ValueError, match="Duplicate barcodes"):
        read_barcodes(barcodes)


def test_matrix_dimension_mismatch(tmp_path: Path) -> None:
    matrix = tmp_path / "matrix.mtx.gz"
    barcodes = tmp_path / "barcodes.tsv.gz"
    features = tmp_path / "features.tsv.gz"
    _write_gzip(matrix, "%%MatrixMarket matrix coordinate integer general\n2 2 0\n")
    _write_gzip(barcodes, "A-1\n")
    _write_gzip(features, "g1\tG1\tGene Expression\n")
    with pytest.raises(ValueError, match="do not match"):
        validate_matrix_dimensions(matrix, read_barcodes(barcodes), read_features(features))
