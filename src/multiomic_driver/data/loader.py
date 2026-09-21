"""Load processed GEO matrices without densifying sparse input."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy import io, sparse


def _require_file(path: str | Path) -> Path:
    resolved = Path(path)
    if not resolved.is_file():
        raise FileNotFoundError(f"Required input file does not exist: {resolved}")
    return resolved


def _read_table(path: Path) -> pd.DataFrame:
    compression = "gzip" if path.suffix == ".gz" else "infer"
    return pd.read_csv(path, sep="\t", header=None, compression=compression)


def load_rna_sample(
    matrix_path: str | Path,
    barcode_path: str | Path,
    feature_path: str | Path,
):
    """Load a Matrix Market RNA sample as AnnData with validated dimensions."""
    from anndata import AnnData

    matrix_file = _require_file(matrix_path)
    barcode_file = _require_file(barcode_path)
    feature_file = _require_file(feature_path)
    matrix = io.mmread(matrix_file).tocsr()
    barcodes = _read_table(barcode_file).iloc[:, 0].astype(str)
    features = _read_table(feature_file)
    # 10x Matrix Market convention is features x cells.
    if matrix.shape != (len(features), len(barcodes)):
        raise ValueError(
            "RNA dimensions are inconsistent: matrix is "
            f"{matrix.shape}, features={len(features)}, barcodes={len(barcodes)}"
        )
    var_names = features.iloc[:, 1] if features.shape[1] > 1 else features.iloc[:, 0]
    return AnnData(
        X=matrix.T.tocsr(),
        obs=pd.DataFrame(index=pd.Index(barcodes, name="cell_barcode")),
        var=pd.DataFrame(index=pd.Index(var_names.astype(str), name="feature")),
    )


def load_atac_sample(matrix_path: str | Path, barcode_path: str | Path):
    """Load an ATAC Matrix Market file and barcodes as a sparse matrix tuple."""
    matrix_file = _require_file(matrix_path)
    barcode_file = _require_file(barcode_path)
    matrix = io.mmread(matrix_file).tocsr()
    barcodes = _read_table(barcode_file).iloc[:, 0].astype(str)
    if matrix.shape[1] == len(barcodes):
        matrix = matrix.T.tocsr()
    elif matrix.shape[0] != len(barcodes):
        raise ValueError(
            f"ATAC matrix {matrix.shape} is incompatible with {len(barcodes)} barcodes"
        )
    return matrix, pd.Index(barcodes, name="cell_barcode")


def load_guide_assignments(guide_path: str | Path) -> pd.DataFrame:
    """Load a tabular guide-assignment file; dataset columns are preserved."""
    path = _require_file(guide_path)
    compression = "gzip" if path.suffix == ".gz" else "infer"
    table = pd.read_csv(path, sep=None, engine="python", compression=compression)
    if table.empty:
        raise ValueError(f"Guide assignment file is empty: {path}")
    return table


def ensure_sparse(matrix) -> sparse.csr_matrix:
    """Return a CSR representation without changing values."""
    return sparse.csr_matrix(matrix)

