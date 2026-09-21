"""Load and stream GSE288996 processed data without dense conversion."""

from __future__ import annotations

import gzip
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import rdata
from anndata import AnnData
from scipy import io, sparse

FEATURE_COLUMNS = [
    "feature_id",
    "feature_name",
    "feature_type",
    "chromosome",
    "start",
    "end",
]


def _require_file(path: str | Path) -> Path:
    resolved = Path(path)
    if not resolved.is_file():
        raise FileNotFoundError(f"Required input file does not exist: {resolved}")
    return resolved


def _read_table(path: Path, **kwargs) -> pd.DataFrame:
    compression = "gzip" if path.suffix == ".gz" else "infer"
    return pd.read_csv(path, compression=compression, **kwargs)


def read_barcodes(path: str | Path) -> pd.Index:
    """Read a one-column barcode file, preserving the 10x suffix."""
    barcode_path = _require_file(path)
    values = _read_table(barcode_path, sep="\t", header=None).iloc[:, 0].astype(str)
    if values.duplicated().any():
        duplicates = values[values.duplicated()].head().tolist()
        raise ValueError(f"Duplicate barcodes in {barcode_path}: {duplicates}")
    return pd.Index(values, name="cell_barcode")


def read_features(path: str | Path) -> pd.DataFrame:
    """Read the six-column Cell Ranger ARC feature annotation."""
    feature_path = _require_file(path)
    features = _read_table(feature_path, sep="\t", header=None)
    if features.shape[1] < 3:
        raise ValueError(f"Feature file has fewer than three columns: {feature_path}")
    features.columns = FEATURE_COLUMNS[: features.shape[1]]
    features["feature_id"] = features["feature_id"].astype(str)
    features["feature_name"] = features["feature_name"].astype(str)
    features["feature_type"] = features["feature_type"].astype(str)
    if features["feature_id"].duplicated().any():
        raise ValueError(f"Feature identifiers are not unique: {feature_path}")
    return features


def read_matrix_market_header(path: str | Path) -> tuple[int, int, int]:
    """Read `(features, cells, nonzero entries)` without loading a matrix."""
    matrix_path = _require_file(path)
    opener = gzip.open if matrix_path.suffix == ".gz" else open
    with opener(matrix_path, "rt") as stream:
        banner = stream.readline().strip()
        if not banner.startswith("%%MatrixMarket matrix coordinate"):
            raise ValueError(f"Unsupported Matrix Market banner in {matrix_path}")
        for line in stream:
            if not line.startswith("%"):
                parts = line.split()
                if len(parts) != 3:
                    raise ValueError(f"Invalid Matrix Market dimensions: {line!r}")
                return tuple(int(value) for value in parts)  # type: ignore[return-value]
    raise ValueError(f"Matrix Market dimensions are missing: {matrix_path}")


def validate_matrix_dimensions(
    matrix_path: str | Path,
    barcodes: pd.Index,
    features: pd.DataFrame,
) -> tuple[int, int, int]:
    """Validate the observed Cell Ranger `features x cells` orientation."""
    shape = read_matrix_market_header(matrix_path)
    if shape[0] != len(features) or shape[1] != len(barcodes):
        raise ValueError(
            f"Matrix dimensions {shape[:2]} do not match "
            f"features={len(features)}, barcodes={len(barcodes)}"
        )
    return shape


def load_rna_sample(
    matrix_path: str | Path,
    barcode_path: str | Path,
    feature_path: str | Path,
) -> AnnData:
    """Load only Gene Expression rows from a small multiome matrix."""
    matrix_file = _require_file(matrix_path)
    barcodes = read_barcodes(barcode_path)
    features = read_features(feature_path)
    validate_matrix_dimensions(matrix_file, barcodes, features)
    matrix = io.mmread(matrix_file).tocsr()
    selected = features["feature_type"].eq("Gene Expression").to_numpy()
    matrix = matrix[selected].T.tocsr()
    selected_features = features.loc[selected].copy()
    selected_features.index = pd.Index(selected_features.pop("feature_id"), name="feature_id")
    return AnnData(X=matrix, obs=pd.DataFrame(index=barcodes), var=selected_features)


def load_atac_sample(
    matrix_path: str | Path,
    barcode_path: str | Path,
    feature_path: str | Path | None = None,
):
    """Load Peaks rows from a small multiome matrix as cells by peaks."""
    matrix_file = _require_file(matrix_path)
    barcodes = read_barcodes(barcode_path)
    matrix = io.mmread(matrix_file).tocsr()
    if feature_path is None:
        if matrix.shape[1] == len(barcodes):
            return matrix.T.tocsr(), barcodes
        if matrix.shape[0] == len(barcodes):
            return matrix, barcodes
        raise ValueError(
            f"ATAC matrix {matrix.shape} is incompatible with {len(barcodes)} barcodes"
        )
    features = read_features(feature_path)
    validate_matrix_dimensions(matrix_file, barcodes, features)
    selected = features["feature_type"].eq("Peaks").to_numpy()
    return matrix[selected].T.tocsr(), barcodes


def load_guide_count_summary(path: str | Path) -> pd.DataFrame:
    """Load GEO's aggregate two-column guide read-count summary."""
    guide_path = _require_file(path)
    table = _read_table(
        guide_path,
        sep=r"\s+",
        header=None,
        names=["read_count", "guide_id"],
    )
    if table.empty or table["guide_id"].isna().any():
        raise ValueError(f"Invalid aggregate guide summary: {guide_path}")
    table["read_count"] = pd.to_numeric(table["read_count"], errors="raise")
    return table


def load_guide_assignments(guide_path: str | Path) -> pd.DataFrame:
    """Load a generic tabular guide-assignment file for compatibility."""
    path = _require_file(guide_path)
    compression = "gzip" if path.suffix == ".gz" else "infer"
    table = pd.read_csv(path, sep=None, engine="python", compression=compression)
    if table.empty:
        raise ValueError(f"Guide assignment file is empty: {path}")
    return table


def load_guide_call_matrix(
    archive_path: str | Path, group_key: str
) -> pd.DataFrame:
    """Read an official boolean guide-call RDS directly from the author archive."""
    path = _require_file(archive_path)
    suffix = f"/{group_key}_calls.rds"
    with ZipFile(path) as archive:
        matches = [name for name in archive.namelist() if name.endswith(suffix)]
        if len(matches) != 1:
            raise ValueError(
                f"Expected one {group_key} calls RDS in {path}, found {len(matches)}"
            )
        with archive.open(matches[0]) as stream:
            calls = rdata.read_rds(stream)
    if not isinstance(calls, pd.DataFrame):
        raise TypeError(f"{group_key} guide calls are not a data frame")
    calls.index = pd.Index(calls.index.astype(str), name="cell_barcode")
    calls.columns = pd.Index(calls.columns.astype(str), name="guide_id")
    if calls.index.duplicated().any():
        raise ValueError(f"Duplicate cell barcodes in {group_key} guide calls")
    return calls.fillna(False).astype(bool)


def stream_multiome_to_h5ad(
    matrix_path: str | Path,
    barcode_path: str | Path,
    feature_path: str | Path,
    obs: pd.DataFrame,
    rna_output: str | Path,
    atac_output: str | Path,
) -> tuple[int, int]:
    """Split an ARC matrix into sparse RNA and ATAC H5AD files.

    SciPy's parallel C++ Matrix Market reader is used one sample at a time.
    Values are downcast to 32-bit integers before CSR conversion to keep the
    peak working set practical on a 16 GiB workstation.
    """
    source = _require_file(matrix_path)
    barcodes = read_barcodes(barcode_path)
    features = read_features(feature_path)
    _n_features, n_cells, expected_nnz = validate_matrix_dimensions(
        source, barcodes, features
    )
    if len(obs) != n_cells:
        raise ValueError(f"obs has {len(obs)} rows but matrix has {n_cells} cells")

    rna_selected = features["feature_type"].eq("Gene Expression").to_numpy()
    atac_selected = features["feature_type"].eq("Peaks").to_numpy()
    if not rna_selected.any() or not atac_selected.any():
        raise ValueError("Expected both Gene Expression and Peaks feature rows")
    if np.any(rna_selected & atac_selected) or not np.all(rna_selected | atac_selected):
        unknown = sorted(set(features.loc[~(rna_selected | atac_selected), "feature_type"]))
        raise ValueError(f"Unsupported feature types: {unknown}")

    rna_path = Path(rna_output)
    atac_path = Path(atac_output)
    rna_path.parent.mkdir(parents=True, exist_ok=True)
    atac_path.parent.mkdir(parents=True, exist_ok=True)
    matrix = io.mmread(source)
    if matrix.nnz != expected_nnz:
        raise ValueError(
            f"Matrix declares {expected_nnz} entries but {matrix.nnz} were read"
        )
    matrix.data = matrix.data.astype(np.int32, copy=False)
    cells_by_features = matrix.T.tocsr()
    del matrix

    rna_matrix = cells_by_features[:, rna_selected].tocsr()
    rna_features = features.loc[rna_selected].copy()
    rna_features.index = pd.Index(
        rna_features.pop("feature_id"), name="feature_id"
    )
    AnnData(X=rna_matrix, obs=obs.copy(), var=rna_features).write_h5ad(
        rna_path, compression="lzf"
    )
    rna_nnz = rna_matrix.nnz
    del rna_matrix

    atac_matrix = cells_by_features[:, atac_selected].tocsr()
    atac_features = features.loc[atac_selected].copy()
    atac_features.index = pd.Index(
        atac_features.pop("feature_id"), name="feature_id"
    )
    AnnData(X=atac_matrix, obs=obs.copy(), var=atac_features).write_h5ad(
        atac_path, compression="lzf"
    )
    atac_nnz = atac_matrix.nnz
    if rna_nnz + atac_nnz != expected_nnz:
        raise ValueError("Feature-type split did not preserve all matrix entries")
    return int(rna_nnz), int(atac_nnz)


def ensure_sparse(matrix) -> sparse.csr_matrix:
    """Return a CSR representation without changing values."""
    return sparse.csr_matrix(matrix)
