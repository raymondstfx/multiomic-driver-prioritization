"""Memory-aware sparse TF-IDF and LSI preprocessing for ATAC data."""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse
from sklearn.decomposition import TruncatedSVD


def collapse_overlapping_peaks(adata: AnnData) -> AnnData:
    """Collapse overlapping called peaks into comparable consensus intervals.

    The input matrices use independently called peak coordinates per sample.
    This operation harmonises those existing intervals without calling new
    peaks or requiring fragment files.
    """
    coordinates = adata.var_names.to_series().str.extract(
        r"^(?P<chromosome>.+):(?P<start>\d+)-(?P<end>\d+)$"
    )
    if coordinates.isna().any(axis=None):
        raise ValueError("ATAC feature IDs must use chromosome:start-end coordinates")
    coordinates[["start", "end"]] = coordinates[["start", "end"]].astype(
        np.int64
    )
    coordinates["original_index"] = np.arange(adata.n_vars)
    ordered = coordinates.sort_values(["chromosome", "start", "end"])

    assignments = np.empty(adata.n_vars, dtype=np.int32)
    regions: list[tuple[str, int, int]] = []
    region_index = -1
    for chromosome, group in ordered.groupby("chromosome", sort=False):
        current_start: int | None = None
        current_end: int | None = None
        for row in group.itertuples():
            if current_start is None or row.start > current_end:
                if current_start is not None:
                    regions.append((chromosome, current_start, current_end))
                region_index += 1
                current_start, current_end = row.start, row.end
            else:
                current_end = max(current_end, row.end)
            assignments[row.original_index] = region_index
        if current_start is not None:
            regions.append((chromosome, current_start, current_end))

    matrix = adata.X.tocsr(copy=False)
    for start in range(0, matrix.nnz, 10_000_000):
        stop = min(start + 10_000_000, matrix.nnz)
        matrix.indices[start:stop] = assignments[matrix.indices[start:stop]]
    matrix = sparse.csr_matrix(
        (matrix.data, matrix.indices, matrix.indptr),
        shape=(adata.n_obs, len(regions)),
        copy=False,
    )
    matrix.sum_duplicates()
    matrix.sort_indices()

    var = pd.DataFrame(regions, columns=["chromosome", "start", "end"])
    var.index = pd.Index(
        [f"{chromosome}:{start}-{end}" for chromosome, start, end in regions],
        name=adata.var_names.name,
    )
    var.insert(0, "feature_type", "Peaks")
    var.insert(0, "feature_name", var.index)
    return AnnData(
        X=matrix,
        obs=adata.obs.copy(),
        var=var,
        uns=adata.uns.copy(),
        obsm=adata.obsm.copy(),
    )


def tfidf_transform(
    matrix: sparse.spmatrix,
    *,
    copy: bool = True,
    value_chunk_size: int = 10_000_000,
) -> sparse.csr_matrix:
    """Apply sparse TF-IDF using ``TF_ij * log(1 + N/(1 + DF_j))``.

    Term frequency is each peak count divided by the total counts in that
    cell. ``N`` is the number of cells and ``DF_j`` is the number of cells
    with a non-zero value for peak ``j``. Scaling is chunked to avoid a dense
    temporary array proportional to the hundreds of millions of non-zero
    entries in the real dataset.
    """
    if not sparse.issparse(matrix):
        raise TypeError("ATAC TF-IDF requires a sparse matrix")
    x = matrix.copy().tocsr() if copy else matrix.tocsr(copy=False)
    if x.dtype != np.float32:
        x.data = x.data.astype(np.float32)
    row_sums = np.asarray(x.sum(axis=1)).ravel().astype(np.float32)
    inverse_totals = np.divide(
        1.0,
        row_sums,
        out=np.zeros_like(row_sums),
        where=row_sums > 0,
    )
    document_frequency = np.bincount(x.indices, minlength=x.shape[1])
    idf = np.log1p(x.shape[0] / (1.0 + document_frequency)).astype(np.float32)

    row_chunk_size = 2_000
    for first_row in range(0, x.shape[0], row_chunk_size):
        last_row = min(first_row + row_chunk_size, x.shape[0])
        first_value = x.indptr[first_row]
        last_value = x.indptr[last_row]
        repeats = np.diff(x.indptr[first_row : last_row + 1])
        x.data[first_value:last_value] *= np.repeat(
            inverse_totals[first_row:last_row], repeats
        )
    for start in range(0, x.nnz, value_chunk_size):
        stop = min(start + value_chunk_size, x.nnz)
        x.data[start:stop] *= idf[x.indices[start:stop]]
    x.eliminate_zeros()
    return x


def preprocess_atac(
    adata: AnnData,
    *,
    min_cells_per_peak: int = 3,
    n_lsi: int = 30,
    random_seed: int = 42,
    copy: bool = True,
    harmonize_overlapping_peaks: bool = True,
) -> AnnData:
    """Filter peaks, calculate sparse float32 TF-IDF, and retain LSI1."""
    result = adata.copy() if copy else adata
    if not result.obs_names.is_unique:
        raise ValueError("ATAC cell IDs must be unique")
    if not sparse.issparse(result.X):
        raise TypeError("ATAC preprocessing requires a sparse count matrix")
    input_features = result.n_vars
    if harmonize_overlapping_peaks:
        result = collapse_overlapping_peaks(result)
    consensus_features = result.n_vars
    matrix = result.X.tocsr(copy=False)
    feature_cells = np.bincount(matrix.indices, minlength=result.n_vars)
    result.var["n_cells_by_counts"] = feature_cells
    keep = feature_cells >= min_cells_per_peak
    if not keep.any():
        raise ValueError("No ATAC peaks remain after filtering")
    if not keep.all():
        result._inplace_subset_var(keep)
    result.X = tfidf_transform(result.X, copy=False)

    max_components = min(result.shape)
    computed_lsi = min(n_lsi, max_components)
    if computed_lsi < 1:
        raise ValueError("ATAC matrix is too small for LSI")
    model = TruncatedSVD(
        n_components=computed_lsi,
        n_iter=5,
        random_state=random_seed,
    )
    result.obsm["X_lsi"] = model.fit_transform(result.X).astype(np.float32)
    result.uns["lsi"] = {
        "singular_values": model.singular_values_.astype(np.float32),
        "explained_variance": model.explained_variance_.astype(np.float32),
        "explained_variance_ratio": model.explained_variance_ratio_.astype(
            np.float32
        ),
    }
    result.uns["preprocessing"] = {
        "tfidf_formula": "TF_ij * log(1 + N / (1 + DF_j))",
        "input_features": input_features,
        "consensus_features": consensus_features,
        "features_after_filtering": result.n_vars,
        "overlapping_peaks_harmonized": harmonize_overlapping_peaks,
        "min_cells_per_peak": min_cells_per_peak,
        "requested_lsi": n_lsi,
        "computed_lsi": computed_lsi,
        "lsi1_retained": True,
        "random_seed": random_seed,
        "matrix_dtype": "float32",
    }
    return result
