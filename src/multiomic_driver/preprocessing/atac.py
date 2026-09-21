"""Sparse TF-IDF and LSI preprocessing for ATAC matrices."""

from __future__ import annotations

import numpy as np
from scipy import sparse
from sklearn.decomposition import TruncatedSVD


def tfidf_transform(matrix: sparse.spmatrix) -> sparse.csr_matrix:
    """Apply sparse term-frequency inverse-document-frequency scaling."""
    x = sparse.csr_matrix(matrix, dtype=float)
    row_sums = np.asarray(x.sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1.0
    tf = sparse.diags(1.0 / row_sums) @ x
    document_frequency = np.asarray((x > 0).sum(axis=0)).ravel()
    idf = np.log1p(x.shape[0] / (1.0 + document_frequency))
    return (tf @ sparse.diags(idf)).tocsr()


def preprocess_atac(
    matrix: sparse.spmatrix,
    *,
    min_cells_per_peak: int = 3,
    n_lsi: int = 30,
    random_seed: int = 42,
) -> tuple[sparse.csr_matrix, np.ndarray]:
    """Filter ATAC peaks, apply TF-IDF, and calculate deterministic LSI."""
    x = sparse.csr_matrix(matrix)
    keep = np.asarray((x > 0).sum(axis=0)).ravel() >= min_cells_per_peak
    filtered = x[:, keep]
    if filtered.shape[1] == 0:
        raise ValueError("No ATAC peaks remain after filtering")
    tfidf = tfidf_transform(filtered)
    max_components = min(tfidf.shape) - 1
    if max_components < 1:
        raise ValueError("ATAC matrix is too small for LSI")
    lsi = TruncatedSVD(
        n_components=min(n_lsi, max_components), random_state=random_seed
    ).fit_transform(tfidf)
    return tfidf, lsi

