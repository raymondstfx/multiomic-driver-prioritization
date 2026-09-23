"""Synthetic tests for Phase 3 preprocessing and metadata retention."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from multiomic_driver.preprocessing.atac import preprocess_atac, tfidf_transform
from multiomic_driver.preprocessing.qc import validate_metadata
from multiomic_driver.preprocessing.rna import preprocess_rna


def _metadata(n_cells: int) -> pd.DataFrame:
    index = [f"DMSO:1:barcode-{index}" for index in range(n_cells)]
    return pd.DataFrame(
        {
            "cell_barcode": [f"barcode-{index}" for index in range(n_cells)],
            "condition": ["DMSO"] * (n_cells // 2)
            + ["Dasatinib"] * (n_cells - n_cells // 2),
            "replicate": [1, 2] * (n_cells // 2),
            "guide_assignment_status": ["singlet"] * n_cells,
            "guide_present": [True] * n_cells,
            "guide": ["HIC2-1"] * n_cells,
            "target_gene": ["HIC2"] * n_cells,
            "is_non_targeting": [False] * n_cells,
        },
        index=index,
    )


def test_rna_preprocessing_filters_normalizes_and_preserves_metadata() -> None:
    counts = sparse.csr_matrix(
        [
            [1, 0, 0, 0, 0, 0],
            [8, 3, 1, 0, 2, 1],
            [2, 7, 1, 3, 0, 2],
            [3, 1, 8, 2, 1, 0],
            [4, 2, 1, 9, 3, 1],
            [1, 4, 3, 2, 8, 2],
        ],
        dtype=np.int32,
    )
    data = ad.AnnData(
        counts,
        obs=_metadata(6),
        var=pd.DataFrame(
            {"feature_name": ["MT-A", "B", "C", "D", "E", "F"]},
            index=[f"gene-{index}" for index in range(6)],
        ),
    )
    result = preprocess_rna(
        data,
        min_genes=2,
        min_cells_per_gene=1,
        target_sum=100,
        n_hvg=4,
        n_pcs=2,
        random_seed=7,
    )

    assert result.n_obs == 5
    np.testing.assert_allclose(
        np.asarray(np.expm1(result.X).sum(axis=1)).ravel(), 100, rtol=1e-5
    )
    assert np.isfinite(result.X.data).all()
    assert result.obsm["X_pca"].shape == (5, 2)
    assert "counts" in result.layers
    assert result.obs_names.tolist() == data.obs_names[1:].tolist()
    for column in ("condition", "replicate", "guide", "target_gene"):
        assert result.obs[column].tolist() == data.obs.iloc[1:][column].tolist()

    repeated = preprocess_rna(
        data,
        min_genes=2,
        min_cells_per_gene=1,
        target_sum=100,
        n_hvg=4,
        n_pcs=2,
        random_seed=7,
    )
    np.testing.assert_array_equal(
        result.var["highly_variable"], repeated.var["highly_variable"]
    )


def test_tfidf_formula_is_sparse_and_explicit() -> None:
    counts = sparse.csr_matrix([[2, 0, 1], [0, 3, 1]], dtype=np.int32)
    result = tfidf_transform(counts)
    expected = np.array(
        [
            [2 / 3 * np.log1p(2 / 2), 0, 1 / 3 * np.log1p(2 / 3)],
            [0, 3 / 4 * np.log1p(2 / 2), 1 / 4 * np.log1p(2 / 3)],
        ]
    )
    assert sparse.isspmatrix_csr(result)
    assert result.dtype == np.float32
    np.testing.assert_allclose(result.toarray(), expected, rtol=1e-6)


def test_atac_preprocessing_filters_peaks_and_preserves_metadata() -> None:
    counts = sparse.csr_matrix(
        [
            [1, 0, 1, 0],
            [1, 2, 0, 0],
            [0, 1, 2, 0],
            [1, 0, 1, 0],
            [0, 1, 1, 0],
            [1, 1, 0, 1],
        ],
        dtype=np.int32,
    )
    data = ad.AnnData(
        counts,
        obs=_metadata(6),
        var=pd.DataFrame(index=["chr1:0-10", "chr1:5-15", "chr1:20-30", "chr1:40-50"]),
    )
    result = preprocess_atac(
        data,
        min_cells_per_peak=2,
        n_lsi=2,
        random_seed=7,
    )

    assert result.shape == (6, 2)
    assert sparse.isspmatrix_csr(result.X)
    assert result.obsm["X_lsi"].shape == (6, 2)
    assert result.uns["preprocessing"]["lsi1_retained"] is True
    assert result.obs_names.tolist() == data.obs_names.tolist()
    assert result.obs["guide"].tolist() == data.obs["guide"].tolist()
    validate_metadata(result, "ATAC")
