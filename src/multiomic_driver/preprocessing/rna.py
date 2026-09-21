"""RNA preprocessing using Scanpy primitives."""

from __future__ import annotations


def preprocess_rna(
    adata,
    *,
    min_genes: int = 200,
    min_cells_per_gene: int = 3,
    target_sum: float = 10_000,
    n_hvg: int = 3_000,
    n_pcs: int = 30,
):
    """Filter, normalise, log-transform, select HVGs, and calculate PCA."""
    import scanpy as sc

    result = adata.copy()
    sc.pp.filter_cells(result, min_genes=min_genes)
    sc.pp.filter_genes(result, min_cells=min_cells_per_gene)
    sc.pp.normalize_total(result, target_sum=target_sum)
    sc.pp.log1p(result)
    sc.pp.highly_variable_genes(result, n_top_genes=min(n_hvg, result.n_vars))
    sc.pp.pca(result, n_comps=min(n_pcs, result.n_obs - 1, result.n_vars - 1))
    return result

