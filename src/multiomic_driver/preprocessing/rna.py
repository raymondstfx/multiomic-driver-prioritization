"""Sparse-safe RNA QC, normalisation, HVG selection, and PCA."""

from __future__ import annotations

import numpy as np
from anndata import AnnData
from scipy import sparse


def detect_mitochondrial_genes(adata: AnnData) -> np.ndarray:
    """Detect mitochondrial genes from symbols and chromosome annotations."""
    if "feature_name" in adata.var:
        symbols = adata.var["feature_name"].astype(str).str.upper()
    else:
        symbols = adata.var_names.astype(str).str.upper()
    detected = np.asarray(symbols.str.startswith("MT-"), dtype=bool).copy()
    if "chromosome" in adata.var:
        chromosomes = adata.var["chromosome"].astype(str).str.upper()
        detected |= np.asarray(chromosomes.isin(["CHRM", "MT", "M"]), dtype=bool)
    return detected


def add_rna_qc_metrics(adata: AnnData) -> None:
    """Add per-cell and per-gene count metrics without dense conversion."""
    if not sparse.issparse(adata.X):
        raise TypeError("RNA preprocessing requires a sparse count matrix")
    matrix = adata.X.tocsr()
    adata.obs["total_counts"] = np.asarray(matrix.sum(axis=1)).ravel()
    adata.obs["n_genes_by_counts"] = np.diff(matrix.indptr)
    adata.var["n_cells_by_counts"] = np.asarray(matrix.getnnz(axis=0)).ravel()
    mitochondrial = detect_mitochondrial_genes(adata)
    adata.var["mt"] = mitochondrial
    if mitochondrial.any():
        mt_counts = np.asarray(matrix[:, mitochondrial].sum(axis=1)).ravel()
        totals = adata.obs["total_counts"].to_numpy(dtype=float)
        adata.obs["pct_counts_mt"] = np.divide(
            mt_counts * 100.0,
            totals,
            out=np.zeros_like(totals),
            where=totals > 0,
        )


def preprocess_rna(
    adata: AnnData,
    *,
    min_genes: int = 200,
    min_cells_per_gene: int = 3,
    target_sum: float = 10_000,
    n_hvg: int = 3_000,
    n_pcs: int = 30,
    random_seed: int = 42,
    copy: bool = True,
    preserve_counts: bool = True,
) -> AnnData:
    """Run conservative RNA preprocessing and retain all filtered genes.

    Counts are stored in ``layers['counts']`` before normalisation when
    ``preserve_counts`` is true. PCA uses the HVG mask, while the normalised
    expression matrix continues to contain every retained gene.
    """
    import scanpy as sc

    result = adata.copy() if copy else adata
    if not result.obs_names.is_unique:
        raise ValueError("RNA cell IDs must be unique")
    add_rna_qc_metrics(result)
    cell_mask = result.obs["n_genes_by_counts"].to_numpy() >= min_genes
    gene_mask = result.var["n_cells_by_counts"].to_numpy() >= min_cells_per_gene
    if not cell_mask.any():
        raise ValueError("RNA cell filtering removed every cell")
    if not gene_mask.any():
        raise ValueError("RNA gene filtering removed every gene")
    if not cell_mask.all():
        result._inplace_subset_obs(cell_mask)
    if not gene_mask.all():
        result._inplace_subset_var(gene_mask)
    if preserve_counts:
        result.layers["counts"] = result.X.copy()

    sc.pp.normalize_total(result, target_sum=target_sum)
    sc.pp.log1p(result)
    if not np.isfinite(result.X.data).all():
        raise ValueError("RNA log-normalisation produced non-finite values")
    requested_hvg = min(n_hvg, result.n_vars)
    sc.pp.highly_variable_genes(
        result,
        n_top_genes=requested_hvg,
        flavor="seurat",
        inplace=True,
    )
    actual_hvg = int(result.var["highly_variable"].sum())
    if actual_hvg == 0:
        raise ValueError("No highly variable genes were selected")
    requested_pcs = min(n_pcs, result.n_obs - 1, actual_hvg - 1)
    if requested_pcs < 1:
        raise ValueError("RNA data are too small for PCA")
    sc.pp.pca(
        result,
        n_comps=requested_pcs,
        mask_var="highly_variable",
        random_state=random_seed,
    )
    result.uns["preprocessing"] = {
        "min_genes": min_genes,
        "min_cells_per_gene": min_cells_per_gene,
        "target_sum": target_sum,
        "requested_hvg": n_hvg,
        "selected_hvg": actual_hvg,
        "requested_pcs": n_pcs,
        "computed_pcs": requested_pcs,
        "random_seed": random_seed,
        "counts_preserved": preserve_counts,
    }
    return result
