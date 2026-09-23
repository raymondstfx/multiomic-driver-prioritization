"""General project visualisations."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from anndata import AnnData


def save_mitochondrial_qc_figures(
    metadata,
    output_dir: str | Path,
) -> None:
    """Save overall, condition, and condition-by-replicate mt% diagnostics."""
    import matplotlib.pyplot as plt

    required = {"pct_counts_mt", "condition", "replicate"}
    missing = required - set(metadata.columns)
    if missing:
        raise ValueError(f"RNA metadata is missing: {sorted(missing)}")
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    values = metadata["pct_counts_mt"].to_numpy(dtype=float)

    fig, axis = plt.subplots(figsize=(7, 5))
    axis.hist(values, bins=80, color="#4c78a8", alpha=0.85)
    axis.axvline(np.median(values), color="#e45756", linestyle="--", label="median")
    axis.set(
        xlabel="Mitochondrial counts (%)",
        ylabel="Cells",
        title="RNA mitochondrial fraction",
    )
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(directory / "rna_mt_distribution.png", dpi=300)
    plt.close(fig)

    conditions = ["DMSO", "Dasatinib"]
    condition_values = [
        metadata.loc[metadata["condition"].astype(str).eq(label), "pct_counts_mt"]
        for label in conditions
    ]
    fig, axis = plt.subplots(figsize=(7, 5))
    axis.violinplot(condition_values, showmedians=True, showextrema=False)
    axis.set_xticks([1, 2], conditions)
    axis.set(
        xlabel="Condition",
        ylabel="Mitochondrial counts (%)",
        title="RNA mitochondrial fraction by condition",
    )
    fig.tight_layout()
    fig.savefig(directory / "rna_mt_by_condition.png", dpi=300)
    plt.close(fig)

    groups = [
        (condition, replicate) for condition in conditions for replicate in (1, 2)
    ]
    group_values = [
        metadata.loc[
            metadata["condition"].astype(str).eq(condition)
            & metadata["replicate"].astype(int).eq(replicate),
            "pct_counts_mt",
        ]
        for condition, replicate in groups
    ]
    labels = [f"{condition}\nR{replicate}" for condition, replicate in groups]
    fig, axis = plt.subplots(figsize=(8, 5))
    axis.boxplot(group_values, tick_labels=labels, showfliers=False)
    axis.set(
        xlabel="Condition × replicate",
        ylabel="Mitochondrial counts (%)",
        title="RNA mitochondrial fraction by experimental group",
    )
    fig.tight_layout()
    fig.savefig(directory / "rna_mt_by_group.png", dpi=300)
    plt.close(fig)


def save_embedding(
    adata: AnnData,
    *,
    basis: str,
    color: str,
    output_path: str | Path,
    title: str,
    axis_prefix: str,
) -> None:
    """Save a compact categorical plot of the first two latent dimensions."""
    import matplotlib.pyplot as plt

    if basis not in adata.obsm:
        raise ValueError(f"Missing embedding: adata.obsm[{basis!r}]")
    if color not in adata.obs:
        raise ValueError(f"Missing observation column: {color}")
    embedding = np.asarray(adata.obsm[basis])
    if embedding.ndim != 2 or embedding.shape[1] < 2:
        raise ValueError(f"{basis} must contain at least two dimensions")

    categories = adata.obs[color].astype(str)
    labels = sorted(categories.unique())
    palette = plt.get_cmap("tab10")
    fig, axis = plt.subplots(figsize=(7, 5.5))
    for index, label in enumerate(labels):
        selected = categories.eq(label).to_numpy()
        axis.scatter(
            embedding[selected, 0],
            embedding[selected, 1],
            s=3,
            alpha=0.45,
            linewidths=0,
            rasterized=True,
            label=label,
            color=palette(index % 10),
        )
    axis.set_xlabel(f"{axis_prefix}1")
    axis.set_ylabel(f"{axis_prefix}2")
    axis.set_title(title)
    axis.legend(title=color, markerscale=3, frameon=False)
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_experimental_design(output_path: str | Path) -> None:
    """Save a compact schematic of the analysis design."""
    import matplotlib.pyplot as plt

    fig, axis = plt.subplots(figsize=(8, 5))
    axis.axis("off")
    text = (
        "DMSO: perturbation + non-targeting\n"
        "Dasatinib: perturbation + non-targeting\n\n"
        "RNA + ATAC\n↓\nBackground correction\n↓\n"
        "Cross-modal scoring\n↓\nCandidate ranking"
    )
    axis.text(0.5, 0.5, text, ha="center", va="center", fontsize=13)
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
