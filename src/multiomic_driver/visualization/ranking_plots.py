"""Plots for modality scores and candidate rankings."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import pandas as pd


def plot_rna_vs_atac(
    ranking: pd.DataFrame,
    output_path: str | Path,
    highlight: Iterable[str] = ("HIC2", "ZFPM2"),
) -> None:
    """Save an RNA-versus-ATAC score scatter plot."""
    import matplotlib.pyplot as plt

    fig, axis = plt.subplots(figsize=(6, 5))
    axis.scatter(ranking["rna_score"], ranking["atac_score"], alpha=0.7)
    for gene in highlight:
        rows = ranking[ranking["candidate"] == gene]
        if not rows.empty:
            row = rows.iloc[0]
            axis.annotate(gene, (row["rna_score"], row["atac_score"]))
    axis.set(xlabel="RNA effect score", ylabel="ATAC effect score")
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)

