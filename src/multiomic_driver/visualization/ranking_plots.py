"""Diagnostic plots for modality scores and candidate rankings."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def _finish(fig, output_path: str | Path) -> None:
    import matplotlib.pyplot as plt

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _scatter(
    ranking: pd.DataFrame,
    *,
    x: str,
    y: str,
    xlabel: str,
    ylabel: str,
    title: str,
    output_path: str | Path,
    reference_lines: bool = False,
) -> None:
    import matplotlib.pyplot as plt

    fig, axis = plt.subplots(figsize=(6.5, 5.5))
    is_hic2 = ranking["candidate"].eq("HIC2")
    axis.scatter(
        ranking.loc[~is_hic2, x],
        ranking.loc[~is_hic2, y],
        color="#4c78a8",
        alpha=0.8,
        label="Other primary candidates",
    )
    if is_hic2.any():
        row = ranking.loc[is_hic2].iloc[0]
        axis.scatter(row[x], row[y], color="#e45756", s=70, label="HIC2")
        axis.annotate("HIC2", (row[x], row[y]), xytext=(5, 5), textcoords="offset points")
    if reference_lines:
        axis.axhline(0, color="0.7", linewidth=1, linestyle="--")
        axis.axvline(0, color="0.7", linewidth=1, linestyle="--")
    axis.set(xlabel=xlabel, ylabel=ylabel, title=title)
    axis.legend(frameon=False)
    _finish(fig, output_path)


def plot_rna_vs_atac(ranking: pd.DataFrame, output_path: str | Path) -> None:
    """Save a raw RNA-versus-ATAC score scatter plot."""
    _scatter(
        ranking,
        x="rna_score",
        y="atac_score",
        xlabel="RNA mean-effect norm",
        ylabel="ATAC mean-effect norm",
        title="RNA versus ATAC perturbation evidence",
        output_path=output_path,
    )


def plot_rna_vs_atac_zscores(
    ranking: pd.DataFrame, output_path: str | Path
) -> None:
    """Save a standardized cross-modality score scatter plot."""
    _scatter(
        ranking,
        x="rna_z",
        y="atac_z",
        xlabel="RNA evidence (z-score)",
        ylabel="ATAC evidence (z-score)",
        title="Standardized RNA and ATAC evidence",
        output_path=output_path,
        reference_lines=True,
    )


def plot_candidate_ranking(
    ranking: pd.DataFrame, output_path: str | Path
) -> None:
    """Plot multi-omic scores without hiding negative values."""
    import matplotlib.pyplot as plt

    ordered = ranking.sort_values("multiomic_score", ascending=True)
    colors = np.where(ordered["candidate"].eq("HIC2"), "#e45756", "#4c78a8")
    fig, axis = plt.subplots(figsize=(7.5, 6))
    axis.barh(ordered["candidate"], ordered["multiomic_score"], color=colors)
    axis.axvline(0, color="0.4", linewidth=1)
    axis.set(
        xlabel="Multi-omic score (mean modality z-score)",
        ylabel="Candidate",
        title="Primary multi-omic candidate ranking",
    )
    _finish(fig, output_path)


def plot_modality_rank_comparison(
    ranking: pd.DataFrame, output_path: str | Path
) -> None:
    """Compare RNA, ATAC, and combined ranks for every candidate."""
    import matplotlib.pyplot as plt

    ordered = ranking.sort_values("multiomic_rank", ascending=False)
    y = np.arange(len(ordered))
    fig, axis = plt.subplots(figsize=(8, 6.5))
    for column, label, color, marker in (
        ("rna_rank", "RNA", "#4c78a8", "o"),
        ("atac_rank", "ATAC", "#f2a541", "s"),
        ("multiomic_rank", "Multi-omic", "#59a14f", "D"),
    ):
        axis.scatter(ordered[column], y, label=label, color=color, marker=marker)
    for index, row in enumerate(ordered.itertuples()):
        axis.plot(
            [row.rna_rank, row.atac_rank, row.multiomic_rank],
            [index, index, index],
            color="0.85",
            linewidth=1,
            zorder=0,
        )
    axis.set_yticks(y, ordered["candidate"])
    axis.invert_xaxis()
    axis.set(
        xlabel="Rank (better →)",
        ylabel="Candidate",
        title="Modality and multi-omic rank comparison",
    )
    axis.legend(frameon=False, ncol=3)
    _finish(fig, output_path)


def plot_replicate_consistency(
    ranking: pd.DataFrame, output_path: str | Path
) -> None:
    """Plot raw replicate cosine values separately from ranking scores."""
    import matplotlib.pyplot as plt

    ordered = ranking.sort_values("multiomic_rank", ascending=False)
    y = np.arange(len(ordered))
    fig, axis = plt.subplots(figsize=(8, 6.5))
    axis.scatter(
        ordered["rna_replicate_cosine"], y, label="RNA", color="#4c78a8"
    )
    axis.scatter(
        ordered["atac_replicate_cosine"], y, label="ATAC", color="#f2a541"
    )
    for index, row in enumerate(ordered.itertuples()):
        axis.plot(
            [row.rna_replicate_cosine, row.atac_replicate_cosine],
            [index, index],
            color="0.85",
            linewidth=1,
            zorder=0,
        )
    axis.axvline(0, color="0.4", linewidth=1, linestyle="--")
    axis.set_yticks(y, ordered["candidate"])
    axis.set_xlim(-1.05, 1.05)
    axis.set(
        xlabel="R1 versus R2 cosine similarity",
        ylabel="Candidate",
        title="Replicate consistency diagnostics",
    )
    axis.legend(frameon=False)
    _finish(fig, output_path)
