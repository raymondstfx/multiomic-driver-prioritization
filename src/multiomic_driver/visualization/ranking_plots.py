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


def plot_top5_candidate_evidence(
    top_candidates: pd.DataFrame, output_path: str | Path
) -> None:
    """Show modality z-scores and replicate diagnostics without combining them."""
    import matplotlib.pyplot as plt

    ordered = top_candidates.sort_values("multiomic_rank")
    x = np.arange(len(ordered))
    width = 0.34
    fig, (score_axis, reliability_axis) = plt.subplots(
        2, 1, figsize=(8.5, 7.5), sharex=True, height_ratios=(1.15, 1)
    )
    score_axis.bar(
        x - width / 2, ordered["rna_z"], width, label="RNA z-score", color="#4c78a8"
    )
    score_axis.bar(
        x + width / 2,
        ordered["atac_z"],
        width,
        label="ATAC z-score",
        color="#f2a541",
    )
    score_axis.axhline(0, color="0.4", linewidth=1)
    for position, rank in zip(x, ordered["multiomic_rank"], strict=True):
        score_axis.text(
            position,
            score_axis.get_ylim()[1] * 0.95,
            f"multi rank {int(rank)}",
            ha="center",
            va="top",
            fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.7},
        )
    score_axis.set_ylabel("Standardised evidence")
    score_axis.set_title("Top-five computational evidence (literature not scored)")
    score_axis.legend(frameon=False, ncol=2)

    reliability_axis.scatter(
        x - 0.08,
        ordered["rna_replicate_cosine"],
        label="RNA replicate cosine",
        color="#4c78a8",
    )
    reliability_axis.scatter(
        x + 0.08,
        ordered["atac_replicate_cosine"],
        label="ATAC replicate cosine",
        color="#f2a541",
    )
    reliability_axis.axhline(0, color="0.4", linewidth=1, linestyle="--")
    reliability_axis.axhline(0.5, color="0.7", linewidth=1, linestyle=":")
    reliability_axis.set(
        ylabel="R1 vs R2 cosine",
        xlabel="Candidate",
        ylim=(-1.05, 1.05),
        xticks=x,
        xticklabels=ordered["candidate"],
    )
    reliability_axis.legend(frameon=False, ncol=2)
    _finish(fig, output_path)


def plot_zfpm2_hic2_vs_ntc(
    expression_summary: pd.DataFrame,
    contrasts: pd.DataFrame,
    output_path: str | Path,
) -> None:
    """Plot descriptive ZFPM2 expression and HIC2-minus-NTC contrasts."""
    import matplotlib.pyplot as plt

    labels = ["DMSO R1", "DMSO R2", "Dasatinib R1", "Dasatinib R2"]
    keys = [("DMSO", "1"), ("DMSO", "2"), ("Dasatinib", "1"), ("Dasatinib", "2")]
    summary = expression_summary.copy()
    summary["replicate"] = summary["replicate"].astype(str)
    contrast_data = contrasts.copy()
    contrast_data["replicate"] = contrast_data["replicate"].astype(str)
    x = np.arange(len(keys))
    width = 0.35
    fig, (expression_axis, contrast_axis) = plt.subplots(
        2, 1, figsize=(9, 7.5), sharex=True, height_ratios=(1.15, 1)
    )
    for offset, group, color in (
        (-width / 2, "HIC2", "#e45756"),
        (width / 2, "NTC", "#72b7b2"),
    ):
        values = []
        for condition, replicate in keys:
            match = summary.loc[
                summary["condition"].eq(condition)
                & summary["replicate"].eq(replicate)
                & summary["group"].eq(group),
                "mean_expression",
            ]
            values.append(float(match.iloc[0]))
        expression_axis.bar(x + offset, values, width, label=group, color=color)
    expression_axis.set(
        ylabel="Mean log-normalised ZFPM2 expression",
        title="ZFPM2 expression in HIC2 and NTC singlet cells",
    )
    expression_axis.legend(frameon=False)

    contrast_values = []
    for condition, replicate in keys:
        match = contrast_data.loc[
            contrast_data["condition"].eq(condition)
            & contrast_data["replicate"].eq(replicate),
            "hic2_minus_ntc",
        ]
        contrast_values.append(float(match.iloc[0]))
    colors = np.where(np.asarray(contrast_values) >= 0, "#59a14f", "#e45756")
    contrast_axis.bar(x, contrast_values, color=colors)
    contrast_axis.axhline(0, color="0.4", linewidth=1)
    contrast_axis.set(
        ylabel="HIC2 minus NTC",
        xlabel="Condition and biological replicate",
        xticks=x,
        xticklabels=labels,
    )
    _finish(fig, output_path)
