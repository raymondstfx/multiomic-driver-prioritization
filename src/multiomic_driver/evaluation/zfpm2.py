"""Targeted downstream ZFPM2 validation outside perturbation ranking."""

from __future__ import annotations

import numpy as np
import pandas as pd

from multiomic_driver.utils.design import DEFAULT_DESIGN

REQUIRED_GROUPS = ("HIC2", "NTC")
REQUIRED_CONDITIONS = DEFAULT_DESIGN.conditions


def summarize_gene_expression(
    expression: np.ndarray,
    metadata: pd.DataFrame,
    *,
    target_candidate: str = "HIC2",
    raw_counts: np.ndarray | None = None,
    library_sizes: np.ndarray | None = None,
) -> pd.DataFrame:
    """Summarise complementary expression metrics in HIC2 and NTC singlets."""
    values = np.asarray(expression, dtype=float).reshape(-1)
    if len(values) != len(metadata):
        raise ValueError("Expression and metadata must have the same number of cells")
    if not np.isfinite(values).all():
        raise ValueError("Expression values must be finite")
    required = {
        "condition",
        "replicate",
        "guide_assignment_status",
        "target_gene",
        "is_non_targeting",
    }
    missing = sorted(required - set(metadata.columns))
    if missing:
        raise ValueError(f"Metadata is missing required columns: {missing}")

    frame = metadata[list(required)].copy()
    frame["expression"] = values
    if raw_counts is not None:
        counts = np.asarray(raw_counts, dtype=float).reshape(-1)
        if len(counts) != len(frame) or (counts < 0).any():
            raise ValueError("Raw counts must be non-negative and match metadata")
        frame["raw_count"] = counts
    if library_sizes is not None:
        libraries = np.asarray(library_sizes, dtype=float).reshape(-1)
        if len(libraries) != len(frame) or (libraries < 0).any():
            raise ValueError("Library sizes must be non-negative and match metadata")
        frame["library_size"] = libraries
    singlet = frame["guide_assignment_status"].astype(str).eq("singlet")
    is_target = frame["target_gene"].astype(str).eq(target_candidate)
    is_ntc = frame["is_non_targeting"].astype(str).str.lower().eq("true")
    frame = frame.loc[singlet & (is_target | is_ntc)].copy()
    frame["group"] = np.where(is_target.loc[frame.index], target_candidate, "NTC")

    aggregation: dict[str, tuple[str, object]] = {
        "n_cells": ("expression", "size"),
        "mean_expression": ("expression", "mean"),
        "median_expression": ("expression", "median"),
        "fraction_expressing": ("expression", lambda value: (value > 0).mean()),
        "mean_positive_expression": (
            "expression",
            lambda value: value[value > 0].mean() if (value > 0).any() else 0.0,
        ),
    }
    if raw_counts is not None:
        aggregation["total_raw_count"] = ("raw_count", "sum")
        aggregation["mean_raw_count"] = ("raw_count", "mean")
    if library_sizes is not None:
        aggregation["total_library_size"] = ("library_size", "sum")
    summary = (
        frame.groupby(["condition", "replicate", "group"], observed=True)
        .agg(**aggregation)
        .reset_index()
    )
    if {"total_raw_count", "total_library_size"}.issubset(summary):
        summary["pseudobulk_cpm"] = np.divide(
            summary["total_raw_count"] * 1_000_000.0,
            summary["total_library_size"],
            out=np.zeros(len(summary), dtype=float),
            where=summary["total_library_size"].to_numpy() > 0,
        )
        summary["pseudobulk_log1p_cpm"] = np.log1p(summary["pseudobulk_cpm"])
    return summary.sort_values(["condition", "replicate", "group"]).reset_index(
        drop=True
    )


def hic2_ntc_contrasts(
    summary: pd.DataFrame,
    *,
    metric: str = "mean_expression",
    target_candidate: str = "HIC2",
) -> pd.DataFrame:
    """Calculate HIC2-minus-NTC contrasts for a selected expression metric."""
    if metric not in summary:
        raise ValueError(f"Expression summary does not contain metric {metric}")
    pivot = summary.pivot(
        index=["condition", "replicate"], columns="group", values=metric
    )
    required_groups = (target_candidate, "NTC")
    missing = sorted(set(required_groups) - set(pivot.columns))
    if missing or pivot[list(required_groups)].isna().any().any():
        raise ValueError(
            "Every condition/replicate requires target and NTC; "
            f"target={target_candidate}, missing={missing}"
        )
    contrasts = pivot.reset_index()
    contrasts["hic2_minus_ntc"] = (
        contrasts[target_candidate] - contrasts["NTC"]
    )
    contrasts["metric"] = metric
    contrasts["target_candidate"] = target_candidate
    return contrasts.rename(
        columns={
            target_candidate: "hic2_mean_expression",
            "NTC": "ntc_mean_expression",
        }
    )[
        [
            "condition",
            "replicate",
            "metric",
            "target_candidate",
            "hic2_mean_expression",
            "ntc_mean_expression",
            "hic2_minus_ntc",
        ]
    ].sort_values(["replicate", "condition"])


def replicate_aware_did(contrasts: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute treatment-specific ZFPM2 DiD separately before summarising replicates."""
    pivot = contrasts.pivot(
        index="replicate", columns="condition", values="hic2_minus_ntc"
    )
    missing = sorted(set(REQUIRED_CONDITIONS) - set(pivot.columns))
    if missing or pivot[list(REQUIRED_CONDITIONS)].isna().any().any():
        raise ValueError(
            f"Each replicate requires DMSO and Dasatinib; missing {missing}"
        )
    effects = pivot.reset_index()
    effects["replicate_specific_did"] = effects["Dasatinib"] - effects["DMSO"]
    effects = effects.rename(
        columns={
            "DMSO": "dmso_hic2_minus_ntc",
            "Dasatinib": "dasatinib_hic2_minus_ntc",
        }
    )[
        [
            "replicate",
            "dmso_hic2_minus_ntc",
            "dasatinib_hic2_minus_ntc",
            "replicate_specific_did",
        ]
    ].sort_values("replicate")
    metric = (
        str(contrasts["metric"].iloc[0])
        if "metric" in contrasts and len(contrasts)
        else "mean_expression"
    )
    target_candidate = (
        str(contrasts["target_candidate"].iloc[0])
        if "target_candidate" in contrasts and len(contrasts)
        else "HIC2"
    )
    effects.insert(0, "metric", metric)
    did = effects["replicate_specific_did"].to_numpy(dtype=float)
    signs = np.sign(did)
    agreement = bool(len(signs) > 1 and np.all(signs == signs[0]) and signs[0] != 0)
    summary = pd.DataFrame(
        [
            {
                "gene": "ZFPM2",
                "target_perturbation": target_candidate,
                "metric": metric,
                "n_replicates": len(effects),
                "mean_replicate_did": did.mean(),
                "replicate_direction_agreement": agreement,
                "interpretation_scope": (
                    "descriptive downstream validation; not a perturbation rank "
                    "or proof of direct regulation"
                ),
            }
        ]
    )
    return effects.reset_index(drop=True), summary
