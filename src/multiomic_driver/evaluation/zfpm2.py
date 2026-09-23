"""Targeted downstream ZFPM2 validation outside perturbation ranking."""

from __future__ import annotations

import numpy as np
import pandas as pd

REQUIRED_GROUPS = ("HIC2", "NTC")
REQUIRED_CONDITIONS = ("DMSO", "Dasatinib")


def summarize_gene_expression(
    expression: np.ndarray,
    metadata: pd.DataFrame,
    *,
    target_candidate: str = "HIC2",
) -> pd.DataFrame:
    """Summarise log-normalised expression in singlet HIC2 and NTC cells."""
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
    singlet = frame["guide_assignment_status"].astype(str).eq("singlet")
    is_target = frame["target_gene"].astype(str).eq(target_candidate)
    is_ntc = frame["is_non_targeting"].astype(str).str.lower().eq("true")
    frame = frame.loc[singlet & (is_target | is_ntc)].copy()
    frame["group"] = np.where(is_target.loc[frame.index], target_candidate, "NTC")

    summary = (
        frame.groupby(["condition", "replicate", "group"], observed=True)
        .agg(
            n_cells=("expression", "size"),
            mean_expression=("expression", "mean"),
            median_expression=("expression", "median"),
            fraction_expressing=("expression", lambda value: (value > 0).mean()),
        )
        .reset_index()
    )
    return summary.sort_values(["condition", "replicate", "group"]).reset_index(
        drop=True
    )


def hic2_ntc_contrasts(summary: pd.DataFrame) -> pd.DataFrame:
    """Calculate HIC2-minus-NTC mean-expression contrasts by condition/replicate."""
    pivot = summary.pivot(
        index=["condition", "replicate"], columns="group", values="mean_expression"
    )
    missing = sorted(set(REQUIRED_GROUPS) - set(pivot.columns))
    if missing or pivot[list(REQUIRED_GROUPS)].isna().any().any():
        raise ValueError(f"Every condition/replicate requires HIC2 and NTC; missing {missing}")
    contrasts = pivot.reset_index()
    contrasts["hic2_minus_ntc"] = contrasts["HIC2"] - contrasts["NTC"]
    return contrasts.rename(
        columns={"HIC2": "hic2_mean_expression", "NTC": "ntc_mean_expression"}
    )[
        [
            "condition",
            "replicate",
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
        raise ValueError(f"Each replicate requires DMSO and Dasatinib; missing {missing}")
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
    did = effects["replicate_specific_did"].to_numpy(dtype=float)
    signs = np.sign(did)
    agreement = bool(len(signs) > 1 and np.all(signs == signs[0]) and signs[0] != 0)
    summary = pd.DataFrame(
        [
            {
                "gene": "ZFPM2",
                "target_perturbation": "HIC2",
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
