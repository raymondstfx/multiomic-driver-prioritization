"""Replicate-aware background-corrected effects for Phase 4."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd

from multiomic_driver.effects.background_corrected import difference_in_differences
from multiomic_driver.utils.design import DEFAULT_DESIGN

CONDITIONS = DEFAULT_DESIGN.conditions
REPLICATES = DEFAULT_DESIGN.replicates
COUNT_COLUMNS = {
    ("DMSO", 1): "dmso_r1_cells",
    ("DMSO", 2): "dmso_r2_cells",
    ("Dasatinib", 1): "dasatinib_r1_cells",
    ("Dasatinib", 2): "dasatinib_r2_cells",
}


class MissingExperimentalGroupError(ValueError):
    """A required candidate or NTC group has no observed cells."""


class InsufficientGroupSizeError(ValueError):
    """A present group is below the configured cell-count threshold."""


def _boolean(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values.dtype):
        return values.fillna(False).astype(bool)
    return values.astype(str).str.lower().eq("true")


def validate_eligibility_for_effects(
    eligibility: pd.DataFrame, *, min_cells_per_group: int
) -> None:
    """Validate that missing and low-count exclusions remain distinguishable."""
    required = {
        "target_gene",
        *COUNT_COLUMNS.values(),
        "min_group_cells",
        "all_four_groups_present",
        "phase4_eligible",
        "exclusion_reason",
    }
    missing_columns = required - set(eligibility.columns)
    if missing_columns:
        raise ValueError(f"Eligibility table is missing: {sorted(missing_columns)}")
    for _, row in eligibility.iterrows():
        missing_groups = [
            group for group, column in COUNT_COLUMNS.items() if pd.isna(row[column])
        ]
        reason = (
            str(row["exclusion_reason"]) if pd.notna(row["exclusion_reason"]) else ""
        )
        if missing_groups:
            if pd.notna(row["min_group_cells"]):
                raise ValueError(
                    "Missing candidate groups must keep min_group_cells as NA"
                )
            if bool(row["phase4_eligible"]):
                raise ValueError("A candidate with a missing group cannot be eligible")
            for condition, replicate in missing_groups:
                token = f"missing_{condition}_R{replicate}"
                if token not in reason:
                    raise ValueError(f"Missing-group exclusion lacks reason: {token}")
        else:
            counts = [int(row[column]) for column in COUNT_COLUMNS.values()]
            low_groups = [
                (group, count)
                for group, count in zip(COUNT_COLUMNS, counts, strict=True)
                if count < min_cells_per_group
            ]
            for (condition, replicate), count in low_groups:
                token = (
                    f"below_min_{condition}_R{replicate}({count}<{min_cells_per_group})"
                )
                if token not in reason:
                    raise ValueError(f"Low-count exclusion lacks reason: {token}")
            if low_groups and bool(row["phase4_eligible"]):
                raise ValueError("A below-threshold candidate cannot be eligible")


def _group_mean(
    matrix: np.ndarray,
    mask: np.ndarray,
    *,
    label: str,
    min_cells_per_group: int,
) -> tuple[np.ndarray, int]:
    count = int(mask.sum())
    if count == 0:
        raise MissingExperimentalGroupError(f"Missing experimental group: {label}")
    if count < min_cells_per_group:
        raise InsufficientGroupSizeError(
            f"Group {label} has {count} cells; requires {min_cells_per_group}"
        )
    return np.asarray(matrix[mask].mean(axis=0)).ravel(), count


def candidate_replicate_effects(
    matrix: np.ndarray,
    metadata: pd.DataFrame,
    candidate: str,
    *,
    min_cells_per_group: int,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Calculate R1 and R2 DiD vectors without merging biological replicates."""
    if matrix.shape[0] != len(metadata):
        raise ValueError("Matrix rows must match metadata rows")
    singlet = metadata["guide_assignment_status"].astype(str).eq("singlet")
    ntc = _boolean(metadata["is_non_targeting"])
    target = metadata["target_gene"].astype(str).eq(candidate)
    conditions = metadata["condition"].astype(str)
    replicates = metadata["replicate"].astype(int)
    effects = []
    rows = []
    for replicate in REPLICATES:
        vectors: dict[tuple[str, str], np.ndarray] = {}
        counts: dict[tuple[str, str], int] = {}
        for condition in CONDITIONS:
            base = singlet & conditions.eq(condition) & replicates.eq(replicate)
            for population, selected in (
                ("candidate", base & target & ~ntc),
                ("NTC", base & ntc),
            ):
                label = f"{candidate}:{condition}:R{replicate}:{population}"
                vector, count = _group_mean(
                    matrix,
                    selected.to_numpy(),
                    label=label,
                    min_cells_per_group=min_cells_per_group,
                )
                vectors[(condition, population)] = vector
                counts[(condition, population)] = count
        effect = difference_in_differences(
            treated_perturbed=vectors[("Dasatinib", "candidate")],
            treated_control=vectors[("Dasatinib", "NTC")],
            untreated_perturbed=vectors[("DMSO", "candidate")],
            untreated_control=vectors[("DMSO", "NTC")],
        )
        effects.append(effect)
        rows.append(
            {
                "target_gene": candidate,
                "replicate": replicate,
                "dmso_candidate_cells": counts[("DMSO", "candidate")],
                "dmso_ntc_cells": counts[("DMSO", "NTC")],
                "dasatinib_candidate_cells": counts[("Dasatinib", "candidate")],
                "dasatinib_ntc_cells": counts[("Dasatinib", "NTC")],
            }
        )
    return np.vstack(effects), pd.DataFrame(rows)


def replicate_consistency(first: np.ndarray, second: np.ndarray) -> float:
    """Return cosine similarity, or NA when either effect has zero magnitude."""
    denominator = float(np.linalg.norm(first) * np.linalg.norm(second))
    if denominator == 0:
        return float("nan")
    return float(np.dot(first, second) / denominator)


def compute_modality_effects(
    matrix: np.ndarray,
    metadata: pd.DataFrame,
    candidates: Iterable[str],
    *,
    min_cells_per_group: int,
    dimension_prefix: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute replicate-specific and mean effects for one modality."""
    replicate_tables = []
    summary_rows = []
    dimensions = [f"{dimension_prefix}{index + 1}" for index in range(matrix.shape[1])]
    for candidate in candidates:
        effects, labels = candidate_replicate_effects(
            matrix,
            metadata,
            candidate,
            min_cells_per_group=min_cells_per_group,
        )
        effect_table = pd.DataFrame(effects, columns=dimensions)
        replicate_tables.append(pd.concat([labels, effect_table], axis=1))
        mean_effect = effects.mean(axis=0)
        summary = {
            "target_gene": candidate,
            "r1_effect_norm": float(np.linalg.norm(effects[0])),
            "r2_effect_norm": float(np.linalg.norm(effects[1])),
            "mean_effect_norm": float(np.linalg.norm(mean_effect)),
            "replicate_cosine_similarity": replicate_consistency(
                effects[0], effects[1]
            ),
        }
        summary.update(dict(zip(dimensions, mean_effect, strict=True)))
        summary_rows.append(summary)
    return (
        pd.concat(replicate_tables, ignore_index=True),
        pd.DataFrame(summary_rows),
    )
