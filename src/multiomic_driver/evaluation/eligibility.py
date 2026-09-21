"""Phase 4 readiness summaries and perturbation eligibility rules."""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

CONDITIONS = ("DMSO", "Dasatinib")
REPLICATES = (1, 2)
GROUP_COLUMNS = {
    ("DMSO", 1): "dmso_r1_cells",
    ("DMSO", 2): "dmso_r2_cells",
    ("Dasatinib", 1): "dasatinib_r1_cells",
    ("Dasatinib", 2): "dasatinib_r2_cells",
}


def mitochondrial_qc_summary(metadata: pd.DataFrame) -> pd.DataFrame:
    """Summarize the overall mitochondrial-count distribution."""
    if "pct_counts_mt" not in metadata:
        raise ValueError("RNA metadata does not contain pct_counts_mt")
    values = pd.to_numeric(metadata["pct_counts_mt"], errors="raise")
    quantiles = values.quantile([0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1])
    names = ("min", "q25", "median", "q75", "q90", "q95", "q99", "max")
    return pd.DataFrame({"metric": names, "value": quantiles.to_numpy()})


def mitochondrial_qc_by_group(metadata: pd.DataFrame) -> pd.DataFrame:
    """Summarize mitochondrial percentages by condition and replicate."""
    required = {"pct_counts_mt", "condition", "replicate"}
    missing = required - set(metadata.columns)
    if missing:
        raise ValueError(f"RNA metadata is missing: {sorted(missing)}")
    rows = []
    for (condition, replicate), group in metadata.groupby(
        ["condition", "replicate"], observed=True, sort=True
    ):
        values = pd.to_numeric(group["pct_counts_mt"], errors="raise")
        rows.append(
            {
                "condition": str(condition),
                "replicate": int(replicate),
                "n_cells": len(values),
                "mean_pct_counts_mt": values.mean(),
                "median_pct_counts_mt": values.median(),
                "q75_pct_counts_mt": values.quantile(0.75),
                "q90_pct_counts_mt": values.quantile(0.90),
                "q95_pct_counts_mt": values.quantile(0.95),
            }
        )
    return pd.DataFrame(rows).sort_values(["condition", "replicate"])


def _normalised_boolean(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values.dtype):
        return values.fillna(False).astype(bool)
    return values.astype(str).str.lower().eq("true")


def targeting_singlets(metadata: pd.DataFrame) -> pd.DataFrame:
    """Select singlet targeting cells while excluding NTC and unknown calls."""
    required = {
        "guide_assignment_status",
        "target_gene",
        "is_non_targeting",
        "condition",
        "replicate",
    }
    missing = required - set(metadata.columns)
    if missing:
        raise ValueError(f"Perturbation metadata is missing: {sorted(missing)}")
    singlet = metadata["guide_assignment_status"].astype(str).eq("singlet")
    ntc_values = metadata["is_non_targeting"]
    if pd.api.types.is_bool_dtype(ntc_values.dtype):
        targeting = ntc_values.notna() & ~ntc_values.fillna(False).astype(bool)
    else:
        targeting = ntc_values.astype(str).str.lower().eq("false")
    return metadata.loc[
        singlet & targeting & metadata["target_gene"].notna()
    ].copy()


def count_cells_by_candidate_group(metadata: pd.DataFrame) -> pd.DataFrame:
    """Return targeting-singlet group counts, preserving absent groups as NA."""
    cells = targeting_singlets(metadata)
    if cells.empty:
        return pd.DataFrame(columns=["target_gene", *GROUP_COLUMNS.values()])
    counts = cells.groupby(
        ["target_gene", "condition", "replicate"], observed=True
    ).size()
    wide = counts.unstack(["condition", "replicate"])
    wide = wide.reindex(columns=pd.MultiIndex.from_tuples(GROUP_COLUMNS))
    wide.columns = [GROUP_COLUMNS[group] for group in wide.columns]
    wide = wide.astype("Int64").reset_index()
    return wide.sort_values("target_gene").reset_index(drop=True)


def _missing_reason(row: pd.Series) -> str:
    labels = {
        "dmso_r1_cells": "DMSO_R1",
        "dmso_r2_cells": "DMSO_R2",
        "dasatinib_r1_cells": "Dasatinib_R1",
        "dasatinib_r2_cells": "Dasatinib_R2",
    }
    missing = [labels[column] for column in labels if pd.isna(row[column])]
    return ";".join(f"missing_{label}" for label in missing)


def _low_count_reason(row: pd.Series, minimum: int) -> str:
    labels = {
        "dmso_r1_cells": "DMSO_R1",
        "dmso_r2_cells": "DMSO_R2",
        "dasatinib_r1_cells": "Dasatinib_R1",
        "dasatinib_r2_cells": "Dasatinib_R2",
    }
    low = []
    for column, label in labels.items():
        value = row[column]
        if pd.notna(value) and int(value) < minimum:
            low.append(f"below_min_{label}({int(value)}<{minimum})")
    return ";".join(low)


def evaluate_candidate_eligibility(
    metadata: pd.DataFrame,
    *,
    min_cells_per_group: int,
    require_all_groups: bool = True,
) -> pd.DataFrame:
    """Evaluate targeting perturbations without imputing absent groups as zero."""
    if min_cells_per_group < 1:
        raise ValueError("min_cells_per_group must be positive")
    result = count_cells_by_candidate_group(metadata)
    count_columns = list(GROUP_COLUMNS.values())
    if result.empty:
        return result.assign(
            min_group_cells=pd.Series(dtype="Int64"),
            all_four_groups_present=pd.Series(dtype=bool),
            meets_min_cells=pd.Series(dtype=bool),
            phase4_eligible=pd.Series(dtype=bool),
            exclusion_reason=pd.Series(dtype=str),
        )
    result["all_four_groups_present"] = result[count_columns].notna().all(axis=1)
    result["min_group_cells"] = result[count_columns].min(
        axis=1, skipna=not require_all_groups
    )
    result["min_group_cells"] = result["min_group_cells"].astype("Int64")
    result["meets_min_cells"] = result["min_group_cells"].ge(
        min_cells_per_group
    ).fillna(False)
    coverage_ok = (
        result["all_four_groups_present"] if require_all_groups else True
    )
    result["phase4_eligible"] = coverage_ok & result["meets_min_cells"]
    reasons = []
    for _, row in result.iterrows():
        missing_reason = _missing_reason(row) if require_all_groups else ""
        low_reason = _low_count_reason(row, min_cells_per_group)
        reasons.append(";".join(filter(None, [missing_reason, low_reason])))
    result["exclusion_reason"] = reasons
    return result[
        [
            "target_gene",
            *count_columns,
            "min_group_cells",
            "all_four_groups_present",
            "meets_min_cells",
            "phase4_eligible",
            "exclusion_reason",
        ]
    ]


def ntc_group_coverage(
    metadata: pd.DataFrame, *, min_cells_per_group: int
) -> pd.DataFrame:
    """Count pooled NTC singlets and assess reference coverage."""
    singlet = metadata["guide_assignment_status"].astype(str).eq("singlet")
    ntc = _normalised_boolean(metadata["is_non_targeting"])
    counts = metadata.loc[singlet & ntc].groupby(
        ["condition", "replicate"], observed=True
    ).size()
    index = pd.MultiIndex.from_product(
        [CONDITIONS, REPLICATES], names=["condition", "replicate"]
    )
    result = counts.reindex(index).rename("cell_count").reset_index()
    result["cell_count"] = result["cell_count"].astype("Int64")
    result["meets_min_cells"] = result["cell_count"].ge(
        min_cells_per_group
    ).fillna(False)
    return result


def readiness_summary(
    eligibility: pd.DataFrame,
    ntc_coverage: pd.DataFrame,
    *,
    required_candidate: str = "HIC2",
    unresolved_major_rna_qc_problem: bool = False,
) -> pd.DataFrame:
    """Build the explicit Phase 4 gate summary."""
    hic2 = eligibility.loc[eligibility["target_gene"].eq(required_candidate)]
    hic2_eligible = bool(
        len(hic2) == 1 and bool(hic2.iloc[0]["phase4_eligible"])
    )
    ntc_ready = bool(
        len(ntc_coverage) == len(GROUP_COLUMNS)
        and ntc_coverage["meets_min_cells"].all()
    )
    complete = eligibility["all_four_groups_present"]
    eligible = eligibility["phase4_eligible"]
    values: Sequence[tuple[str, object]] = (
        ("targeting_perturbations", len(eligibility)),
        ("primary_phase4_candidates", int(eligible.sum())),
        ("excluded_incomplete_coverage", int((~complete).sum())),
        ("excluded_insufficient_cells", int((complete & ~eligible).sum())),
        ("hic2_eligible", hic2_eligible),
        ("ntc_all_groups_ready", ntc_ready),
        ("unresolved_major_rna_qc_problem", unresolved_major_rna_qc_problem),
        (
            "phase4_ready",
            bool(
                ntc_ready
                and eligible.any()
                and hic2_eligible
                and not unresolved_major_rna_qc_problem
            ),
        ),
    )
    return pd.DataFrame(values, columns=["metric", "value"])
