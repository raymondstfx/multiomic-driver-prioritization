"""Reusable Phase 3 quality-control and metadata summaries."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import pandas as pd
from anndata import AnnData

REQUIRED_METADATA = (
    "cell_barcode",
    "condition",
    "replicate",
    "guide_assignment_status",
    "guide_present",
    "guide",
    "target_gene",
    "is_non_targeting",
)


def validate_metadata(adata: AnnData, modality: str) -> None:
    """Validate identity and guide metadata required by downstream grouping."""
    if not adata.obs_names.is_unique:
        raise ValueError(f"{modality} cell IDs are not unique")
    missing = set(REQUIRED_METADATA) - set(adata.obs.columns)
    if missing:
        raise ValueError(f"{modality} metadata is missing: {sorted(missing)}")
    if set(adata.obs["condition"].astype(str)) != {"DMSO", "Dasatinib"}:
        raise ValueError(f"{modality} contains unexpected condition values")
    if set(adata.obs["replicate"].astype(int)) != {1, 2}:
        raise ValueError(f"{modality} contains unexpected replicate values")
    observed_groups = set(
        zip(
            adata.obs["condition"].astype(str),
            adata.obs["replicate"].astype(int),
            strict=True,
        )
    )
    expected_groups = {
        ("DMSO", 1),
        ("DMSO", 2),
        ("Dasatinib", 1),
        ("Dasatinib", 2),
    }
    if observed_groups != expected_groups:
        raise ValueError(f"{modality} is missing a condition/replicate group")
    guide_present = adata.obs["guide_present"]
    if not pd.api.types.is_bool_dtype(guide_present.dtype):
        guide_present = guide_present.astype(str).str.lower().eq("true")
    singlet = adata.obs["guide_assignment_status"].astype(str).eq("singlet")
    if not guide_present.astype(bool).equals(singlet):
        raise ValueError(f"{modality} guide-present and singlet flags disagree")


def summary_table(metrics: Mapping[str, object]) -> pd.DataFrame:
    """Convert named QC metrics to a stable two-column table."""
    return pd.DataFrame(
        [{"metric": metric, "value": value} for metric, value in metrics.items()]
    )


def _is_non_targeting(metadata: pd.DataFrame) -> pd.Series:
    values = metadata["is_non_targeting"]
    if pd.api.types.is_bool_dtype(values.dtype):
        return values.fillna(False).astype(bool)
    return values.astype(str).str.lower().eq("true")


def group_counts(adata: AnnData, modality: str) -> pd.DataFrame:
    """Count retained cells by modality, condition, and replicate."""
    return (
        adata.obs.assign(modality=modality)
        .groupby(["modality", "condition", "replicate"], observed=True)
        .size()
        .rename("cell_count")
        .reset_index()
    )


def perturbation_counts(adata: AnnData, modality: str) -> pd.DataFrame:
    """Count singlet perturbations by target, condition, and replicate."""
    metadata = adata.obs.loc[
        adata.obs["guide_assignment_status"].astype(str).eq("singlet")
        & adata.obs["target_gene"].notna()
    ].copy()
    metadata["is_non_targeting"] = _is_non_targeting(metadata)
    metadata["modality"] = modality
    return (
        metadata.groupby(
            [
                "modality",
                "target_gene",
                "is_non_targeting",
                "condition",
                "replicate",
            ],
            observed=True,
        )
        .size()
        .rename("cell_count")
        .reset_index()
    )


def hic2_ntc_counts(adata: AnnData, modality: str) -> pd.DataFrame:
    """Report HIC2 and pooled NTC singlet retention in all four groups."""
    singlets = adata.obs.loc[
        adata.obs["guide_assignment_status"].astype(str).eq("singlet")
    ].copy()
    is_ntc = _is_non_targeting(singlets)
    singlets["category"] = "other"
    singlets.loc[is_ntc, "category"] = "NTC"
    singlets.loc[singlets["target_gene"].astype(str).eq("HIC2"), "category"] = (
        "HIC2"
    )
    observed = (
        singlets[singlets["category"].isin(["HIC2", "NTC"])]
        .groupby(["category", "condition", "replicate"], observed=True)
        .size()
    )
    index = pd.MultiIndex.from_product(
        [["HIC2", "NTC"], ["DMSO", "Dasatinib"], [1, 2]],
        names=["category", "condition", "replicate"],
    )
    result = observed.reindex(index, fill_value=0).rename("cell_count").reset_index()
    result.insert(0, "modality", modality)
    return result


def _replace_modality_rows(path: Path, rows: pd.DataFrame, modality: str) -> None:
    if path.exists():
        existing = pd.read_csv(path)
        existing = existing[existing["modality"] != modality]
        rows = pd.concat([existing, rows], ignore_index=True)
    rows.to_csv(path, index=False)


def update_cross_modality_tables(
    adata: AnnData, modality: str, tables_dir: str | Path
) -> pd.DataFrame:
    """Update shared Phase 3 group, perturbation, and HIC2/NTC tables."""
    validate_metadata(adata, modality)
    directory = Path(tables_dir)
    directory.mkdir(parents=True, exist_ok=True)
    _replace_modality_rows(
        directory / "post_qc_group_counts.csv", group_counts(adata, modality), modality
    )
    _replace_modality_rows(
        directory / "perturbation_cell_counts.csv",
        perturbation_counts(adata, modality),
        modality,
    )
    retention = hic2_ntc_counts(adata, modality)
    _replace_modality_rows(
        directory / "hic2_ntc_retention.csv", retention, modality
    )
    return retention


def summarize_cells(metadata: pd.DataFrame) -> pd.DataFrame:
    """Summarise cell counts by available experimental metadata columns."""
    dimensions = [
        column
        for column in ("condition", "replicate", "guide", "target_gene")
        if column in metadata.columns
    ]
    rows = [{"metric": "total_cells", "group": "all", "value": len(metadata)}]
    for column in dimensions:
        for group, count in metadata[column].value_counts(dropna=False).items():
            rows.append(
                {"metric": f"cells_by_{column}", "group": str(group), "value": count}
            )
    return pd.DataFrame(rows)
