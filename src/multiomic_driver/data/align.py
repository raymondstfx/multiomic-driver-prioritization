"""Guide classification and safe within-sample modality alignment."""

from __future__ import annotations

import re
from collections.abc import Iterable

import pandas as pd

_GUIDE_PATTERN = re.compile(r"^(?P<target>.+)-(?P<number>[12])$")


def normalize_barcode(barcode: str) -> str:
    """Strip surrounding whitespace while preserving the observed 10x suffix."""
    return str(barcode).strip()


def guide_target(guide_id: str) -> str:
    """Map an official `TARGET-1/2` guide call to its vector target."""
    match = _GUIDE_PATTERN.fullmatch(str(guide_id))
    if match is None:
        raise ValueError(f"Unsupported guide identifier: {guide_id}")
    return match.group("target")


def classify_guide_calls(calls: pd.DataFrame) -> pd.DataFrame:
    """Apply the authors' singlet rule to a boolean cell-by-guide matrix.

    One called guide is a singlet. Two calls are a singlet only when both
    target the same vector; all other positive combinations are multiplets.
    """
    if calls.index.duplicated().any():
        raise ValueError("Guide calls contain duplicate cell barcodes")
    targets = {column: guide_target(column) for column in calls.columns}
    rows: list[dict[str, object]] = []
    for barcode, row in calls.astype(bool).iterrows():
        detected = [column for column, value in row.items() if value]
        detected_targets = sorted({targets[guide] for guide in detected})
        if not detected:
            status = "no_call"
        elif len(detected) <= 2 and len(detected_targets) == 1:
            status = "singlet"
        else:
            status = "multiplet"
        target = detected_targets[0] if status == "singlet" else pd.NA
        is_non_targeting = (
            str(target).startswith("NTC") if status == "singlet" else pd.NA
        )
        rows.append(
            {
                "cell_barcode": normalize_barcode(barcode),
                "guide": ";".join(detected) if detected else pd.NA,
                "target_gene": target,
                "is_non_targeting": is_non_targeting,
                "guide_assignment_status": status,
                "n_guides_called": len(detected),
            }
        )
    result = pd.DataFrame(rows).set_index("cell_barcode")
    result["is_non_targeting"] = result["is_non_targeting"].astype("boolean")
    return result


def barcode_overlap_summary(
    rna_barcodes: Iterable[str],
    atac_barcodes: Iterable[str],
    guide_barcodes: Iterable[str],
) -> dict[str, int]:
    """Count all pairwise and three-way barcode intersections."""
    rna = {normalize_barcode(value) for value in rna_barcodes}
    atac = {normalize_barcode(value) for value in atac_barcodes}
    guide = {normalize_barcode(value) for value in guide_barcodes}
    return {
        "rna_cells": len(rna),
        "atac_cells": len(atac),
        "guide_cells": len(guide),
        "rna_atac": len(rna & atac),
        "rna_guide": len(rna & guide),
        "atac_guide": len(atac & guide),
        "rna_atac_guide": len(rna & atac & guide),
    }


def build_cell_metadata(
    rna_barcodes: Iterable[str],
    atac_barcodes: Iterable[str],
    guide_assignments: pd.DataFrame,
    *,
    condition: str,
    replicate: int,
    rna_sample_id: str,
    atac_sample_id: str,
    guide_sample_id: str,
) -> pd.DataFrame:
    """Build one metadata row for every RNA or ATAC cell in one experiment."""

    def checked(values: Iterable[str], modality: str) -> list[str]:
        normalized = [normalize_barcode(value) for value in values]
        if len(normalized) != len(set(normalized)):
            raise ValueError(f"Duplicate normalized {modality} barcodes")
        return normalized

    rna = checked(rna_barcodes, "RNA")
    atac = checked(atac_barcodes, "ATAC")
    if guide_assignments.index.duplicated().any():
        raise ValueError("Guide assignments contain duplicate barcodes")
    guide = guide_assignments.copy()
    guide.index = pd.Index(
        [normalize_barcode(value) for value in guide.index], name="cell_barcode"
    )
    if guide.index.duplicated().any():
        raise ValueError("Guide assignments contain duplicate normalized barcodes")

    all_barcodes = sorted(set(rna) | set(atac))
    metadata = pd.DataFrame({"cell_barcode": all_barcodes})
    metadata["cell_id"] = [
        f"{condition}:{replicate}:{barcode}" for barcode in all_barcodes
    ]
    metadata["condition"] = condition
    metadata["replicate"] = replicate
    metadata["rna_sample_id"] = rna_sample_id
    metadata["atac_sample_id"] = atac_sample_id
    metadata["guide_sample_id"] = guide_sample_id
    metadata["rna_present"] = metadata["cell_barcode"].isin(rna)
    metadata["atac_present"] = metadata["cell_barcode"].isin(atac)
    metadata = metadata.join(guide, on="cell_barcode")
    metadata["guide_capture_present"] = metadata["guide_assignment_status"].notna()
    metadata["guide_assignment_status"] = metadata[
        "guide_assignment_status"
    ].fillna("not_in_guide_calls")
    metadata["n_guides_called"] = metadata["n_guides_called"].fillna(0).astype(int)
    metadata["guide_present"] = metadata["guide_assignment_status"].eq("singlet")
    metadata["is_non_targeting"] = metadata["is_non_targeting"].astype("boolean")
    return metadata


def align_modalities(
    rna_barcodes: Iterable[str],
    atac_barcodes: Iterable[str],
    guide_assignments: pd.DataFrame,
    *,
    guide_barcode_column: str = "cell_barcode",
) -> pd.DataFrame:
    """Compatibility helper returning the three-way barcode intersection."""
    if guide_barcode_column not in guide_assignments.columns:
        raise ValueError(f"Missing guide barcode column: {guide_barcode_column}")
    rna = {normalize_barcode(value) for value in rna_barcodes}
    atac = {normalize_barcode(value) for value in atac_barcodes}
    guides = guide_assignments.copy()
    guides["cell_barcode"] = guides[guide_barcode_column].map(normalize_barcode)
    if guides["cell_barcode"].duplicated().any():
        raise ValueError("Guide assignments contain duplicate normalized barcodes")
    common = sorted(rna & atac & set(guides["cell_barcode"]))
    result = guides.set_index("cell_barcode").loc[common].reset_index()
    result["rna_present"] = True
    result["atac_present"] = True
    result["guide_present"] = True
    return result


def metadata_for_modality(
    metadata: pd.DataFrame, barcodes: Iterable[str], modality: str
) -> pd.DataFrame:
    """Order metadata exactly like a modality's barcode file for H5AD writing."""
    ordered = [normalize_barcode(value) for value in barcodes]
    indexed = metadata.set_index("cell_barcode", drop=False)
    missing = set(ordered) - set(indexed.index)
    if missing:
        raise ValueError(f"Metadata is missing {len(missing)} {modality} barcodes")
    result = indexed.loc[ordered].copy()
    result.index = pd.Index(result.pop("cell_id"), name="cell_id")
    # Nullable booleans do not have a stable HDF5 representation across AnnData versions.
    result["is_non_targeting"] = result["is_non_targeting"].map(
        {True: "true", False: "false"}, na_action="ignore"
    )
    result["is_non_targeting"] = result["is_non_targeting"].fillna("unknown")
    return result

