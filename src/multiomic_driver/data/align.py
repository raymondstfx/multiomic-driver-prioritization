"""Barcode normalisation and cross-modality alignment."""

from __future__ import annotations

from collections.abc import Iterable

import pandas as pd


def normalize_barcode(barcode: str) -> str:
    """Normalise whitespace and a terminal 10x lane suffix only."""
    value = str(barcode).strip()
    return value.removesuffix("-1")


def align_modalities(
    rna_barcodes: Iterable[str],
    atac_barcodes: Iterable[str],
    guide_assignments: pd.DataFrame,
    *,
    guide_barcode_column: str = "cell_barcode",
) -> pd.DataFrame:
    """Return one row per barcode present in RNA, ATAC, and guide assignments."""
    if guide_barcode_column not in guide_assignments.columns:
        raise ValueError(f"Missing guide barcode column: {guide_barcode_column}")
    rna = {normalize_barcode(x) for x in rna_barcodes}
    atac = {normalize_barcode(x) for x in atac_barcodes}
    guides = guide_assignments.copy()
    guides["cell_barcode"] = guides[guide_barcode_column].map(normalize_barcode)
    if guides["cell_barcode"].duplicated().any():
        raise ValueError("Guide assignments contain duplicate normalised barcodes")
    common = sorted(rna & atac & set(guides["cell_barcode"]))
    result = guides.set_index("cell_barcode").loc[common].reset_index()
    result["rna_present"] = True
    result["atac_present"] = True
    result["guide_present"] = True
    return result

