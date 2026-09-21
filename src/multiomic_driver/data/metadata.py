"""Parse structured sample metadata from GEO filenames."""

from __future__ import annotations

import re
from pathlib import Path
from typing import TypedDict


class SampleMetadata(TypedDict):
    """Metadata encoded in a GSE288996 K562 sample name."""

    sample_id: str
    cell_line: str
    condition: str
    modality: str
    replicate: int


_SAMPLE_PATTERN = re.compile(
    r"(?P<sample_id>GSM\d+)_(?P<cell_line>[^_]+)_"
    r"(?P<condition>Dasatinib|DMSO)_(?P<modality>RNA|ATAC|guideRNA)_"
    r"(?P<replicate>\d+)(?:[_.].*)?$",
    flags=re.IGNORECASE,
)


def parse_sample_name(path: str | Path) -> SampleMetadata:
    """Parse a supported K562 sample name, including optional file suffixes."""
    name = Path(path).name
    match = _SAMPLE_PATTERN.search(name)
    if match is None:
        raise ValueError(f"Could not parse GEO sample metadata from: {name}")
    values = match.groupdict()
    canonical = {"rna": "RNA", "atac": "ATAC", "guiderna": "guideRNA"}
    return SampleMetadata(
        sample_id=values["sample_id"].upper(),
        cell_line=values["cell_line"],
        condition=(
            "Dasatinib" if values["condition"].lower() == "dasatinib" else "DMSO"
        ),
        modality=canonical[values["modality"].lower()],
        replicate=int(values["replicate"]),
    )

