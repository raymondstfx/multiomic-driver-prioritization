"""Reusable data-quality summaries."""

from __future__ import annotations

import pandas as pd


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

