"""Post-hoc evaluation against externally supported candidates."""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

KNOWN_VALIDATED = ("HIC2", "ZFPM2")


def get_candidate_rank(
    ranking_df: pd.DataFrame, gene: str, rank_column: str = "multiomic_rank"
) -> int | None:
    """Return a candidate rank, or None when the candidate is absent."""
    match = ranking_df.loc[ranking_df["candidate"] == gene, rank_column]
    return None if match.empty else int(match.iloc[0])


def evaluate_known_candidates(
    ranking_df: pd.DataFrame, known_genes: Sequence[str] = KNOWN_VALIDATED
) -> pd.DataFrame:
    """Extract available modality and multi-omic ranks for known genes."""
    rank_columns = [
        name
        for name in ("rna_rank", "atac_rank", "multiomic_rank")
        if name in ranking_df.columns
    ]
    indexed = ranking_df.set_index("candidate")
    rows = []
    for gene in known_genes:
        row = {"gene": gene, "present": gene in indexed.index}
        for column in rank_columns:
            row[column] = indexed.at[gene, column] if gene in indexed.index else pd.NA
        rows.append(row)
    return pd.DataFrame(rows)

