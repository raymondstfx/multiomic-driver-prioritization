"""Descriptive Phase 6 evaluation of fixed Phase 5 rankings."""

from __future__ import annotations

from collections.abc import Iterable, Mapping

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REQUIRED_COLUMNS = (
    "candidate",
    "rna_score",
    "atac_score",
    "rna_z",
    "atac_z",
    "multiomic_score",
    "rna_rank",
    "atac_rank",
    "multiomic_rank",
    "rna_replicate_cosine",
    "atac_replicate_cosine",
)


def validate_ranking(
    ranking: pd.DataFrame,
    *,
    expected_count: int = 12,
    required_candidates: Iterable[str] = ("HIC2",),
    excluded_candidates: Iterable[str] = ("ZFPM2",),
) -> None:
    """Fail if a fixed ranking is malformed or contains external validation genes."""
    missing_columns = sorted(set(REQUIRED_COLUMNS) - set(ranking.columns))
    if missing_columns:
        raise ValueError(f"Ranking is missing required columns: {missing_columns}")
    if len(ranking) != expected_count:
        raise ValueError(
            f"Expected {expected_count} primary candidates, found {len(ranking)}"
        )
    if ranking["candidate"].isna().any() or not ranking["candidate"].is_unique:
        raise ValueError("Candidate names must be present and unique")

    candidates = set(ranking["candidate"].astype(str))
    absent = sorted(set(required_candidates) - candidates)
    forbidden = sorted(set(excluded_candidates) & candidates)
    if absent:
        raise ValueError(f"Required candidates are absent: {absent}")
    if forbidden:
        raise ValueError(f"External validation genes must not be ranked: {forbidden}")

    numeric = [column for column in REQUIRED_COLUMNS if column != "candidate"]
    values = ranking[numeric].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Ranking scores, ranks, and replicate diagnostics must be finite")
    expected_ranks = list(range(1, len(ranking) + 1))
    for column in ("rna_rank", "atac_rank", "multiomic_rank"):
        ranks = ranking[column].to_numpy(dtype=float)
        if not np.equal(ranks, np.floor(ranks)).all():
            raise ValueError(f"{column} must contain integer ranks")
        if sorted(ranks.astype(int).tolist()) != expected_ranks:
            raise ValueError(f"{column} must be a permutation of 1..{len(ranking)}")


def assign_cross_modal_profile(rna_z: float, atac_z: float) -> str:
    """Assign a descriptive profile from the signs of fixed modality z-scores."""
    if not np.isfinite([rna_z, atac_z]).all():
        raise ValueError("Cross-modal profiles require finite z-scores")
    if rna_z > 0 and atac_z > 0:
        return "both_high"
    if rna_z > 0 and atac_z <= 0:
        return "rna_dominant"
    if rna_z <= 0 and atac_z > 0:
        return "atac_dominant"
    return "both_below_average"


def rank_comparison_table(ranking: pd.DataFrame) -> pd.DataFrame:
    """Calculate descriptive shifts without changing any Phase 5 value."""
    result = ranking[
        [
            "candidate",
            "rna_rank",
            "atac_rank",
            "multiomic_rank",
            "rna_z",
            "atac_z",
            "multiomic_score",
        ]
    ].copy()
    result["rna_to_multiomic_shift"] = result["rna_rank"] - result["multiomic_rank"]
    result["atac_to_multiomic_shift"] = (
        result["atac_rank"] - result["multiomic_rank"]
    )
    result["cross_modal_profile"] = [
        assign_cross_modal_profile(row.rna_z, row.atac_z)
        for row in result.itertuples()
    ]
    return result.sort_values("multiomic_rank").reset_index(drop=True)


def rank_correlation_summary(ranking: pd.DataFrame) -> pd.DataFrame:
    """Return descriptive Spearman correlations for ranks and modality z-scores."""
    comparisons = (
        ("RNA rank vs ATAC rank", "rna_rank", "atac_rank"),
        ("RNA rank vs multi-omic rank", "rna_rank", "multiomic_rank"),
        ("ATAC rank vs multi-omic rank", "atac_rank", "multiomic_rank"),
        ("RNA z-score vs ATAC z-score", "rna_z", "atac_z"),
    )
    rows = []
    for label, left, right in comparisons:
        rho = float(spearmanr(ranking[left], ranking[right]).statistic)
        rows.append(
            {
                "comparison": label,
                "left_variable": left,
                "right_variable": right,
                "spearman_rho": rho,
                "n_candidates": len(ranking),
                "interpretation_scope": "descriptive; small candidate set",
            }
        )
    return pd.DataFrame(rows)


def consistency_category(
    value: float, *, moderate_threshold: float = 0.5, positive_threshold: float = 0.0
) -> str:
    """Categorise a cosine diagnostic without changing the ranking score."""
    if not np.isfinite(value):
        return "not_estimable"
    if value >= moderate_threshold:
        return "moderate/high agreement"
    if value >= positive_threshold:
        return "weak positive agreement"
    return "directional inconsistency"


def replicate_reliability_table(
    ranking: pd.DataFrame,
    *,
    moderate_threshold: float = 0.5,
    positive_threshold: float = 0.0,
) -> pd.DataFrame:
    """Annotate replicate cosines separately from score calculation."""
    result = ranking[
        ["candidate", "rna_replicate_cosine", "atac_replicate_cosine"]
    ].copy()
    result["rna_consistency_category"] = result["rna_replicate_cosine"].map(
        lambda value: consistency_category(
            value,
            moderate_threshold=moderate_threshold,
            positive_threshold=positive_threshold,
        )
    )
    result["atac_consistency_category"] = result["atac_replicate_cosine"].map(
        lambda value: consistency_category(
            value,
            moderate_threshold=moderate_threshold,
            positive_threshold=positive_threshold,
        )
    )
    return result.merge(
        ranking[["candidate", "multiomic_rank"]], on="candidate"
    ).sort_values("multiomic_rank")


def hic2_validation_table(ranking: pd.DataFrame) -> pd.DataFrame:
    """Extract HIC2 only after the label-free ranking has been fixed."""
    match = ranking.loc[ranking["candidate"].eq("HIC2")].copy()
    if len(match) != 1:
        raise ValueError(f"Expected exactly one HIC2 row, found {len(match)}")
    columns = [
        "candidate",
        "rna_rank",
        "atac_rank",
        "multiomic_rank",
        "rna_score",
        "atac_score",
        "multiomic_score",
        "rna_replicate_cosine",
        "atac_replicate_cosine",
    ]
    result = match[columns].copy()
    result["validation_role"] = "experimentally supported pooled perturbation"
    return result


def top_candidate_comparison(ranking: pd.DataFrame, top_n: int = 5) -> pd.DataFrame:
    """Return fixed evidence and consistency diagnostics for the top candidates."""
    columns = [
        "candidate",
        "rna_score",
        "atac_score",
        "rna_z",
        "atac_z",
        "multiomic_score",
        "multiomic_rank",
        "rna_replicate_cosine",
        "atac_replicate_cosine",
    ]
    result = ranking.nsmallest(top_n, "multiomic_rank")[columns].copy()
    result["cross_modal_profile"] = [
        assign_cross_modal_profile(row.rna_z, row.atac_z)
        for row in result.itertuples()
    ]
    return result


def literature_evidence_table(
    top_candidates: pd.DataFrame, evidence: Mapping[str, Mapping[str, object]]
) -> pd.DataFrame:
    """Attach qualitative, source-linked literature notes without scoring them."""
    rows = []
    for row in top_candidates.sort_values("multiomic_rank").itertuples():
        item = evidence.get(row.candidate, {})
        rows.append(
            {
                "candidate": row.candidate,
                "multiomic_rank": row.multiomic_rank,
                "rna_rank": row.rna_rank,
                "atac_rank": row.atac_rank,
                "replicate_summary": (
                    f"RNA cosine={row.rna_replicate_cosine:.3f}; "
                    f"ATAC cosine={row.atac_replicate_cosine:.3f}"
                ),
                "evidence_category": item.get(
                    "evidence_category", "not manually reviewed"
                ),
                "literature_summary": item.get("literature_summary", ""),
                "direct_relevance_to_dasatinib": item.get(
                    "direct_relevance_to_dasatinib", "unknown"
                ),
                "direct_relevance_to_k562_or_cml": item.get(
                    "direct_relevance_to_k562_or_cml", "unknown"
                ),
                "source": item.get("source", ""),
                "notes": item.get("notes", ""),
            }
        )
    return pd.DataFrame(rows)


def final_candidate_evaluation(
    ranking: pd.DataFrame,
    rank_comparison: pd.DataFrame,
    literature: pd.DataFrame,
) -> pd.DataFrame:
    """Combine fixed computational results with explicitly qualitative annotations."""
    columns = [
        "candidate",
        "multiomic_rank",
        "rna_rank",
        "atac_rank",
        "rna_z",
        "atac_z",
        "multiomic_score",
        "rna_replicate_cosine",
        "atac_replicate_cosine",
    ]
    result = ranking[columns].merge(
        rank_comparison[["candidate", "cross_modal_profile"]], on="candidate"
    )
    result["experimental_validation_status"] = np.where(
        result["candidate"].eq("HIC2"),
        "experimentally supported pooled perturbation",
        "not experimentally validated in this analysis",
    )
    result = result.merge(
        literature[["candidate", "literature_summary"]],
        on="candidate",
        how="left",
    )
    result["literature_summary"] = result["literature_summary"].fillna(
        "Not included in the focused top-five literature review."
    )
    result["interpretation"] = (
        "Computational prioritisation: "
        + result["cross_modal_profile"].str.replace("_", " ")
        + "; rank is not causal proof."
    )
    return result.sort_values("multiomic_rank").reset_index(drop=True)
