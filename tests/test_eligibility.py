"""Synthetic tests for the Phase 4 readiness gate."""

from __future__ import annotations

import pandas as pd

from multiomic_driver.evaluation.eligibility import (
    evaluate_candidate_eligibility,
    mitochondrial_qc_by_group,
    mitochondrial_qc_summary,
    ntc_group_coverage,
    readiness_summary,
)


def _cells(
    target: str,
    condition: str,
    replicate: int,
    count: int,
    *,
    ntc: bool = False,
) -> list[dict[str, object]]:
    return [
        {
            "target_gene": target,
            "condition": condition,
            "replicate": replicate,
            "guide_assignment_status": "singlet",
            "is_non_targeting": ntc,
        }
        for _ in range(count)
    ]


def _coverage_metadata() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    groups = [("DMSO", 1), ("DMSO", 2), ("Dasatinib", 1), ("Dasatinib", 2)]
    for condition, replicate in groups:
        rows.extend(_cells("HIC2", condition, replicate, 25))
        rows.extend(_cells("LOW", condition, replicate, 10 if replicate == 2 else 25))
        rows.extend(_cells("NTC1", condition, replicate, 40, ntc=True))
    for condition, replicate in groups[:-1]:
        rows.extend(_cells("MISSING", condition, replicate, 25))
    rows.append(
        {
            "target_gene": pd.NA,
            "condition": "DMSO",
            "replicate": 1,
            "guide_assignment_status": "multiplet",
            "is_non_targeting": pd.NA,
        }
    )
    return pd.DataFrame(rows)


def test_candidate_eligibility_preserves_missing_groups() -> None:
    result = evaluate_candidate_eligibility(
        _coverage_metadata(), min_cells_per_group=20
    ).set_index("target_gene")

    assert bool(result.at["HIC2", "phase4_eligible"])
    assert not bool(result.at["LOW", "phase4_eligible"])
    assert "below_min" in result.at["LOW", "exclusion_reason"]
    assert not bool(result.at["MISSING", "phase4_eligible"])
    assert pd.isna(result.at["MISSING", "dasatinib_r2_cells"])
    assert pd.isna(result.at["MISSING", "min_group_cells"])
    assert result.at["MISSING", "exclusion_reason"] == "missing_Dasatinib_R2"
    assert "NTC1" not in result.index


def test_ntc_and_readiness_gate() -> None:
    metadata = _coverage_metadata()
    eligibility = evaluate_candidate_eligibility(
        metadata, min_cells_per_group=20
    )
    ntc = ntc_group_coverage(metadata, min_cells_per_group=20)
    summary = readiness_summary(eligibility, ntc).set_index("metric")["value"]

    assert ntc["meets_min_cells"].all()
    assert summary["targeting_perturbations"] == 3
    assert summary["primary_phase4_candidates"] == 1
    assert summary["excluded_incomplete_coverage"] == 1
    assert summary["excluded_insufficient_cells"] == 1
    assert bool(summary["hic2_eligible"])
    assert bool(summary["phase4_ready"])


def test_mitochondrial_summaries() -> None:
    metadata = pd.DataFrame(
        {
            "condition": ["DMSO", "DMSO", "Dasatinib", "Dasatinib"],
            "replicate": [1, 2, 1, 2],
            "pct_counts_mt": [10.0, 20.0, 30.0, 40.0],
        }
    )
    overall = mitochondrial_qc_summary(metadata).set_index("metric")["value"]
    by_group = mitochondrial_qc_by_group(metadata)

    assert overall["min"] == 10.0
    assert overall["median"] == 25.0
    assert overall["max"] == 40.0
    assert len(by_group) == 4
    assert set(by_group["n_cells"]) == {1}
