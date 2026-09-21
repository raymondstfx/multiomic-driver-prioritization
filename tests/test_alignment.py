import pandas as pd
import pytest

from multiomic_driver.data.align import (
    align_modalities,
    classify_guide_calls,
    guide_target,
    normalize_barcode,
)


def test_barcode_intersection() -> None:
    guides = pd.DataFrame(
        {"cell_barcode": ["AA-1", "BB-1", "DD-1"], "guide": ["g1", "g2", "g4"]}
    )
    result = align_modalities(
        ["AA-1", "BB-1", "CC-1"], ["AA-1", "BB-1", "EE-1"], guides
    )
    assert result["cell_barcode"].tolist() == ["AA-1", "BB-1"]


def test_real_barcode_suffix_is_preserved() -> None:
    assert normalize_barcode(" AAACAGCCACTAGCGT-1\n") == "AAACAGCCACTAGCGT-1"


def test_official_guide_singlet_rule() -> None:
    calls = pd.DataFrame(
        {
            "HIC2-1": [True, True, True, False],
            "HIC2-2": [False, True, False, False],
            "NTC28-1": [False, False, True, False],
        },
        index=["one", "same_target_pair", "mixed", "none"],
    )
    result = classify_guide_calls(calls)
    assert result.loc["one", "guide_assignment_status"] == "singlet"
    assert result.loc["same_target_pair", "target_gene"] == "HIC2"
    assert result.loc["mixed", "guide_assignment_status"] == "multiplet"
    assert result.loc["none", "guide_assignment_status"] == "no_call"


def test_non_targeting_guide_detection() -> None:
    calls = pd.DataFrame({"NTC28-1": [True]}, index=["cell-1"])
    result = classify_guide_calls(calls)
    assert bool(result.loc["cell-1", "is_non_targeting"])


def test_invalid_guide_column_is_rejected() -> None:
    with pytest.raises(ValueError, match="Unsupported guide identifier"):
        guide_target("HIC2_unknown")
