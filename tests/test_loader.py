from pathlib import Path

import pytest

from multiomic_driver.data.loader import load_guide_assignments
from multiomic_driver.data.metadata import parse_sample_name


def test_metadata_parsing() -> None:
    result = parse_sample_name("GSM8780564_K562_Dasatinib_RNA_1")
    assert result == {
        "sample_id": "GSM8780564",
        "cell_line": "K562",
        "condition": "Dasatinib",
        "modality": "RNA",
        "replicate": 1,
    }


def test_missing_file_error() -> None:
    with pytest.raises(FileNotFoundError, match="Required input file"):
        load_guide_assignments(Path("definitely_missing.tsv"))

