import pandas as pd

from multiomic_driver.data.align import align_modalities


def test_barcode_intersection() -> None:
    guides = pd.DataFrame(
        {"cell_barcode": ["AA-1", "BB-1", "DD-1"], "guide": ["g1", "g2", "g4"]}
    )
    result = align_modalities(["AA-1", "BB-1", "CC-1"], ["AA", "BB", "EE"], guides)
    assert result["cell_barcode"].tolist() == ["AA", "BB"]

