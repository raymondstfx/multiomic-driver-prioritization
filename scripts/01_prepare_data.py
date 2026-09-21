"""Discover GSE288996 K562 inputs and prepare aligned interim objects."""

from pathlib import Path

from _common import run_stage

from multiomic_driver.data.metadata import parse_sample_name
from multiomic_driver.utils.io import require_path


def main(config: dict) -> None:
    raw = require_path(config["paths"]["raw"], "raw GSE288996 directory")
    parseable = []
    for path in Path(raw).iterdir():
        try:
            parseable.append((path, parse_sample_name(path)))
        except ValueError:
            continue
    if not parseable:
        raise FileNotFoundError(
            f"No recognised K562 sample files found in {raw}; downloads may be incomplete."
        )
    raise NotImplementedError(
        "Files were discovered, but sample-specific matrix/feature/guide mapping must "
        "be confirmed from the completed GEO download before merging."
    )


if __name__ == "__main__":
    run_stage("01_prepare_data", main)

