"""Tests for canonical artifact routing and checksums."""

import hashlib

from multiomic_driver.utils.artifacts import file_sha256, summary_dir, table_dir


def test_artifact_paths_and_sha256(tmp_path) -> None:
    config = {"paths": {"results": str(tmp_path / "results")}}
    assert table_dir(config, "phase07").name == "phase07_sensitivity"
    assert summary_dir(config).name == "summary"
    artifact = summary_dir(config) / "example.txt"
    artifact.write_text("stable", encoding="utf-8")
    assert file_sha256(artifact) == hashlib.sha256(b"stable").hexdigest()
