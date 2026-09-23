"""Canonical phase-specific result paths and artifact metadata."""

from __future__ import annotations

import hashlib
from pathlib import Path

PHASE_DIRECTORIES = {
    "phase03": "phase03_qc",
    "phase035": "phase035_readiness",
    "phase04": "phase04_effects",
    "phase05": "phase05_ranking",
    "phase06": "phase06_validation",
    "phase07": "phase07_sensitivity",
}


def results_root(config: dict) -> Path:
    """Return the configured generated-results root."""
    return Path(config["paths"]["results"])


def table_dir(config: dict, phase: str) -> Path:
    """Return and create a phase-specific table directory."""
    directory = results_root(config) / "tables" / PHASE_DIRECTORIES[phase]
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def figure_dir(config: dict, phase: str) -> Path:
    """Return and create a phase-specific figure directory."""
    directory = results_root(config) / "figures" / PHASE_DIRECTORIES[phase]
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def summary_dir(config: dict) -> Path:
    """Return and create the canonical human-facing summary directory."""
    directory = results_root(config) / "summary"
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def file_sha256(path: str | Path) -> str:
    """Calculate a stable SHA-256 digest for a generated artifact."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
