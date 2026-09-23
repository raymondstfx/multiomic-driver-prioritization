"""Configuration and filesystem utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def load_yaml(path: str | Path) -> dict[str, Any]:
    """Load a YAML mapping and raise clear errors for invalid input."""
    config_path = Path(path)
    if not config_path.is_file():
        raise FileNotFoundError(f"Configuration file does not exist: {config_path}")
    with config_path.open(encoding="utf-8") as stream:
        content = yaml.safe_load(stream) or {}
    if not isinstance(content, dict):
        raise TypeError(f"Configuration root must be a mapping: {config_path}")
    return content


def load_project_config(path: str | Path) -> dict[str, Any]:
    """Resolve and merge the data and analysis YAML files in default.yaml."""
    root_config_path = Path(path)
    root = load_yaml(root_config_path)
    base = root_config_path.resolve().parent.parent
    merged: dict[str, Any] = {}
    for key in ("data_config", "analysis_config"):
        if key not in root:
            raise ValueError(f"Missing '{key}' in {root_config_path}")
        child = Path(root[key])
        if not child.is_absolute():
            child = base / child
        merged.update(load_yaml(child))
    merged["_project_root"] = str(base)
    return merged


def require_path(path: str | Path, description: str = "input") -> Path:
    """Require an existing path with an actionable missing-data message."""
    resolved = Path(path)
    if not resolved.exists():
        raise FileNotFoundError(
            f"Missing {description}: {resolved}. Download or generate it before "
            "running this pipeline stage."
        )
    return resolved
