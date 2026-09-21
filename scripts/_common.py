"""Shared command-line plumbing for pipeline entry points."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Callable
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from multiomic_driver.utils.io import load_project_config
from multiomic_driver.utils.logging import configure_logging


def run_stage(name: str, action: Callable[[dict], None]) -> None:
    """Load configuration, execute a stage, and report actionable failures."""
    parser = argparse.ArgumentParser(description=name.replace("_", " ").title())
    parser.add_argument("--config", default="configs/default.yaml")
    arguments = parser.parse_args()
    logger = configure_logging(name)
    try:
        logger.info("Starting stage with config %s", arguments.config)
        action(load_project_config(arguments.config))
        logger.info("Stage completed")
    except (FileNotFoundError, ValueError, NotImplementedError) as error:
        logger.error("%s", error)
        raise SystemExit(2) from error

