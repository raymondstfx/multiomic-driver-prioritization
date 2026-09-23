"""Consistent console and file logging."""

from __future__ import annotations

import logging
from pathlib import Path


def configure_logging(
    name: str, log_dir: str | Path = "results/logs"
) -> logging.Logger:
    """Create an idempotent logger that writes to console and a project log."""
    logger = logging.getLogger(name)
    if logger.handlers:
        return logger
    logger.setLevel(logging.INFO)
    directory = Path(log_dir)
    directory.mkdir(parents=True, exist_ok=True)
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    handlers = [logging.StreamHandler(), logging.FileHandler(directory / f"{name}.log")]
    for handler in handlers:
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger
