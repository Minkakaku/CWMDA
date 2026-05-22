"""Logging helpers shared across workflows."""
from __future__ import annotations

import logging
from pathlib import Path


DEFAULT_FORMAT = "[%(asctime)s] %(levelname)s %(name)s: %(message)s"


def setup_logger(name: str, log_dir: str | Path, level: str = "INFO") -> logging.Logger:
    log_path = Path(log_dir)
    log_path.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(name)
    logger.setLevel(level.upper())

    if not logger.handlers:
        formatter = logging.Formatter(DEFAULT_FORMAT)
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        file_handler = logging.FileHandler(log_path / f"{name}.log", encoding="utf-8")
        file_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
        logger.addHandler(file_handler)

    return logger
