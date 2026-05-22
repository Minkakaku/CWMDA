"""File helpers for reusable validations and output layout."""
from __future__ import annotations

from pathlib import Path


def ensure_dir(path: str | Path) -> Path:
    target = Path(path)
    target.mkdir(parents=True, exist_ok=True)
    return target


def validate_file(path: str | Path, label: str) -> Path:
    target = Path(path)
    if not target.exists():
        raise FileNotFoundError(f"Missing {label}: {target}")
    return target
