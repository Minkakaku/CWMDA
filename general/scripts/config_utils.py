"""Config utilities for multi-omics workflows."""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any, Dict

import yaml


def _deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    merged = deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = deepcopy(value)
    return merged


def load_yaml(path: str | Path) -> Dict[str, Any]:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def load_config(common_path: str | Path, specific_path: str | Path) -> Dict[str, Any]:
    common = load_yaml(common_path)
    specific = load_yaml(specific_path)
    return _deep_merge(common, specific)


def resolve_path(value: str | Path, base_dir: str | Path | None = None) -> Path:
    path = Path(value)
    if not path.is_absolute() and base_dir:
        path = Path(base_dir) / path
    return path.resolve()
