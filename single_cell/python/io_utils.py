from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import pandas as pd


PathLike = str | os.PathLike


def as_path(path_value: PathLike) -> Path:
    return Path(path_value).expanduser().resolve()


def ensure_dir(path_value: PathLike) -> Path:
    path = as_path(path_value)
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_table(path_value: PathLike) -> pd.DataFrame:
    path = as_path(path_value)
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt", ".xls"} else ","
    return pd.read_csv(path, sep=sep)


def save_fig(
    plt_obj: object = plt,
    out_dir: PathLike = ".",
    file_name: str | None = None,
    fig_size: dict[str, float] | None = None,
) -> None:
    if not file_name:
        raise ValueError("file_name is required.")

    out_path = ensure_dir(out_dir)
    figure = plt.gcf()
    if fig_size:
        width = fig_size.get("w", fig_size.get("width", 5))
        height = fig_size.get("h", fig_size.get("height", 5))
        figure.set_size_inches(width, height)
    figure.tight_layout()
    figure.savefig(out_path / f"{file_name}.png", bbox_inches="tight", pad_inches=0.1)
    figure.savefig(out_path / f"{file_name}.pdf", bbox_inches="tight", pad_inches=0.1)
    plt.clf()


def parse_csv(raw_value: str | None) -> list[str]:
    if not raw_value:
        return []
    return [item.strip() for item in str(raw_value).split(",") if item.strip()]


def first_existing_column(columns: Iterable[str], candidates: Sequence[str]) -> str | None:
    column_set = set(columns)
    for name in candidates:
        if name in column_set:
            return name
    return None
