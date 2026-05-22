"""Shared QC report utilities."""
from __future__ import annotations

from pathlib import Path
from typing import Dict

from .file_utils import ensure_dir


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="zh">
<head>
  <meta charset="UTF-8" />
  <title>QC Report</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 24px; }
    table { border-collapse: collapse; width: 60%; }
    th, td { border: 1px solid #ddd; padding: 8px; }
    th { background-color: #f4f4f4; }
  </style>
</head>
<body>
  <h1>质控报告</h1>
  <table>
    <tr><th>指标</th><th>数值</th></tr>
    {rows}
  </table>
</body>
</html>
"""


def write_qc_report(stats: Dict[str, str], output_dir: str | Path, tsv_name: str, html_name: str) -> None:
    output_path = ensure_dir(output_dir)
    tsv_path = output_path / tsv_name
    html_path = output_path / html_name

    with tsv_path.open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        for key, value in stats.items():
            handle.write(f"{key}\t{value}\n")

    rows = "\n".join(
        f"<tr><td>{key}</td><td>{value}</td></tr>" for key, value in stats.items()
    )
    html_path.write_text(HTML_TEMPLATE.format(rows=rows), encoding="utf-8")
