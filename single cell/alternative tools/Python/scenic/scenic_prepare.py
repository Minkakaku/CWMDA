from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import anndata as ad
import loompy as lp
import numpy as np

TOOLS_DIR = Path(__file__).resolve().parents[3] / "Python single cell" / "tools"
sys.path.insert(0, str(TOOLS_DIR))

from io_utils import as_path, ensure_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare loom input for pySCENIC and optionally launch the shell workflow.")
    parser.add_argument("--input-h5ad", required=True, help="Input h5ad file.")
    parser.add_argument("--output-dir", required=True, help="Output directory.")
    parser.add_argument("--loom-name", default="sce.loom", help="Output loom filename.")
    parser.add_argument("--run-shell", action="store_true", help="Run a pySCENIC shell workflow after loom export.")
    parser.add_argument("--shell-script", default=None, help="Shell script used when --run-shell is enabled.")
    parser.add_argument("--grn-output", default="grn_result.csv")
    parser.add_argument("--ctx-output", default="ctx_result.csv")
    parser.add_argument("--auc-output", default="auc_result.loom")
    return parser.parse_args()


def create_loom(input_h5ad: Path, output_dir: Path, loom_name: str) -> Path:
    adata = ad.read_h5ad(str(input_h5ad))
    matrix = adata.raw.X if adata.raw is not None else adata.X
    gene_names = np.array(adata.raw.var_names if adata.raw is not None else adata.var_names)
    cell_ids = np.array(adata.raw.obs_names if adata.raw is not None else adata.obs_names)

    loom_path = output_dir / loom_name
    row_attrs = {"Gene": gene_names}
    col_attrs = {"CellID": cell_ids}
    lp.create(str(loom_path), matrix.transpose(), row_attrs, col_attrs)
    return loom_path


def run_shell_workflow(args: argparse.Namespace, output_dir: Path, loom_path: Path) -> None:
    if not args.shell_script:
        raise ValueError("--shell-script is required when --run-shell is enabled.")

    command = [
        "bash",
        str(as_path(args.shell_script)),
        str(loom_path),
        str(output_dir / args.grn_output),
        str(output_dir / args.ctx_output),
        str(output_dir / args.auc_output),
    ]
    subprocess.run(command, check=True)


def main() -> None:
    args = parse_args()
    output_dir = ensure_dir(args.output_dir)
    loom_path = create_loom(as_path(args.input_h5ad), output_dir, args.loom_name)
    if args.run_shell:
        run_shell_workflow(args, output_dir, loom_path)
    print("SCENIC preparation finished: {}".format(output_dir))


if __name__ == "__main__":
    main()
