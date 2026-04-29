from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import ktplotspy as kpy
import pandas as pd

TOOLS_DIR = Path(__file__).resolve().parents[3] / "Python single cell" / "tools"
sys.path.insert(0, str(TOOLS_DIR))

from io_utils import as_path, ensure_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate CellPhoneDB chord diagrams from result tables.")
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--means", required=True)
    parser.add_argument("--pvalues", required=True)
    parser.add_argument("--deconvoluted", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--celltype-key", default="cell_type")
    parser.add_argument("--genes", default=None, help="Optional comma-separated genes.")
    parser.add_argument("--remove-self", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = ensure_dir(args.out_dir)
    adata = ad.read_h5ad(str(as_path(args.input_h5ad)))
    means = pd.read_csv(as_path(args.means), sep="\t")
    pvals = pd.read_csv(as_path(args.pvalues), sep="\t")
    decon = pd.read_csv(as_path(args.deconvoluted), sep="\t")
    genes = [item.strip() for item in args.genes.split(",") if item.strip()] if args.genes else None

    cell_types = sorted({str(value) for value in adata.obs[args.celltype_key].dropna().tolist()})
    for cell_type in cell_types:
        plot = kpy.plot_cpdb_chord(
            adata,
            means,
            pvals,
            decon,
            celltype_key=args.celltype_key,
            cell_type1=cell_type,
            cell_type2=".",
            remove_self=args.remove_self,
            genes=genes,
            figsize=(6, 6),
            gap=0,
            labelposition=50,
        )
        plot.save(out_dir / f"{cell_type}_chord_diagram.png", format="png", dpi=300)
        plot.save(out_dir / f"{cell_type}_chord_diagram.pdf", dpi=300)


if __name__ == "__main__":
    main()
