from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import anndata
import ktplotspy as kpy
import matplotlib.pyplot as plt
import pandas as pd
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
from cellphonedb.utils import search_utils

TOOLS_DIR = Path(__file__).resolve().parents[3] / "Python single cell" / "tools"
sys.path.insert(0, str(TOOLS_DIR))

from io_utils import as_path, ensure_dir


warnings.filterwarnings("ignore")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generic CellPhoneDB pipeline for h5ad input.")
    parser.add_argument("--input-h5ad", required=True, help="Input h5ad file.")
    parser.add_argument("--out-dir", default="cellphonedb_result", help="Output directory.")
    parser.add_argument("--celltype-key", default="cell_type", help="Cell type column in adata.obs.")
    parser.add_argument("--cpdb-zip", required=True, help="CellPhoneDB database zip file.")
    parser.add_argument("--counts-file", default=None, help="Optional counts file override. Defaults to input-h5ad.")
    parser.add_argument("--counts-data", default="hgnc_symbol", help="Gene ID type used by counts file.")
    parser.add_argument("--iterations", type=int, default=1000)
    parser.add_argument("--threshold", type=float, default=0.01)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--pvalue", type=float, default=0.05)
    parser.add_argument("--subsampling", action="store_true", help="Enable CellPhoneDB subsampling.")
    parser.add_argument("--subsampling-log", action="store_true", help="Enable log transform during subsampling.")
    parser.add_argument("--subsampling-num-pc", type=int, default=100)
    parser.add_argument("--subsampling-num-cells", type=int, default=15000)
    parser.add_argument("--top-genes-per-celltype", type=int, default=10, help="Bubble plot genes per cell type.")
    parser.add_argument("--output-suffix", default="", help="Suffix used by CellPhoneDB outputs.")
    return parser.parse_args()


def write_metadata(adata: anndata.AnnData, celltype_key: str, out_dir: Path) -> Path:
    if celltype_key not in adata.obs.columns:
        raise KeyError("Column '{}' was not found in adata.obs.".format(celltype_key))
    meta_path = out_dir / "meta.tsv"
    adata.obs[[celltype_key]].to_csv(meta_path, index=True, sep="\t")
    return meta_path


def top_interaction_genes(
    cell_type: str,
    all_cell_types: list[str],
    significant_means: pd.DataFrame,
    deconvoluted: pd.DataFrame,
    top_n: int,
) -> list[str]:
    result = search_utils.search_analysis_results(
        query_cell_types_1=[cell_type],
        query_cell_types_2=all_cell_types,
        significant_means=significant_means,
        deconvoluted=deconvoluted,
        separator="|",
        long_format=False,
    )
    if result.empty:
        return []

    result = result.copy()
    result["Total"] = result.iloc[:, 5:].sum(axis=1)
    top_rows = result.sort_values(by="Total", ascending=False).head(top_n)
    genes = top_rows[["gene_a", "gene_b"]].values.flatten().tolist()
    return sorted({gene for gene in genes if pd.notna(gene)})


def plot_heatmap(adata: anndata.AnnData, pvalues: pd.DataFrame, celltype_key: str, out_dir: Path) -> None:
    kpy.plot_cpdb_heatmap(
        adata=adata,
        pvals=pvalues,
        celltype_key=celltype_key,
        figsize=(15, 15),
        title="Sum of significant interactions",
        symmetrical=False,
    )
    plt.savefig(out_dir / "heatmap.png", facecolor="white", bbox_inches="tight", dpi=300)
    plt.savefig(out_dir / "heatmap.pdf", facecolor="white", bbox_inches="tight")
    plt.close(plt.gcf())


def plot_bubbles(
    adata: anndata.AnnData,
    celltype_key: str,
    means: pd.DataFrame,
    pvalues: pd.DataFrame,
    deconvoluted: pd.DataFrame,
    significant_means: pd.DataFrame,
    out_dir: Path,
    top_genes_per_celltype: int,
) -> None:
    cell_types = sorted({str(value) for value in adata.obs[celltype_key].dropna().tolist()})
    for cell_type in cell_types:
        genes = top_interaction_genes(
            cell_type=cell_type,
            all_cell_types=cell_types,
            significant_means=significant_means,
            deconvoluted=deconvoluted,
            top_n=top_genes_per_celltype,
        )
        plot = kpy.plot_cpdb(
            adata=adata,
            cell_type1=cell_type,
            cell_type2=".",
            means=means,
            pvals=pvalues,
            genes=genes if genes else None,
            celltype_key=celltype_key,
            figsize=(20, 10),
            title="Interactions",
            keep_significant_only=True,
            max_size=6,
            highlight_size=0.75,
            standard_scale=True,
        )
        suffix = "Top{}_bubble".format(top_genes_per_celltype) if genes else "bubble"
        plot.save(out_dir / "{}_{}.png".format(cell_type, suffix), facecolor="white", bbox_inches="tight", dpi=300)
        plot.save(out_dir / "{}_{}.pdf".format(cell_type, suffix), facecolor="white")


def run_pipeline(args: argparse.Namespace) -> None:
    input_h5ad = as_path(args.input_h5ad)
    out_dir = ensure_dir(args.out_dir)
    counts_file = as_path(args.counts_file) if args.counts_file else input_h5ad
    cpdb_zip = as_path(args.cpdb_zip)

    adata = anndata.read_h5ad(str(input_h5ad))
    meta_path = write_metadata(adata, args.celltype_key, out_dir)

    deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
        cpdb_file_path=str(cpdb_zip),
        meta_file_path=str(meta_path),
        counts_file_path=str(counts_file),
        counts_data=args.counts_data,
        iterations=args.iterations,
        threshold=args.threshold,
        threads=args.threads,
        debug_seed=42,
        result_precision=3,
        pvalue=args.pvalue,
        subsampling=args.subsampling,
        subsampling_log=args.subsampling_log,
        subsampling_num_pc=args.subsampling_num_pc,
        subsampling_num_cells=args.subsampling_num_cells,
        separator="|",
        debug=False,
        output_path=str(out_dir),
        output_suffix=args.output_suffix,
    )

    plot_heatmap(adata, pvalues, args.celltype_key, out_dir)
    plot_bubbles(
        adata=adata,
        celltype_key=args.celltype_key,
        means=means,
        pvalues=pvalues,
        deconvoluted=deconvoluted,
        significant_means=significant_means,
        out_dir=out_dir,
        top_genes_per_celltype=args.top_genes_per_celltype,
    )

    print("CellPhoneDB finished: {}".format(out_dir))


def main() -> None:
    run_pipeline(parse_args())


if __name__ == "__main__":
    main()
