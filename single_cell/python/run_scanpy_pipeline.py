from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Optional

import pandas as pd
import scanpy as sc

from scanpy_pipeline import (
    cell_filtering,
    celltype_annotation,
    configure_scanpy,
    data_qc,
    dim_reduction_cluster,
    find_marker_gene,
    get_10x_mtx,
    make_report_dir,
)


STEP_OUTPUTS = {
    "load": Path("01-CellRanger") / "adata_raw.h5ad",
    "qc": Path("02-Data_QC") / "adata_qc.h5ad",
    "filter": Path("03-Cell_filter") / "adata_counts.h5ad",
    "cluster": Path("04-Dimension_reduction_cluster") / "adata_cluster.h5ad",
    "marker": Path("05-Marker_gene_analysis") / "adata_marker.h5ad",
    "annotate": Path("06-Cell_type_annotation") / "result.h5ad",
}


def _as_path(path_value: str) -> Path:
    return Path(path_value).expanduser().resolve()


def _resolve_out_root(args: argparse.Namespace) -> Path:
    if getattr(args, "out_root", None):
        return _as_path(args.out_root)
    if getattr(args, "sample_sheet", None):
        return _as_path(args.sample_sheet).parent
    if getattr(args, "input_h5ad", None):
        return _as_path(args.input_h5ad).parents[1]
    raise ValueError("Please provide --out-root, or pass --sample-sheet / --input-h5ad located under the output directory.")


def _default_input_h5ad(out_root: Path, step_name: str) -> Path:
    step_map = {
        "qc": STEP_OUTPUTS["load"],
        "filter": STEP_OUTPUTS["qc"],
        "cluster": STEP_OUTPUTS["filter"],
        "marker": STEP_OUTPUTS["cluster"],
        "annotate": STEP_OUTPUTS["marker"],
    }
    return out_root / step_map[step_name]


def _read_h5ad(path_value: Path) -> sc.AnnData:
    if not path_value.exists():
        raise FileNotFoundError("Input h5ad not found: {}".format(path_value))
    return sc.read_h5ad(path_value)


def _read_celltype_map(path_value: Optional[str]) -> Optional[Dict[str, str]]:
    if not path_value:
        return None

    path = _as_path(path_value)
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt", ".xls"} else ","
    df = pd.read_csv(path, sep=sep)
    required_columns = {"cluster", "cell_type"}
    missing_columns = required_columns.difference(df.columns)
    if missing_columns:
        raise ValueError(
            "celltype map is missing required columns: {}".format(", ".join(sorted(missing_columns)))
        )
    return {
        str(cluster): str(cell_type)
        for cluster, cell_type in zip(df["cluster"], df["cell_type"])
    }


def _configure_runtime(out_root: Path, cache_dir: Optional[str]) -> None:
    runtime_cache = _as_path(cache_dir) if cache_dir else out_root / ".scanpy_cache"
    configure_scanpy(runtime_cache)


def run_load(args: argparse.Namespace) -> Path:
    out_root = _resolve_out_root(args)
    _configure_runtime(out_root, args.cache_dir)
    make_report_dir(out_root)

    adata = get_10x_mtx(sample_sheet=args.sample_sheet, out_dir=out_root)
    output_path = out_root / STEP_OUTPUTS["load"]
    adata.write(output_path)
    print("load finished: {}".format(output_path))
    return output_path


def run_qc(args: argparse.Namespace) -> Path:
    out_root = _resolve_out_root(args)
    _configure_runtime(out_root, args.cache_dir)

    input_h5ad = _as_path(args.input_h5ad) if args.input_h5ad else _default_input_h5ad(out_root, "qc")
    adata = _read_h5ad(input_h5ad)
    adata = data_qc(adata, out_root / "02-Data_QC")
    output_path = out_root / STEP_OUTPUTS["qc"]
    adata.write(output_path)
    print("qc finished: {}".format(output_path))
    return output_path


def run_filter(args: argparse.Namespace) -> Path:
    out_root = _resolve_out_root(args)
    _configure_runtime(out_root, args.cache_dir)

    input_h5ad = _as_path(args.input_h5ad) if args.input_h5ad else _default_input_h5ad(out_root, "filter")
    adata = _read_h5ad(input_h5ad)
    adata = cell_filtering(
        _adata=adata,
        min_genes=args.min_genes,
        min_cells=args.min_cells,
        pct_counts_mt=args.pct_counts_mt,
        pct_counts_ribo=args.pct_counts_ribo,
        upper_lim=args.upper_lim,
        total_counts=args.total_counts,
        expected_doublet_rate=args.expected_doublet_rate,
        threshold=args.threshold,
        out_dir=out_root / "03-Cell_filter",
    )
    output_path = out_root / STEP_OUTPUTS["filter"]
    print("filter finished: {}".format(output_path))
    return output_path


def run_cluster(args: argparse.Namespace) -> Path:
    out_root = _resolve_out_root(args)
    _configure_runtime(out_root, args.cache_dir)

    input_h5ad = _as_path(args.input_h5ad) if args.input_h5ad else _default_input_h5ad(out_root, "cluster")
    adata = _read_h5ad(input_h5ad)
    adata = dim_reduction_cluster(
        _adata=adata,
        method=args.method,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        resolution=args.resolution,
        batch_effect=args.batch_effect,
        out_dir=out_root / "04-Dimension_reduction_cluster",
    )
    output_path = out_root / STEP_OUTPUTS["cluster"]
    adata.write(output_path)
    print("cluster finished: {}".format(output_path))
    return output_path


def run_marker(args: argparse.Namespace) -> Path:
    out_root = _resolve_out_root(args)
    _configure_runtime(out_root, args.cache_dir)

    input_h5ad = _as_path(args.input_h5ad) if args.input_h5ad else _default_input_h5ad(out_root, "marker")
    adata = _read_h5ad(input_h5ad)
    adata = find_marker_gene(
        _adata=adata,
        method=args.marker_method,
        out_dir=out_root / "05-Marker_gene_analysis",
        logfc=args.logfc,
        fdr=args.fdr,
        pct=args.pct,
    )
    output_path = out_root / STEP_OUTPUTS["marker"]
    adata.write(output_path)
    print("marker finished: {}".format(output_path))
    return output_path


def run_annotate(args: argparse.Namespace) -> Path:
    out_root = _resolve_out_root(args)
    _configure_runtime(out_root, args.cache_dir)

    input_h5ad = _as_path(args.input_h5ad) if args.input_h5ad else _default_input_h5ad(out_root, "annotate")
    adata = _read_h5ad(input_h5ad)
    celltype_map = _read_celltype_map(args.celltype_map)
    adata, annotation_dict = celltype_annotation(
        _adata=adata,
        database=args.database,
        organism=args.organism,
        out_dir=out_root / "06-Cell_type_annotation",
        cell_type=celltype_map,
        celltypist_model=args.celltypist_model,
    )
    output_path = out_root / STEP_OUTPUTS["annotate"]
    print("annotate finished: {}".format(output_path))
    print("cluster->cell_type mapping: {}".format(annotation_dict))
    return output_path


def run_all(args: argparse.Namespace) -> Path:
    out_root = _resolve_out_root(args)
    _configure_runtime(out_root, args.cache_dir)
    make_report_dir(out_root)

    adata = get_10x_mtx(sample_sheet=args.sample_sheet, out_dir=out_root)
    adata.write(out_root / STEP_OUTPUTS["load"])

    adata = data_qc(adata, out_root / "02-Data_QC")
    adata.write(out_root / STEP_OUTPUTS["qc"])

    adata = cell_filtering(
        _adata=adata,
        min_genes=args.min_genes,
        min_cells=args.min_cells,
        pct_counts_mt=args.pct_counts_mt,
        pct_counts_ribo=args.pct_counts_ribo,
        upper_lim=args.upper_lim,
        total_counts=args.total_counts,
        expected_doublet_rate=args.expected_doublet_rate,
        threshold=args.threshold,
        out_dir=out_root / "03-Cell_filter",
    )

    adata = dim_reduction_cluster(
        _adata=adata,
        method=args.method,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        resolution=args.resolution,
        batch_effect=args.batch_effect,
        out_dir=out_root / "04-Dimension_reduction_cluster",
    )
    adata.write(out_root / STEP_OUTPUTS["cluster"])

    adata = find_marker_gene(
        _adata=adata,
        method=args.marker_method,
        out_dir=out_root / "05-Marker_gene_analysis",
        logfc=args.logfc,
        fdr=args.fdr,
        pct=args.pct,
    )
    adata.write(out_root / STEP_OUTPUTS["marker"])

    celltype_map = _read_celltype_map(args.celltype_map)
    adata, annotation_dict = celltype_annotation(
        _adata=adata,
        database=args.database,
        organism=args.organism,
        out_dir=out_root / "06-Cell_type_annotation",
        cell_type=celltype_map,
        celltypist_model=args.celltypist_model,
    )

    print("all finished: {}".format(out_root / STEP_OUTPUTS["annotate"]))
    print("cluster->cell_type mapping: {}".format(annotation_dict))
    return out_root / STEP_OUTPUTS["annotate"]


def add_common_runtime_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--out-root", default=None, help="Output root. Defaults to the parent of sample_sheet or input_h5ad.")
    parser.add_argument("--cache-dir", default=None, help="Optional scanpy cache directory. Defaults to out_root/.scanpy_cache.")


def add_filter_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--min-cells", type=int, default=3)
    parser.add_argument("--pct-counts-mt", type=float, default=20)
    parser.add_argument("--pct-counts-ribo", type=float, default=40)
    parser.add_argument("--upper-lim", type=float, default=8000)
    parser.add_argument("--total-counts", type=float, default=40000)
    parser.add_argument("--expected-doublet-rate", type=float, default=5)
    parser.add_argument("--threshold", type=float, default=0.25)


def add_cluster_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--method", default="leiden", help="leiden or louvain.")
    parser.add_argument("--n-pcs", type=int, default=20)
    parser.add_argument("--n-neighbors", type=int, default=20)
    parser.add_argument("--resolution", type=float, default=0.5)
    parser.add_argument("--batch-effect", default="harmony", help="harmony, bbknn, or empty.")


def add_marker_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--marker-method", default="wilcoxon")
    parser.add_argument("--logfc", type=float, default=0.5)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--pct", type=float, default=0.25)


def add_annotation_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--database", default="CellTypist", help="CellTypist only.")
    parser.add_argument("--organism", default="human", help="human or mouse.")
    parser.add_argument("--celltypist-model", default=None, help="Optional CellTypist model name. Defaults are chosen from organism.")
    parser.add_argument(
        "--celltype-map",
        default=None,
        help="Optional manual annotation mapping file with cluster and cell_type columns.",
    )


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Server-friendly Scanpy single-cell pipeline entrypoint.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    load_parser = subparsers.add_parser("load", help="Load sample_sheet and merge all samples.")
    add_common_runtime_args(load_parser)
    load_parser.add_argument("--sample-sheet", required=True, help="sample_sheet generated by sample_discovery.py assign.")

    qc_parser = subparsers.add_parser("qc", help="Run QC metrics.")
    add_common_runtime_args(qc_parser)
    qc_parser.add_argument("--input-h5ad", default=None, help="Optional. Defaults to 01-CellRanger/adata_raw.h5ad.")

    filter_parser = subparsers.add_parser("filter", help="Run cell filtering and doublet detection.")
    add_common_runtime_args(filter_parser)
    filter_parser.add_argument("--input-h5ad", default=None, help="Optional. Defaults to 02-Data_QC/adata_qc.h5ad.")
    add_filter_args(filter_parser)

    cluster_parser = subparsers.add_parser("cluster", help="Run dimensionality reduction, clustering, and UMAP.")
    add_common_runtime_args(cluster_parser)
    cluster_parser.add_argument("--input-h5ad", default=None, help="Optional. Defaults to 03-Cell_filter/adata_counts.h5ad.")
    add_cluster_args(cluster_parser)

    marker_parser = subparsers.add_parser("marker", help="Run marker gene analysis.")
    add_common_runtime_args(marker_parser)
    marker_parser.add_argument("--input-h5ad", default=None, help="Optional. Defaults to 04-Dimension_reduction_cluster/adata_cluster.h5ad.")
    add_marker_args(marker_parser)

    annotate_parser = subparsers.add_parser("annotate", help="Run cell type annotation.")
    add_common_runtime_args(annotate_parser)
    annotate_parser.add_argument("--input-h5ad", default=None, help="Optional. Defaults to 05-Marker_gene_analysis/adata_marker.h5ad.")
    add_annotation_args(annotate_parser)

    all_parser = subparsers.add_parser("all", help="Run the full workflow from sample_sheet.")
    add_common_runtime_args(all_parser)
    all_parser.add_argument("--sample-sheet", required=True, help="sample_sheet generated by sample_discovery.py assign.")
    add_filter_args(all_parser)
    add_cluster_args(all_parser)
    add_marker_args(all_parser)
    add_annotation_args(all_parser)

    return parser


def main() -> None:
    parser = build_argument_parser()
    args = parser.parse_args()

    if args.command == "load":
        run_load(args)
    elif args.command == "qc":
        run_qc(args)
    elif args.command == "filter":
        run_filter(args)
    elif args.command == "cluster":
        run_cluster(args)
    elif args.command == "marker":
        run_marker(args)
    elif args.command == "annotate":
        run_annotate(args)
    elif args.command == "all":
        run_all(args)
    else:
        raise ValueError("Unsupported command: {}".format(args.command))


if __name__ == "__main__":
    main()
