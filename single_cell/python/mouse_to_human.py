from __future__ import annotations

import argparse

import scanpy as sc

from io_utils import as_path, ensure_dir, read_table


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert mouse genes in an h5ad object to human symbols.")
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--ortholog-table", required=True, help="TSV with mouse gene as index and human_gene column.")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--output-name", default="human.h5ad")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    adata = sc.read_h5ad(str(as_path(args.input_h5ad)))
    mapping = read_table(args.ortholog_table)
    if "human_gene" not in mapping.columns:
        raise KeyError("ortholog table must contain a human_gene column.")

    mapping = mapping[~mapping["human_gene"].duplicated(keep="first")].copy()
    valid_genes = [gene for gene in adata.raw.var_names if gene in mapping.index] if adata.raw is not None else [gene for gene in adata.var_names if gene in mapping.index]
    converted = adata.raw.to_adata() if adata.raw is not None else adata.copy()
    converted = converted[:, valid_genes].copy()
    merged = converted.var.merge(mapping, left_index=True, right_index=True, how="left", validate="1:1")
    converted.var = merged.set_index("human_gene")
    converted.raw = converted

    output_dir = ensure_dir(args.output_dir)
    converted.write_h5ad(output_dir / args.output_name)
    print("Converted h5ad written to: {}".format(output_dir / args.output_name))


if __name__ == "__main__":
    main()
