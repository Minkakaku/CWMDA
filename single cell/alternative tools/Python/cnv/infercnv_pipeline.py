#!/usr/bin/env python3
"""Integrated infercnvpy pipeline for .h5ad and Seurat .rds inputs.

This script combines:
1. The infercnvpy workflow used in infercnv.ipynb
2. The Seurat-to-AnnData preparation idea from Toh5ad.R

Examples
--------
Run from a .h5ad file:

    python3 infercnv.py \
        --input 01.sc_annoted.h5ad \
        --species mouse \
        --group-column cell_type \
        --reference-types Macrophage,T,iNeutrophil,DC,B,Neutrophil,VSMCs,Endothelial,Plasma,Mesothelial

Run from a Seurat .rds file:

    python3 infercnv.py \
        --input 04_rds/01.sc_annoted.rds \
        --species human \
        --group-column celltype_major \
        --reference-types T,B,Myeloid \
        --output-dir infercnv_out
"""

import argparse
import importlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


SPECIES_MAP = {
    "human": "hsapiens",
    "mouse": "mmusculus",
}

DEFAULT_FALLBACK_GROUP_COLUMNS = ["cell_type", "celltype"]
DEFAULT_REFERENCE_TYPES = ["T", "B", "Myeloid"]
DEFAULT_KEEP_REDUCTIONS = ["pca", "umap"]


R_CONVERTER_FUNCTION = r"""
convert_rds_to_files <- function(
    rds_path,
    assay,
    counts_layer,
    group_column,
    fallback_columns,
    keep_reductions,
    outdir
) {
    suppressPackageStartupMessages({
        library(Seurat)
        library(Matrix)
    })

    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    obj <- readRDS(rds_path)
    if (!(assay %in% names(obj@assays))) {
        stop(
            sprintf(
                "Assay '%s' was not found. Available assays: %s",
                assay,
                paste(names(obj@assays), collapse = ", ")
            )
        )
    }

    available_reductions <- character(0)
    if (!is.null(obj@reductions)) {
        available_reductions <- names(obj@reductions)
    }
    keep_reductions <- keep_reductions[keep_reductions != ""]
    keep_reductions <- intersect(keep_reductions, available_reductions)

    obj_rna <- tryCatch(
        DietSeurat(obj, assays = assay, dimreducs = keep_reductions),
        error = function(e) obj
    )
    DefaultAssay(obj_rna) <- assay

    if (exists("JoinLayers")) {
        obj_rna <- tryCatch(
            JoinLayers(obj_rna),
            error = function(e) obj_rna
        )
    }

    get_counts <- function(x, assay, counts_layer) {
        tryCatch(
            LayerData(x, assay = assay, layer = counts_layer),
            error = function(e1) {
                tryCatch(
                    GetAssayData(x, assay = assay, layer = counts_layer),
                    error = function(e2) {
                        GetAssayData(x, assay = assay, slot = counts_layer)
                    }
                )
            }
        )
    }

    mat <- get_counts(obj_rna, assay, counts_layer)
    if (!inherits(mat, "dgCMatrix")) {
        mat <- as(mat, "dgCMatrix")
    }
    mat <- mat[Matrix::rowSums(mat) > 0, , drop = FALSE]

    meta <- obj_rna@meta.data[colnames(mat), , drop = FALSE]
    chosen_col <- group_column
    fallback_columns <- fallback_columns[fallback_columns != ""]
    if (!(chosen_col %in% colnames(meta))) {
        matched <- fallback_columns[fallback_columns %in% colnames(meta)]
        if (length(matched) == 0) {
            stop(
                sprintf(
                    "Could not find metadata column '%s'. Available columns: %s",
                    group_column,
                    paste(colnames(meta), collapse = ", ")
                )
            )
        }
        chosen_col <- matched[[1]]
    }

    meta[[chosen_col]] <- as.character(meta[[chosen_col]])
    meta$cell_id <- rownames(meta)
    meta <- meta[, c("cell_id", setdiff(colnames(meta), "cell_id")), drop = FALSE]

    var_df <- data.frame(
        gene_name = rownames(mat),
        gene_id = rownames(mat),
        stringsAsFactors = FALSE
    )

    Matrix::writeMM(t(mat), file.path(outdir, "counts.mtx"))
    write.table(
        meta,
        file = file.path(outdir, "obs.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    write.table(
        var_df,
        file = file.path(outdir, "var.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    writeLines(chosen_col, file.path(outdir, "group_column_used.txt"))

    if ("pca" %in% keep_reductions) {
        pca_mat <- Embeddings(obj_rna, "pca")[colnames(mat), , drop = FALSE]
        pca_df <- data.frame(cell_id = rownames(pca_mat), pca_mat, check.names = FALSE)
        write.table(
            pca_df,
            file = file.path(outdir, "obsm_pca.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )
    }

    if ("umap" %in% keep_reductions) {
        umap_mat <- Embeddings(obj_rna, "umap")[colnames(mat), , drop = FALSE]
        umap_df <- data.frame(cell_id = rownames(umap_mat), umap_mat, check.names = FALSE)
        write.table(
            umap_df,
            file = file.path(outdir, "obsm_umap.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )
    }
}
"""


R_CONVERTER_SCRIPT = (
    R_CONVERTER_FUNCTION
    + r"""
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
    stop("Expected 7 arguments: input assay layer group_column fallback keep_reductions outdir")
}
fallback_columns <- strsplit(args[[5]], ",", fixed = TRUE)[[1]]
keep_reductions <- strsplit(args[[6]], ",", fixed = TRUE)[[1]]
convert_rds_to_files(
    rds_path = args[[1]],
    assay = args[[2]],
    counts_layer = args[[3]],
    group_column = args[[4]],
    fallback_columns = fallback_columns,
    keep_reductions = keep_reductions,
    outdir = args[[7]]
)
"""
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run infercnvpy from .h5ad or Seurat .rds input."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input file path. Supported suffixes: .h5ad and .rds",
    )
    parser.add_argument(
        "--species",
        required=True,
        choices=sorted(SPECIES_MAP),
        help="Species for infercnvpy biomart annotation.",
    )
    parser.add_argument(
        "--output-dir",
        default="infercnv_out",
        help="Directory used for all outputs.",
    )
    parser.add_argument(
        "--assay",
        default="RNA",
        help="Seurat assay name used for .rds input.",
    )
    parser.add_argument(
        "--layer",
        default="counts",
        help="Counts layer name for both .h5ad and .rds input.",
    )
    parser.add_argument(
        "--group-column",
        default="celltype_major",
        help="Primary metadata column used to define reference and observation cells.",
    )
    parser.add_argument(
        "--fallback-group-columns",
        default=",".join(DEFAULT_FALLBACK_GROUP_COLUMNS),
        help="Comma-separated fallback metadata columns.",
    )
    parser.add_argument(
        "--reference-types",
        default=",".join(DEFAULT_REFERENCE_TYPES),
        help="Comma-separated normal/reference cell types.",
    )
    parser.add_argument(
        "--observation-mode",
        choices=["all_non_reference", "epithelial_only"],
        default="all_non_reference",
        help="How to define observation cells.",
    )
    parser.add_argument(
        "--epithelial-label",
        default="Epithelial",
        help="Observation label used when observation-mode is epithelial_only.",
    )
    parser.add_argument(
        "--window-size",
        type=int,
        default=250,
        help="Window size passed to infercnvpy.",
    )
    parser.add_argument(
        "--gene-id-column",
        default="gene_ids",
        help="Column name in adata.var used for biomart gene mapping.",
    )
    parser.add_argument(
        "--biomart-gene-id",
        default="external_gene_name",
        help="Gene identifier type used in biomart.",
    )
    parser.add_argument(
        "--heatmap-groupby",
        default=None,
        help="Column used for the main chromosome heatmap. Defaults to the chosen group column.",
    )
    parser.add_argument(
        "--keep-reductions",
        default=",".join(DEFAULT_KEEP_REDUCTIONS),
        help="Comma-separated reductions to preserve when converting .rds input.",
    )
    parser.add_argument(
        "--r-loader",
        choices=["auto", "rscript", "rpy2"],
        default="auto",
        help="Method used to read Seurat .rds input.",
    )
    parser.add_argument(
        "--rscript-bin",
        default="Rscript",
        help="Rscript executable used when r-loader is rscript or auto.",
    )
    parser.add_argument(
        "--skip-run",
        action="store_true",
        help="Only prepare data and write outputs before infercnvpy analysis.",
    )
    return parser.parse_args()


def csv_to_list(raw):
    return [item.strip() for item in raw.split(",") if item.strip()]


def ensure_file(path):
    if not path.exists():
        raise SystemExit("Input file does not exist: %s" % path)
    if path.suffix.lower() not in [".h5ad", ".rds"]:
        raise SystemExit("Only .h5ad and .rds inputs are supported.")


def optional_import(module_name):
    try:
        return importlib.import_module(module_name)
    except ImportError:
        return None


def import_or_exit(module_name, pip_name=None):
    try:
        return importlib.import_module(module_name)
    except ImportError as exc:
        package_name = pip_name or module_name
        raise SystemExit(
            "Missing Python dependency: %s\nOriginal error: %s"
            % (package_name, exc)
        )


def pick_group_column(columns, preferred, fallbacks):
    if preferred in columns:
        return preferred
    for name in fallbacks:
        if name in columns:
            return name
    raise SystemExit(
        "Could not find a usable group column. Tried: %s. Available columns: %s"
        % ([preferred] + list(fallbacks), ", ".join(columns))
    )


def choose_r_loader(preferred_loader, rscript_bin):
    has_rscript = shutil.which(rscript_bin) is not None
    has_rpy2 = optional_import("rpy2.robjects") is not None

    if preferred_loader == "rscript":
        if not has_rscript:
            raise SystemExit("Rscript was not found in PATH: %s" % rscript_bin)
        return "rscript"

    if preferred_loader == "rpy2":
        if not has_rpy2:
            raise SystemExit("rpy2 is not installed in the current Python environment.")
        return "rpy2"

    if has_rscript:
        return "rscript"
    if has_rpy2:
        return "rpy2"

    raise SystemExit(
        "Reading .rds input requires either Rscript in PATH or the Python package rpy2."
    )


def run_r_converter_with_rscript(
    input_path,
    assay,
    counts_layer,
    group_column,
    fallback_group_columns,
    keep_reductions,
    temp_dir,
    rscript_bin,
):
    script_path = temp_dir / "convert_rds_to_files.R"
    script_path.write_text(R_CONVERTER_SCRIPT, encoding="utf-8")

    command = [
        rscript_bin,
        str(script_path),
        str(input_path),
        assay,
        counts_layer,
        group_column,
        ",".join(fallback_group_columns),
        ",".join(keep_reductions),
        str(temp_dir),
    ]

    process = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    if process.returncode != 0:
        raise SystemExit(
            "Rscript conversion failed.\nSTDOUT:\n%s\nSTDERR:\n%s"
            % (process.stdout, process.stderr)
        )


def run_r_converter_with_rpy2(
    input_path,
    assay,
    counts_layer,
    group_column,
    fallback_group_columns,
    keep_reductions,
    temp_dir,
):
    ro = import_or_exit("rpy2.robjects", pip_name="rpy2")
    function = ro.r(R_CONVERTER_FUNCTION + "\nconvert_rds_to_files")

    try:
        function(
            str(input_path),
            assay,
            counts_layer,
            group_column,
            ro.StrVector(list(fallback_group_columns)),
            ro.StrVector(list(keep_reductions)),
            str(temp_dir),
        )
    except Exception as exc:
        raise SystemExit("rpy2 conversion failed: %s" % exc)


def read_r_conversion_outputs(temp_dir, gene_id_column):
    pd = import_or_exit("pandas")
    anndata = import_or_exit("anndata")
    scipy_io = import_or_exit("scipy.io", pip_name="scipy")

    counts_path = temp_dir / "counts.mtx"
    obs_path = temp_dir / "obs.tsv"
    var_path = temp_dir / "var.tsv"
    group_col_path = temp_dir / "group_column_used.txt"

    if not counts_path.exists():
        raise SystemExit("R conversion did not produce counts.mtx")
    if not obs_path.exists():
        raise SystemExit("R conversion did not produce obs.tsv")
    if not var_path.exists():
        raise SystemExit("R conversion did not produce var.tsv")
    if not group_col_path.exists():
        raise SystemExit("R conversion did not produce group_column_used.txt")

    matrix = scipy_io.mmread(str(counts_path)).tocsr()
    obs = pd.read_csv(str(obs_path), sep="\t")
    var = pd.read_csv(str(var_path), sep="\t")

    if "cell_id" not in obs.columns:
        raise SystemExit("obs.tsv is missing the required cell_id column.")
    if "gene_name" not in var.columns:
        raise SystemExit("var.tsv is missing the required gene_name column.")

    obs = obs.set_index("cell_id")
    var = var.set_index("gene_name")
    adata = anndata.AnnData(X=matrix, obs=obs, var=var)

    if gene_id_column not in adata.var.columns:
        if "gene_id" in adata.var.columns:
            adata.var[gene_id_column] = adata.var["gene_id"].astype(str)
        else:
            adata.var[gene_id_column] = adata.var_names.astype(str)

    pca_path = temp_dir / "obsm_pca.tsv"
    umap_path = temp_dir / "obsm_umap.tsv"

    if pca_path.exists():
        pca_df = pd.read_csv(str(pca_path), sep="\t").set_index("cell_id")
        pca_df = pca_df.loc[adata.obs_names]
        adata.obsm["X_pca"] = pca_df.to_numpy()

    if umap_path.exists():
        umap_df = pd.read_csv(str(umap_path), sep="\t").set_index("cell_id")
        umap_df = umap_df.loc[adata.obs_names]
        adata.obsm["X_umap"] = umap_df.to_numpy()

    used_group_column = group_col_path.read_text(encoding="utf-8").strip()
    return adata, used_group_column


def load_from_rds(
    input_path,
    assay,
    counts_layer,
    group_column,
    fallback_group_columns,
    keep_reductions,
    r_loader,
    rscript_bin,
    gene_id_column,
):
    selected_loader = choose_r_loader(r_loader, rscript_bin)

    with tempfile.TemporaryDirectory(prefix="infercnv_rds_") as temp_name:
        temp_dir = Path(temp_name)
        if selected_loader == "rscript":
            run_r_converter_with_rscript(
                input_path=input_path,
                assay=assay,
                counts_layer=counts_layer,
                group_column=group_column,
                fallback_group_columns=fallback_group_columns,
                keep_reductions=keep_reductions,
                temp_dir=temp_dir,
                rscript_bin=rscript_bin,
            )
        else:
            run_r_converter_with_rpy2(
                input_path=input_path,
                assay=assay,
                counts_layer=counts_layer,
                group_column=group_column,
                fallback_group_columns=fallback_group_columns,
                keep_reductions=keep_reductions,
                temp_dir=temp_dir,
            )

        adata, used_group_column = read_r_conversion_outputs(
            temp_dir=temp_dir,
            gene_id_column=gene_id_column,
        )

    return adata, used_group_column, selected_loader


def load_from_h5ad(input_path, counts_layer, group_column, fallback_group_columns):
    scanpy = import_or_exit("scanpy")

    adata = scanpy.read_h5ad(str(input_path))
    used_group_column = pick_group_column(
        list(adata.obs.columns),
        group_column,
        fallback_group_columns,
    )
    if counts_layer and counts_layer in adata.layers:
        adata.X = adata.layers[counts_layer].copy()

    return adata, used_group_column


def prepare_adata(adata, counts_layer, gene_id_column):
    numpy = import_or_exit("numpy")

    if counts_layer and counts_layer in adata.layers:
        adata.X = adata.layers[counts_layer].copy()

    original_var_names = [str(name) for name in adata.var_names]
    adata.var_names = original_var_names

    if gene_id_column in adata.var.columns:
        raw_gene_ids = [str(value) for value in adata.var[gene_id_column].tolist()]
    else:
        raw_gene_ids = list(original_var_names)

    adata.var_names_make_unique()

    cleaned_gene_ids = []
    for idx, gene_id in enumerate(raw_gene_ids):
        if gene_id is None:
            cleaned_gene_ids.append(original_var_names[idx])
            continue
        gene_id = str(gene_id)
        if gene_id == "" or gene_id.lower() == "nan":
            cleaned_gene_ids.append(original_var_names[idx])
        else:
            cleaned_gene_ids.append(gene_id)
    adata.var[gene_id_column] = cleaned_gene_ids

    gene_mask = numpy.asarray(adata.X.sum(axis=0)).ravel() > 0
    adata = adata[:, gene_mask].copy()

    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    return adata


def subset_for_infercnv(
    adata,
    group_column,
    reference_types,
    observation_mode,
    epithelial_label,
):
    numpy = import_or_exit("numpy")

    group_values = adata.obs[group_column].astype(str)
    ref_mask = group_values.isin(reference_types)
    if not ref_mask.any():
        raise SystemExit(
            "No reference cells found. Requested: %s. Available: %s"
            % (list(reference_types), sorted(group_values.unique().tolist()))
        )

    if observation_mode == "epithelial_only":
        obs_mask = group_values.eq(epithelial_label)
        if not obs_mask.any():
            raise SystemExit(
                "No epithelial observation cells found for label '%s'. Available: %s"
                % (epithelial_label, sorted(group_values.unique().tolist()))
            )
    else:
        obs_mask = ~ref_mask
        if not obs_mask.any():
            raise SystemExit(
                "No observation cells remained after excluding reference cells."
            )

    use_mask = ref_mask | obs_mask
    adata = adata[use_mask].copy()
    adata.obs["infercnv_role"] = numpy.where(
        adata.obs[group_column].astype(str).isin(reference_types),
        "reference",
        "observation",
    )
    return adata


def write_annotation_file(adata, output_dir):
    pandas = import_or_exit("pandas")

    annotation = pandas.DataFrame(
        {
            "cell": adata.obs_names.astype(str),
            "group": adata.obs["infercnv_role"].astype(str).tolist(),
        }
    )
    annotation_path = output_dir / "cell_annotations.txt"
    annotation.to_csv(
        str(annotation_path),
        sep="\t",
        header=False,
        index=False,
    )
    return annotation_path


def save_input_umap_if_available(adata, group_column, output_dir, plt, scanpy):
    if "X_umap" not in adata.obsm:
        return None

    if "X_umap_input" not in adata.obsm:
        adata.obsm["X_umap_input"] = adata.obsm["X_umap"].copy()

    plot_path = output_dir / "input_umap_by_group.pdf"
    try:
        scanpy.pl.umap(adata, color=group_column, show=False)
        plt.savefig(str(plot_path), bbox_inches="tight", dpi=300)
        plt.close()
        return plot_path
    except Exception:
        plt.close()
        return None


def run_infercnv_analysis(
    adata,
    species,
    gene_id_column,
    biomart_gene_id,
    group_column,
    window_size,
    heatmap_groupby,
    output_dir,
):
    scanpy = import_or_exit("scanpy")
    infercnvpy = import_or_exit("infercnvpy")
    matplotlib = import_or_exit("matplotlib")
    matplotlib.use("Agg")
    plt = import_or_exit("matplotlib.pyplot", pip_name="matplotlib")

    completed_steps = []
    warnings = []
    output_files = {}

    input_umap_path = save_input_umap_if_available(
        adata=adata,
        group_column=group_column,
        output_dir=output_dir,
        plt=plt,
        scanpy=scanpy,
    )
    if input_umap_path is not None:
        output_files["input_umap_by_group"] = str(input_umap_path)

    infercnvpy.io.genomic_position_from_biomart(
        adata=adata,
        adata_gene_id=gene_id_column,
        biomart_gene_id=biomart_gene_id,
        species=SPECIES_MAP[species],
    )
    completed_steps.append("genomic_position_from_biomart")

    infercnvpy.tl.infercnv(
        adata,
        reference_key="infercnv_role",
        reference_cat=["reference"],
        window_size=window_size,
    )
    completed_steps.append("infercnv")

    main_groupby = heatmap_groupby or group_column
    main_heatmap_path = output_dir / "cnv_chromosome_heatmap_by_group.pdf"
    infercnvpy.pl.chromosome_heatmap(adata, groupby=main_groupby, show=False)
    plt.savefig(str(main_heatmap_path), bbox_inches="tight", dpi=300)
    plt.close()
    output_files["chromosome_heatmap_by_group"] = str(main_heatmap_path)
    completed_steps.append("chromosome_heatmap_by_group")

    neighbors_ready = False
    try:
        infercnvpy.tl.pca(adata)
        infercnvpy.pp.neighbors(adata)
        completed_steps.extend(["cnv_pca", "cnv_neighbors"])
        neighbors_ready = True
    except Exception as exc:
        warnings.append("Skipping pca/neighbors: %s" % exc)

    if neighbors_ready:
        try:
            infercnvpy.tl.leiden(adata)
            completed_steps.append("cnv_leiden")

            leiden_heatmap_path = output_dir / "cnv_chromosome_heatmap_by_leiden.pdf"
            infercnvpy.pl.chromosome_heatmap(
                adata,
                groupby="cnv_leiden",
                dendrogram=True,
                show=False,
            )
            plt.savefig(str(leiden_heatmap_path), bbox_inches="tight", dpi=300)
            plt.close()
            output_files["chromosome_heatmap_by_leiden"] = str(leiden_heatmap_path)
            completed_steps.append("chromosome_heatmap_by_leiden")
        except Exception as exc:
            warnings.append("Skipping leiden heatmap: %s" % exc)

        try:
            infercnvpy.tl.umap(adata)
            infercnvpy.tl.cnv_score(adata)
            completed_steps.extend(["cnv_umap", "cnv_score"])

            if "X_umap" in adata.obsm:
                cnv_umap_path = output_dir / "cnv_umap.pdf"
                fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
                scanpy.pl.umap(adata, color="cnv_score", ax=axes[0], show=False)
                scanpy.pl.umap(adata, color=group_column, ax=axes[1], show=False)
                plt.tight_layout()
                plt.savefig(str(cnv_umap_path), bbox_inches="tight", dpi=300)
                plt.close()
                output_files["cnv_umap"] = str(cnv_umap_path)
        except Exception as exc:
            warnings.append("Skipping cnv UMAP/cnv_score plotting: %s" % exc)

    return adata, completed_steps, warnings, output_files


def write_summary(
    output_dir,
    input_path,
    species,
    group_column,
    reference_types,
    observation_mode,
    epithelial_label,
    loader_used,
    skipped_run,
    annotation_path,
    prepared_h5ad_path,
    result_h5ad_path,
    completed_steps,
    warnings,
    output_files,
    adata,
):
    counts = (
        adata.obs[group_column]
        .astype(str)
        .value_counts()
        .sort_index()
        .to_dict()
    )

    summary = {
        "input": str(input_path),
        "species": species,
        "biomart_species": SPECIES_MAP[species],
        "group_column": group_column,
        "reference_types": list(reference_types),
        "observation_mode": observation_mode,
        "epithelial_label": epithelial_label,
        "loader_used": loader_used,
        "skip_run": skipped_run,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "cells_per_group": counts,
        "annotation_file": str(annotation_path),
        "prepared_h5ad_file": str(prepared_h5ad_path),
        "result_h5ad_file": str(result_h5ad_path),
        "completed_steps": list(completed_steps),
        "warnings": list(warnings),
        "output_files": dict(output_files),
    }

    summary_path = output_dir / "run_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary_path


def main():
    args = parse_args()

    input_path = Path(args.input).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    ensure_file(input_path)

    fallback_group_columns = csv_to_list(args.fallback_group_columns)
    reference_types = csv_to_list(args.reference_types)
    keep_reductions = csv_to_list(args.keep_reductions)

    if not reference_types:
        reference_types = list(DEFAULT_REFERENCE_TYPES)

    if input_path.suffix.lower() == ".h5ad":
        adata, used_group_column = load_from_h5ad(
            input_path=input_path,
            counts_layer=args.layer,
            group_column=args.group_column,
            fallback_group_columns=fallback_group_columns,
        )
        loader_used = "h5ad"
    else:
        adata, used_group_column, loader_used = load_from_rds(
            input_path=input_path,
            assay=args.assay,
            counts_layer=args.layer,
            group_column=args.group_column,
            fallback_group_columns=fallback_group_columns,
            keep_reductions=keep_reductions,
            r_loader=args.r_loader,
            rscript_bin=args.rscript_bin,
            gene_id_column=args.gene_id_column,
        )

    adata = prepare_adata(
        adata=adata,
        counts_layer=args.layer,
        gene_id_column=args.gene_id_column,
    )
    adata = subset_for_infercnv(
        adata=adata,
        group_column=used_group_column,
        reference_types=reference_types,
        observation_mode=args.observation_mode,
        epithelial_label=args.epithelial_label,
    )

    annotation_path = write_annotation_file(adata, output_dir)

    prepared_h5ad_path = output_dir / "prepared_input.h5ad"
    adata.write(str(prepared_h5ad_path))

    completed_steps = []
    warnings = []
    output_files = {}

    if not args.skip_run:
        adata, completed_steps, warnings, output_files = run_infercnv_analysis(
            adata=adata,
            species=args.species,
            gene_id_column=args.gene_id_column,
            biomart_gene_id=args.biomart_gene_id,
            group_column=used_group_column,
            window_size=args.window_size,
            heatmap_groupby=args.heatmap_groupby,
            output_dir=output_dir,
        )

    result_h5ad_path = output_dir / "infercnv_result.h5ad"
    adata.write(str(result_h5ad_path))

    summary_path = write_summary(
        output_dir=output_dir,
        input_path=input_path,
        species=args.species,
        group_column=used_group_column,
        reference_types=reference_types,
        observation_mode=args.observation_mode,
        epithelial_label=args.epithelial_label,
        loader_used=loader_used,
        skipped_run=args.skip_run,
        annotation_path=annotation_path,
        prepared_h5ad_path=prepared_h5ad_path,
        result_h5ad_path=result_h5ad_path,
        completed_steps=completed_steps,
        warnings=warnings,
        output_files=output_files,
        adata=adata,
    )

    print("Finished. Outputs written to: %s" % output_dir)
    print("Loader used: %s" % loader_used)
    print("Annotation file: %s" % annotation_path)
    print("Prepared h5ad: %s" % prepared_h5ad_path)
    print("Result h5ad: %s" % result_h5ad_path)
    print("Summary: %s" % summary_path)
    if warnings:
        print("Warnings:", file=sys.stderr)
        for warning in warnings:
            print("- %s" % warning, file=sys.stderr)


if __name__ == "__main__":
    main()
