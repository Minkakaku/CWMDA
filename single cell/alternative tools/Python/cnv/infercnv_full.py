#!/usr/bin/env python3
"""Integrated infercnvpy pipeline for .h5ad inputs.

This script combines:
1. The infercnvpy workflow used in infercnv.ipynb
2. infercnvpy's CopyKAT wrapper for tumor/normal cell prediction

Examples
--------
Run CopyKAT from a h5ad converted from Seurat:

    python infercnv.py \
        --method copykat \
        --input 04_qs2/05.sc_final.h5ad \
        --species human \
        --group-column cell_type \
        --output-dir 06_CopyKAT

Run infercnvpy's native inferCNV method from a .h5ad file:

    python infercnv.py \
        --method infercnv \
        --input 04_qs2/05.sc_final.h5ad \
        --species mouse \
        --group-column cell_type \
        --reference-types Macrophage,T,iNeutrophil,DC,B,Neutrophil,VSMCs,Endothelial,Plasma,Mesothelial
"""

import argparse
import importlib
import json
import os
import sys
from pathlib import Path


SPECIES_MAP = {
    "human": "hsapiens",
    "mouse": "mmusculus",
}

DEFAULT_FALLBACK_GROUP_COLUMNS = [
    "cell_type",
    "celltype_major",
    "celltype",
    "sub_cell_type",
]
DEFAULT_REFERENCE_TYPES = [
    "T",
    "NK",
    "B",
    "Plasma",
    "Myeloid",
    "MAST",
    "Endothelia",
    "Fibroblast",
    "SMC",
    "Pericyte",
]
DEFAULT_INPUT_H5AD = "04_qs2/05.sc_final.h5ad"
FALLBACK_INPUT_H5ADS = [
    "04_qs2/04.sc_annotated.h5ad",
    "04_qs2/01.sc_annoted.h5ad",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run infercnvpy from a .h5ad input."
    )
    parser.add_argument(
        "--input",
        default=DEFAULT_INPUT_H5AD,
        help="Input .h5ad file path.",
    )
    parser.add_argument(
        "--method",
        choices=["copykat", "infercnv"],
        default="copykat",
        help="CNV method to run. copykat predicts tumor/normal cells; infercnv keeps the original inferCNV workflow.",
    )
    parser.add_argument(
        "--species",
        required=True,
        choices=sorted(SPECIES_MAP),
        help="Species for infercnvpy biomart annotation.",
    )
    parser.add_argument(
        "--output-dir",
        default="06_CopyKAT",
        help="Directory used for all outputs.",
    )
    parser.add_argument(
        "--layer",
        default="counts",
        help="Counts layer name in the .h5ad input.",
    )
    parser.add_argument(
        "--group-column",
        default="cell_type",
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
        default="Epithelia",
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
        "--skip-run",
        action="store_true",
        help="Only prepare data and write outputs before infercnvpy analysis.",
    )
    parser.add_argument(
        "--copykat-key",
        default="copykat",
        help="Key used for CopyKAT results in AnnData.",
    )
    parser.add_argument(
        "--copykat-gene-ids",
        default="S",
        help="Gene annotation style passed to CopyKAT. Use 'S' for gene symbols.",
    )
    parser.add_argument(
        "--copykat-window-size",
        type=int,
        default=25,
        help="Number of genes per sliding window passed to CopyKAT.",
    )
    parser.add_argument(
        "--copykat-sam-name",
        default="copykat",
        help="Sample name prefix used by the R copykat package.",
    )
    parser.add_argument(
        "--copykat-distance",
        default="euclidean",
        help="Distance metric passed to CopyKAT.",
    )
    parser.add_argument(
        "--copykat-segmentation-cut",
        type=float,
        default=0.1,
        help="Segmentation cutoff passed to CopyKAT.",
    )
    parser.add_argument(
        "--copykat-min-genes-chr",
        type=int,
        default=5,
        help="Minimum genes per chromosome passed to CopyKAT.",
    )
    parser.add_argument(
        "--copykat-n-jobs",
        type=int,
        default=1,
        help="Parallel jobs passed to CopyKAT. Use 1 to avoid fork-related memory spikes.",
    )
    parser.add_argument(
        "--copykat-cell-line",
        choices=["yes", "no"],
        default="no",
        help="Whether CopyKAT should run in cell-line mode.",
    )
    parser.add_argument(
        "--copykat-normal-from-reference",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use cells matching --reference-types as known normal cells for CopyKAT when available.",
    )
    parser.add_argument(
        "--include-reference-cells-for-copykat",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Keep sampled reference cells in the CopyKAT input matrix.",
    )
    parser.add_argument(
        "--max-reference-cells",
        type=int,
        default=5000,
        help="Maximum total reference cells kept for CopyKAT. Use 0 to keep all reference cells.",
    )
    parser.add_argument(
        "--max-reference-cells-per-type",
        type=int,
        default=500,
        help="Maximum reference cells kept per reference cell type for CopyKAT. Use 0 to keep all cells per type.",
    )
    parser.add_argument(
        "--reference-random-seed",
        type=int,
        default=2005,
        help="Random seed used when sampling reference cells.",
    )
    return parser.parse_args()


def csv_to_list(raw):
    return [item.strip() for item in raw.split(",") if item.strip()]


def ensure_file(path):
    if not path.exists():
        raise SystemExit("Input file does not exist: %s" % path)
    if path.suffix.lower() != ".h5ad":
        raise SystemExit("Only .h5ad input is supported. Convert the Seurat object before running this script.")


def resolve_input_path(raw_input):
    input_path = Path(raw_input).expanduser().resolve()

    if input_path.exists() or raw_input != DEFAULT_INPUT_H5AD:
        return input_path

    for fallback in FALLBACK_INPUT_H5ADS:
        fallback_path = Path(fallback).expanduser().resolve()
        if fallback_path.exists():
            print(
                "Default input not found (%s); using fallback: %s"
                % (input_path, fallback_path),
                file=sys.stderr,
            )
            return fallback_path

    return input_path


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


def value_counts_as_int_dict(series):
    return {
        str(key): int(value)
        for key, value in series.astype(str).value_counts().sort_index().items()
    }


def format_count_dict(counts, max_items=30):
    items = list(counts.items())
    shown = items[:max_items]
    text = ", ".join("%s=%s" % (key, value) for key, value in shown)
    if len(items) > max_items:
        text += ", ..."
    return text


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
    require_reference=True,
    include_reference_cells=True,
    max_reference_cells=0,
    max_reference_cells_per_type=0,
    reference_random_seed=2005,
):
    numpy = import_or_exit("numpy")

    group_values = adata.obs[group_column].astype(str)
    ref_mask = group_values.isin(reference_types)
    if (require_reference or include_reference_cells) and not ref_mask.any():
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

    ref_mask_to_keep = ref_mask.copy()
    selected_reference_cells = []
    if include_reference_cells and ref_mask.any():
        rng = numpy.random.default_rng(reference_random_seed)

        for ref_type in reference_types:
            type_mask = group_values.eq(ref_type)
            type_cells = numpy.asarray(adata.obs_names[type_mask.to_numpy()].astype(str))

            if len(type_cells) == 0:
                continue

            if max_reference_cells_per_type and len(type_cells) > max_reference_cells_per_type:
                type_cells = rng.choice(
                    type_cells,
                    size=max_reference_cells_per_type,
                    replace=False,
                )

            selected_reference_cells.extend(type_cells.tolist())

        if max_reference_cells and len(selected_reference_cells) > max_reference_cells:
            selected_reference_cells = rng.choice(
                numpy.asarray(selected_reference_cells),
                size=max_reference_cells,
                replace=False,
            ).tolist()

        selected_reference_cells = set(selected_reference_cells)
        if len(selected_reference_cells) == 0:
            raise SystemExit(
                "Reference cell sampling selected 0 cells. Check --reference-types and --group-column."
            )
        ref_mask_to_keep = group_values.index.astype(str).isin(selected_reference_cells)
    elif not include_reference_cells:
        ref_mask_to_keep = ref_mask & False

    use_mask = obs_mask | ref_mask_to_keep
    adata = adata[use_mask].copy()
    adata.obs["infercnv_role"] = numpy.where(
        adata.obs[group_column].astype(str).isin(reference_types),
        "reference",
        "observation",
    )
    adata.uns["reference_cells_selected"] = int(
        (adata.obs["infercnv_role"].astype(str) == "reference").sum()
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


def save_embedding_panel(adata, obsm_key, color_columns, output_path, title_prefix=None):
    if obsm_key not in adata.obsm:
        return None

    numpy = import_or_exit("numpy")
    pandas = import_or_exit("pandas")
    matplotlib = import_or_exit("matplotlib")
    matplotlib.use("Agg")
    plt = import_or_exit("matplotlib.pyplot", pip_name="matplotlib")

    coords = numpy.asarray(adata.obsm[obsm_key])
    if coords.ndim != 2 or coords.shape[1] < 2:
        return None

    color_columns = [col for col in color_columns if col in adata.obs.columns]
    if not color_columns:
        return None

    fig, axes = plt.subplots(
        1,
        len(color_columns),
        figsize=(5.6 * len(color_columns), 5.2),
        squeeze=False,
    )

    x = coords[:, 0]
    y = coords[:, 1]

    for ax, column in zip(axes.ravel(), color_columns):
        values = adata.obs[column]
        title = column if title_prefix is None else "%s: %s" % (title_prefix, column)

        if pandas.api.types.is_numeric_dtype(values):
            numeric_values = pandas.to_numeric(values, errors="coerce")
            sc = ax.scatter(
                x,
                y,
                c=numeric_values,
                s=4,
                linewidths=0,
                cmap="viridis",
                alpha=0.85,
            )
            fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        else:
            cat_values = values.astype("object").where(~values.isna(), "NA").astype(str)
            categories = pandas.Categorical(cat_values).categories.tolist()
            cmap = plt.get_cmap("tab20", max(len(categories), 1))

            for idx, category in enumerate(categories):
                mask = numpy.asarray(cat_values == category)
                ax.scatter(
                    x[mask],
                    y[mask],
                    s=4,
                    linewidths=0,
                    alpha=0.8,
                    color=cmap(idx),
                    label=category,
                )

            if len(categories) <= 24:
                ax.legend(
                    loc="center left",
                    bbox_to_anchor=(1.01, 0.5),
                    frameon=False,
                    fontsize=7,
                    markerscale=2.5,
                )

        ax.set_title(title)
        ax.set_xlabel("%s_1" % obsm_key)
        ax.set_ylabel("%s_2" % obsm_key)
        ax.set_aspect("equal", adjustable="datalim")

    fig.tight_layout()
    fig.savefig(str(output_path), bbox_inches="tight", dpi=300)
    plt.close(fig)
    return output_path


def normalize_copykat_predictions(adata, copykat_key, copykat_work_dir=None):
    pandas = import_or_exit("pandas")

    def clean_label(value):
        if pandas.isna(value):
            return "unknown"
        return str(value).strip()

    def to_status(value):
        clean = clean_label(value).lower()
        compact = (
            clean
            .replace("_", "")
            .replace("-", "")
            .replace(".", "")
            .replace(" ", "")
        )

        tumor_labels = {
            "aneuploid",
            "tumor",
            "tumour",
            "malignant",
            "malign",
            "cancer",
            "cancercell",
            "cancercells",
        }
        normal_labels = {
            "diploid",
            "normal",
            "nontumor",
            "nontumour",
            "benign",
            "normalcell",
            "normalcells",
        }

        if compact in tumor_labels:
            return "Tumor"
        if compact in normal_labels:
            return "Normal"
        return "Unknown"

    def score_series(series):
        labels = series.map(clean_label)
        statuses = labels.map(to_status)
        recognized = int(statuses.isin(["Tumor", "Normal"]).sum())
        non_missing = int(labels.ne("unknown").sum())
        return recognized, non_missing

    def obs_prediction_candidates():
        preferred_columns = [
            copykat_key,
            "%s_prediction" % copykat_key,
            "%s.pred" % copykat_key,
            "copykat.pred",
            "copykat_prediction",
            "prediction",
        ]
        ordered_columns = []

        for column in preferred_columns:
            if column in adata.obs.columns and column not in ordered_columns:
                ordered_columns.append(column)

        for column in adata.obs.columns:
            lower = column.lower()
            if (
                ("copykat" in lower or "prediction" in lower)
                and column not in ordered_columns
                and column not in {"copykat_tumor_status", "copykat_is_tumor"}
            ):
                ordered_columns.append(column)

        return [
            ("obs:%s" % column, adata.obs[column].astype("object"))
            for column in ordered_columns
        ]

    def raw_prediction_candidates():
        if copykat_work_dir is None:
            return []

        raw_dir = Path(copykat_work_dir)
        if not raw_dir.exists():
            return []

        candidates = []
        raw_files = sorted(raw_dir.glob("*prediction*.txt"))
        raw_files.extend(sorted(raw_dir.glob("*prediction*.tsv")))
        raw_files.extend(sorted(raw_dir.glob("*prediction*.csv")))

        obs_index = pandas.Index([str(value) for value in adata.obs_names])

        for raw_file in raw_files:
            try:
                table = pandas.read_csv(str(raw_file), sep=None, engine="python")
            except Exception:
                continue

            if table.empty:
                continue

            lower_to_column = {
                str(column).strip().lower(): column
                for column in table.columns
            }
            cell_column = None
            for candidate in [
                "cell.names",
                "cell.name",
                "cell_names",
                "cell_id",
                "cell",
                "barcode",
            ]:
                if candidate in lower_to_column:
                    cell_column = lower_to_column[candidate]
                    break

            prediction_column = None
            for lower, column in lower_to_column.items():
                compact = lower.replace("_", "").replace(".", "").replace("-", "")
                if compact in {"copykatpred", "prediction", "pred"}:
                    prediction_column = column
                    break

            if cell_column is None or prediction_column is None:
                continue

            table = table.loc[:, [cell_column, prediction_column]].dropna(
                subset=[cell_column]
            )
            table[cell_column] = table[cell_column].astype(str)
            table = table.drop_duplicates(subset=[cell_column], keep="first")

            predictions = table.set_index(cell_column)[prediction_column].reindex(obs_index)
            candidates.append(
                (
                    "raw:%s:%s" % (raw_file.name, prediction_column),
                    predictions.astype("object"),
                )
            )

        return candidates

    candidates = obs_prediction_candidates() + raw_prediction_candidates()

    if not candidates:
        raise SystemExit(
            "CopyKAT did not write a prediction column. Available obs columns: %s"
            % ", ".join(adata.obs.columns)
        )

    prediction_source, prediction_series = max(
        candidates,
        key=lambda item: score_series(item[1]),
    )
    predictions = prediction_series.map(clean_label).astype(str)

    adata.obs["copykat_prediction"] = predictions
    adata.obs["copykat_tumor_status"] = [to_status(value) for value in predictions]
    adata.obs["copykat_is_tumor"] = adata.obs["copykat_tumor_status"].eq("Tumor")
    return prediction_source


def write_copykat_tables(adata, output_dir, group_column):
    pandas = import_or_exit("pandas")

    preferred_columns = [
        group_column,
        "infercnv_role",
        "samples",
        "orig.ident",
        "cell_type",
        "copykat_prediction",
        "copykat_tumor_status",
        "copykat_is_tumor",
    ]
    export_columns = []
    for column in preferred_columns:
        if column in adata.obs.columns and column not in export_columns:
            export_columns.append(column)

    calls = adata.obs.loc[:, export_columns].copy()
    calls.insert(0, "cell_id", adata.obs_names.astype(str))

    calls_path = output_dir / "copykat_cell_calls.csv"
    calls.to_csv(str(calls_path), index=False)

    prediction_counts = (
        adata.obs["copykat_tumor_status"]
        .astype(str)
        .value_counts()
        .rename_axis("copykat_tumor_status")
        .reset_index(name="n_cells")
    )
    prediction_counts_path = output_dir / "copykat_prediction_counts.csv"
    prediction_counts.to_csv(str(prediction_counts_path), index=False)

    output_files = {
        "copykat_cell_calls": str(calls_path),
        "copykat_prediction_counts": str(prediction_counts_path),
    }

    for column in [group_column, "samples", "orig.ident"]:
        if column not in adata.obs.columns:
            continue

        crosstab = pandas.crosstab(
            adata.obs[column].astype(str),
            adata.obs["copykat_tumor_status"].astype(str),
        )
        csv_path = output_dir / ("copykat_status_by_%s.csv" % column.replace(".", "_"))
        crosstab.to_csv(str(csv_path))
        output_files["copykat_status_by_%s" % column.replace(".", "_")] = str(csv_path)

    return output_files


def save_copykat_barplots(adata, output_dir, group_column):
    pandas = import_or_exit("pandas")
    matplotlib = import_or_exit("matplotlib")
    matplotlib.use("Agg")
    plt = import_or_exit("matplotlib.pyplot", pip_name="matplotlib")

    output_files = {}

    for column in [group_column, "samples", "orig.ident"]:
        if column not in adata.obs.columns:
            continue

        crosstab = pandas.crosstab(
            adata.obs[column].astype(str),
            adata.obs["copykat_tumor_status"].astype(str),
        )
        if crosstab.empty:
            continue

        fig, ax = plt.subplots(figsize=(max(6, 0.45 * len(crosstab.index)), 4.8))
        crosstab.plot(kind="bar", stacked=True, ax=ax, width=0.82)
        ax.set_xlabel(column)
        ax.set_ylabel("Cell number")
        ax.legend(title="CopyKAT", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
        ax.tick_params(axis="x", labelrotation=45)
        fig.tight_layout()

        pdf_path = output_dir / ("copykat_status_by_%s.pdf" % column.replace(".", "_"))
        png_path = output_dir / ("copykat_status_by_%s.png" % column.replace(".", "_"))
        fig.savefig(str(pdf_path), bbox_inches="tight")
        fig.savefig(str(png_path), bbox_inches="tight", dpi=300)
        plt.close(fig)

        output_files["copykat_barplot_by_%s_pdf" % column.replace(".", "_")] = str(pdf_path)
        output_files["copykat_barplot_by_%s_png" % column.replace(".", "_")] = str(png_path)

    return output_files


def run_copykat_analysis(
    adata,
    species,
    group_column,
    output_dir,
    copykat_key,
    copykat_gene_ids,
    copykat_window_size,
    copykat_sam_name,
    copykat_distance,
    copykat_segmentation_cut,
    copykat_min_genes_chr,
    copykat_n_jobs,
    copykat_cell_line,
    copykat_normal_from_reference,
):
    scanpy = import_or_exit("scanpy")
    infercnvpy = import_or_exit("infercnvpy")
    numpy = import_or_exit("numpy")
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

    normal_cell_names = ""
    if copykat_normal_from_reference and "infercnv_role" in adata.obs.columns:
        reference_mask = adata.obs["infercnv_role"].astype(str).eq("reference")
        if reference_mask.any():
            normal_cell_names_list = (
                adata.obs_names[reference_mask.to_numpy()]
                .astype(str)
                .tolist()
            )
            try:
                ro = importlib.import_module("rpy2.robjects")
                normal_cell_names = ro.StrVector(normal_cell_names_list)
            except ImportError:
                normal_cell_names = normal_cell_names_list

    copykat_work_dir = output_dir / "copykat_raw"
    copykat_work_dir.mkdir(parents=True, exist_ok=True)

    copykat_kwargs = {
        "gene_ids": copykat_gene_ids,
        "organism": species,
        "layer": "counts" if "counts" in adata.layers else None,
        "key_added": copykat_key,
        "inplace": True,
        "min_genes_chr": copykat_min_genes_chr,
        "window_size": copykat_window_size,
        "segmentation_cut": copykat_segmentation_cut,
        "s_name": copykat_sam_name,
        "distance": copykat_distance,
        "norm_cell_names": normal_cell_names,
        "cell_line": copykat_cell_line,
    }
    if copykat_n_jobs is not None:
        copykat_kwargs["n_jobs"] = copykat_n_jobs

    copykat_kwargs = {
        key: value for key, value in copykat_kwargs.items()
        if value is not None
    }

    cwd = Path.cwd()
    try:
        os.chdir(copykat_work_dir)
        infercnvpy.tl.copykat(adata, **copykat_kwargs)
    finally:
        os.chdir(cwd)

    completed_steps.append("copykat")

    prediction_column = normalize_copykat_predictions(
        adata,
        copykat_key,
        copykat_work_dir=copykat_work_dir,
    )
    completed_steps.append("copykat_prediction")

    output_files.update(write_copykat_tables(adata, output_dir, group_column))
    output_files.update(save_copykat_barplots(adata, output_dir, group_column))

    if "X_umap" in adata.obsm:
        copykat_input_umap = save_embedding_panel(
            adata=adata,
            obsm_key="X_umap",
            color_columns=[
                "copykat_tumor_status",
                "copykat_prediction",
                group_column,
            ],
            output_path=output_dir / "copykat_on_input_umap.pdf",
        )
        if copykat_input_umap is not None:
            output_files["copykat_on_input_umap"] = str(copykat_input_umap)
            completed_steps.append("copykat_on_input_umap")

    try:
        heatmap_path = output_dir / "copykat_chromosome_heatmap_by_status.pdf"
        infercnvpy.pl.chromosome_heatmap(
            adata,
            groupby="copykat_tumor_status",
            use_rep=copykat_key,
            show=False,
        )
        plt.savefig(str(heatmap_path), bbox_inches="tight", dpi=300)
        plt.close()
        output_files["copykat_chromosome_heatmap_by_status"] = str(heatmap_path)
        completed_steps.append("copykat_chromosome_heatmap_by_status")
    except Exception as exc:
        plt.close()
        warnings.append("Skipping CopyKAT chromosome heatmap: %s" % exc)

    try:
        heatmap_group_path = output_dir / "copykat_chromosome_heatmap_by_group.pdf"
        infercnvpy.pl.chromosome_heatmap(
            adata,
            groupby=group_column,
            use_rep=copykat_key,
            show=False,
        )
        plt.savefig(str(heatmap_group_path), bbox_inches="tight", dpi=300)
        plt.close()
        output_files["copykat_chromosome_heatmap_by_group"] = str(heatmap_group_path)
        completed_steps.append("copykat_chromosome_heatmap_by_group")
    except Exception as exc:
        plt.close()
        warnings.append("Skipping CopyKAT group heatmap: %s" % exc)

    if hasattr(infercnvpy.pl, "chromosome_heatmap_summary"):
        try:
            summary_path = output_dir / "copykat_chromosome_heatmap_summary.pdf"
            infercnvpy.pl.chromosome_heatmap_summary(
                adata,
                groupby="copykat_tumor_status",
                use_rep=copykat_key,
                show=False,
            )
            plt.savefig(str(summary_path), bbox_inches="tight", dpi=300)
            plt.close()
            output_files["copykat_chromosome_heatmap_summary"] = str(summary_path)
            completed_steps.append("copykat_chromosome_heatmap_summary")
        except Exception as exc:
            plt.close()
            warnings.append("Skipping CopyKAT summary heatmap: %s" % exc)

    try:
        pca_key = "%s_pca" % copykat_key
        neighbors_key = "%s_neighbors" % copykat_key
        umap_key = "%s_umap" % copykat_key

        infercnvpy.tl.pca(adata, use_rep=copykat_key, key_added=pca_key)
        infercnvpy.pp.neighbors(adata, use_rep=pca_key, key_added=neighbors_key)
        infercnvpy.tl.umap(adata, neighbors_key=neighbors_key, key_added=umap_key)
        completed_steps.extend(["copykat_pca", "copykat_neighbors", "copykat_umap"])

        embedding_key = umap_key
        if embedding_key not in adata.obsm and ("X_%s" % umap_key) in adata.obsm:
            embedding_key = "X_%s" % umap_key

        cnv_umap_path = save_embedding_panel(
            adata=adata,
            obsm_key=embedding_key,
            color_columns=[
                "copykat_tumor_status",
                "copykat_prediction",
                group_column,
            ],
            output_path=output_dir / "copykat_cnv_umap.pdf",
            title_prefix="CNV UMAP",
        )
        if cnv_umap_path is not None:
            output_files["copykat_cnv_umap"] = str(cnv_umap_path)
    except Exception as exc:
        warnings.append("Skipping CopyKAT CNV UMAP: %s" % exc)

    output_files["copykat_raw_dir"] = str(copykat_work_dir)
    output_files["copykat_prediction_column_raw"] = prediction_column

    return adata, completed_steps, warnings, output_files


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
    method,
    species,
    group_column,
    reference_types,
    observation_mode,
    epithelial_label,
    include_reference_cells,
    max_reference_cells,
    max_reference_cells_per_type,
    skipped_run,
    annotation_path,
    prepared_h5ad_path,
    result_h5ad_path,
    completed_steps,
    warnings,
    output_files,
    adata,
):
    counts = value_counts_as_int_dict(adata.obs[group_column])

    summary = {
        "input": str(input_path),
        "method": method,
        "species": species,
        "biomart_species": SPECIES_MAP[species],
        "group_column": group_column,
        "reference_types": list(reference_types),
        "observation_mode": observation_mode,
        "epithelial_label": epithelial_label,
        "include_reference_cells": bool(include_reference_cells),
        "max_reference_cells": int(max_reference_cells),
        "max_reference_cells_per_type": int(max_reference_cells_per_type),
        "reference_cells_selected": int(adata.uns.get("reference_cells_selected", 0)),
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

    if "copykat_tumor_status" in adata.obs.columns:
        summary["copykat_tumor_status_counts"] = value_counts_as_int_dict(
            adata.obs["copykat_tumor_status"]
        )

    summary_path = output_dir / "run_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary_path


def main():
    args = parse_args()

    input_path = resolve_input_path(args.input)
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    ensure_file(input_path)

    fallback_group_columns = csv_to_list(args.fallback_group_columns)
    reference_types = csv_to_list(args.reference_types)

    if not reference_types:
        reference_types = list(DEFAULT_REFERENCE_TYPES)

    adata, used_group_column = load_from_h5ad(
        input_path=input_path,
        counts_layer=args.layer,
        group_column=args.group_column,
        fallback_group_columns=fallback_group_columns,
    )

    adata = prepare_adata(
        adata=adata,
        counts_layer=args.layer,
        gene_id_column=args.gene_id_column,
    )
    group_counts_before_subset = value_counts_as_int_dict(adata.obs[used_group_column])
    print(
        "Available %s counts before subsetting: %s"
        % (used_group_column, format_count_dict(group_counts_before_subset))
    )
    print("Requested reference types: %s" % ", ".join(reference_types))

    include_reference_cells = (
        args.method == "infercnv" or args.include_reference_cells_for_copykat
    )
    adata = subset_for_infercnv(
        adata=adata,
        group_column=used_group_column,
        reference_types=reference_types,
        observation_mode=args.observation_mode,
        epithelial_label=args.epithelial_label,
        require_reference=args.method == "infercnv",
        include_reference_cells=include_reference_cells,
        max_reference_cells=args.max_reference_cells if args.method == "copykat" else 0,
        max_reference_cells_per_type=(
            args.max_reference_cells_per_type if args.method == "copykat" else 0
        ),
        reference_random_seed=args.reference_random_seed,
    )

    print(
        "Cells after subsetting: %d (%s)"
        % (adata.n_obs, value_counts_as_int_dict(adata.obs["infercnv_role"]))
    )
    if include_reference_cells and adata.uns.get("reference_cells_selected", 0) == 0:
        raise SystemExit("No reference cells were kept after subsetting.")

    annotation_path = write_annotation_file(adata, output_dir)

    prepared_h5ad_path = output_dir / "prepared_input.h5ad"
    adata.write(str(prepared_h5ad_path))

    completed_steps = []
    warnings = []
    output_files = {}

    if not args.skip_run:
        if args.method == "copykat":
            adata, completed_steps, warnings, output_files = run_copykat_analysis(
                adata=adata,
                species=args.species,
                group_column=used_group_column,
                output_dir=output_dir,
                copykat_key=args.copykat_key,
                copykat_gene_ids=args.copykat_gene_ids,
                copykat_window_size=args.copykat_window_size,
                copykat_sam_name=args.copykat_sam_name,
                copykat_distance=args.copykat_distance,
                copykat_segmentation_cut=args.copykat_segmentation_cut,
                copykat_min_genes_chr=args.copykat_min_genes_chr,
                copykat_n_jobs=args.copykat_n_jobs,
                copykat_cell_line=args.copykat_cell_line,
                copykat_normal_from_reference=args.copykat_normal_from_reference,
            )
        else:
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
        method=args.method,
        species=args.species,
        group_column=used_group_column,
        reference_types=reference_types,
        observation_mode=args.observation_mode,
        epithelial_label=args.epithelial_label,
        include_reference_cells=include_reference_cells,
        max_reference_cells=args.max_reference_cells if args.method == "copykat" else 0,
        max_reference_cells_per_type=(
            args.max_reference_cells_per_type if args.method == "copykat" else 0
        ),
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
