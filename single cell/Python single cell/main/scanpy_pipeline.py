from __future__ import annotations

import os
import shutil
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from PIL import Image
from matplotlib import font_manager as fm
from matplotlib.pyplot import rc_context


REPORT_SUBDIRS = [
    "01-CellRanger",
    "02-Data_QC",
    "03-Cell_filter",
    "04-Dimension_reduction_cluster",
    "05-Marker_gene_analysis",
    "06-Cell_type_annotation",
]
OUTS_DIR_NAMES = ("outs", "Outs")
OPTIONAL_CELLRANGER_FILES = (
    "metrics_summary.csv",
    "cloupe.cloupe",
    "web_summary.html",
)

warnings.filterwarnings("ignore")
sc.settings.verbosity = 3
sc.settings.set_figure_params(figsize=(5, 5), facecolor="white")


PathLike = Union[str, os.PathLike]


def configure_scanpy(cache_dir: Optional[PathLike] = None) -> None:
    if cache_dir:
        cache_path = Path(cache_dir).expanduser().resolve()
        cache_path.mkdir(parents=True, exist_ok=True)
        sc.settings.cachedir = str(cache_path)


def _as_path(path_value: PathLike) -> Path:
    return Path(path_value).expanduser().resolve()


def _ensure_dir(path_value: PathLike) -> Path:
    path = _as_path(path_value)
    path.mkdir(parents=True, exist_ok=True)
    return path


def _normalize_method_name(value: Optional[str]) -> str:
    return (value or "").replace("*", "").strip().lower()


def _read_table(table_path: PathLike) -> pd.DataFrame:
    path = _as_path(table_path)
    suffix = path.suffix.lower()
    sep = "\t" if suffix in {".tsv", ".txt", ".xls"} else ","
    return pd.read_csv(path, sep=sep)


def make_report_dir(out_dir: PathLike) -> str:
    out_path = _ensure_dir(out_dir)
    for subdir in REPORT_SUBDIRS:
        (out_path / subdir).mkdir(parents=True, exist_ok=True)
    return str(out_path)


def _coerce_sample_table(
    sample_sheet: Union[pd.DataFrame, PathLike]
) -> pd.DataFrame:
    if isinstance(sample_sheet, pd.DataFrame):
        sample_df = sample_sheet.copy()
    else:
        sample_df = _read_table(sample_sheet)

    required_columns = {"sample", "group", "matrix_dir", "outs_dir"}
    missing_columns = required_columns.difference(sample_df.columns)
    if missing_columns:
        raise ValueError(
            "sample_sheet is missing required columns: {}".format(", ".join(sorted(missing_columns)))
        )

    sample_df = sample_df.copy()
    sample_df["sample"] = sample_df["sample"].astype(str)
    sample_df["group"] = sample_df["group"].astype(str)

    duplicate_samples = sample_df["sample"][sample_df["sample"].duplicated()].unique()
    if len(duplicate_samples) > 0:
        raise ValueError(
            "sample_sheet contains duplicated sample names: {}".format(", ".join(map(str, duplicate_samples)))
        )

    return sample_df


def _copy_cellranger_outputs(outs_dir: Path, destination_dir: Path) -> None:
    destination_dir.mkdir(parents=True, exist_ok=True)

    for filename in OPTIONAL_CELLRANGER_FILES:
        source_file = outs_dir / filename
        if source_file.exists():
            shutil.copy2(source_file, destination_dir / filename)

    matrix_dir = outs_dir / "filtered_feature_bc_matrix"
    if matrix_dir.exists():
        matrix_target = destination_dir / "filtered_feature_bc_matrix"
        if matrix_target.exists():
            shutil.rmtree(matrix_target)
        shutil.copytree(matrix_dir, matrix_target, dirs_exist_ok=True)

    matrix_h5 = outs_dir / "filtered_feature_bc_matrix.h5"
    if matrix_h5.exists():
        shutil.copy2(matrix_h5, destination_dir / "filtered_feature_bc_matrix.h5")


def get_10x_mtx(
    sample_sheet: Union[pd.DataFrame, PathLike],
    out_dir: PathLike,
) -> sc.AnnData:
    sample_df = _coerce_sample_table(sample_sheet)
    report_root = _as_path(out_dir)
    adatas_list = []

    for row in sample_df.itertuples(index=False):
        matrix_dir = _as_path(row.matrix_dir)
        outs_dir = _as_path(row.outs_dir)
        if not matrix_dir.exists():
            raise FileNotFoundError("Matrix directory not found: {}".format(matrix_dir))
        if not outs_dir.exists():
            raise FileNotFoundError("outs/Outs directory not found: {}".format(outs_dir))

        matrix_type = getattr(row, "matrix_type", None)
        if matrix_dir.is_file() or matrix_type == "h5":
            data = sc.read_10x_h5(str(matrix_dir))
        else:
            data = sc.read_10x_mtx(
                str(matrix_dir),
                var_names="gene_symbols",
                cache=True,
            )
        data.var_names_make_unique()
        data.obs["sample"] = row.sample
        data.obs["group"] = row.group
        adatas_list.append(data)

        sample_report_dir = report_root / "01-CellRanger" / row.sample
        _copy_cellranger_outputs(outs_dir, sample_report_dir)

    if not adatas_list:
        raise ValueError("No samples were loaded. Please check sample_sheet and input paths.")

    if len(adatas_list) == 1:
        adata = adatas_list[0].copy()
    else:
        adata = adatas_list[0].concatenate(*adatas_list[1:], index_unique=None)
        if "batch" in adata.obs:
            del adata.obs["batch"]

    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    return adata


def save_fig(
    figure_or_plot: object,
    out_dir: PathLike,
    file_name: str,
    fig_size: Optional[Mapping[str, float]] = None,
) -> None:
    out_path = _ensure_dir(out_dir)

    if hasattr(figure_or_plot, "savefig") and hasattr(figure_or_plot, "tight_layout"):
        figure = figure_or_plot
    else:
        figure = plt.gcf()

    if fig_size:
        width = fig_size.get("width", 5)
        height = fig_size.get("height", 5)
        figure.set_size_inches(width, height)

    figure.tight_layout()
    png_path = out_path / "{}.png".format(file_name)
    pdf_path = out_path / "{}.pdf".format(file_name)
    figure.savefig(png_path, bbox_inches="tight", pad_inches=0.1)
    figure.savefig(pdf_path, bbox_inches="tight", pad_inches=0.1)
    plt.close(figure)


def _write_obs_table(adata: sc.AnnData, output_path: PathLike) -> None:
    output = _as_path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(adata.obs).to_csv(output, sep="\t")


def _plot_qc_panels(adata: sc.AnnData, out_dir: PathLike, prefix: str) -> None:
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"],
        jitter=0,
        groupby="sample",
        rotation=45,
        multi_panel=True,
        show=False,
    )
    axes = plt.gcf().axes
    for axis in axes:
        axis.set_xticklabels(axis.get_xticklabels(), horizontalalignment="right")
    save_fig(plt.gcf(), out_dir, file_name="{}_nCount_nFeature_voilon".format(prefix), fig_size={"width": 25, "height": 10})

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))
    sc.pl.scatter(
        adata,
        x="total_counts",
        y="pct_counts_mt",
        color="sample",
        legend_loc="none",
        ax=ax1,
        show=False,
    )
    sc.pl.scatter(
        adata,
        x="total_counts",
        y="n_genes_by_counts",
        color="sample",
        ax=ax2,
        show=False,
    )
    save_fig(fig, out_dir, file_name="{}_nCount_nFeature_scatter".format(prefix))


def data_qc(_adata: sc.AnnData, out_dir: PathLike) -> sc.AnnData:
    out_path = _ensure_dir(out_dir)
    _adata.var["mt"] = _adata.var_names.str.upper().str.startswith("MT-")
    _adata.var["ribo"] = _adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(
        _adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    _plot_qc_panels(_adata, out_path, prefix="orig")
    _write_obs_table(_adata, out_path / "all_cell_features_raw.xls")

    for sample_name in _adata.obs["sample"].drop_duplicates().tolist():
        sample_dir = out_path / sample_name
        sample_subset = _adata[_adata.obs["sample"] == sample_name].copy()
        _plot_qc_panels(sample_subset, sample_dir, prefix="orig")
        _write_obs_table(sample_subset, sample_dir / "cell_features_raw.xls")

    return _adata


def cell_filtering(
    _adata: sc.AnnData,
    min_genes: int,
    min_cells: int,
    pct_counts_mt: float,
    pct_counts_ribo: float,
    upper_lim: float,
    total_counts: float,
    expected_doublet_rate: float,
    threshold: float,
    out_dir: PathLike,
) -> sc.AnnData:
    out_path = _ensure_dir(out_dir)

    sc.pp.filter_cells(_adata, min_genes=min_genes)
    sc.pp.filter_genes(_adata, min_cells=min_cells)
    _adata = _adata[_adata.obs["n_genes_by_counts"] < upper_lim].copy()
    _adata = _adata[_adata.obs["pct_counts_mt"] < pct_counts_mt, :].copy()
    _adata = _adata[_adata.obs["pct_counts_ribo"] < pct_counts_ribo, :].copy()
    _adata = _adata[_adata.obs["total_counts"] < total_counts, :].copy()

    if _adata.n_obs == 0:
        raise ValueError("No cells remain after filtering. Please loosen the filtering thresholds.")

    sc.external.pp.scrublet(
        _adata,
        expected_doublet_rate=expected_doublet_rate / 100.0,
        threshold=threshold,
        batch_key="sample",
        random_state=1234,
    )

    samples = _adata.obs["sample"].drop_duplicates().tolist()
    if len(samples) > 1:
        sc.external.pl.scrublet_score_distribution(_adata)
        save_fig(plt.gcf(), out_path, file_name="scrublet_score_distribution", fig_size={"width": 8, "height": 6})
    _adata = _adata[_adata.obs["predicted_doublet"] == False, :].copy()  # noqa: E712

    _write_obs_table(_adata, out_path / "all_cell_features_filtered.xls")
    _plot_qc_panels(_adata, out_path, prefix="filtered")

    for sample_name in _adata.obs["sample"].drop_duplicates().tolist():
        sample_dir = out_path / sample_name
        sample_subset = _adata[_adata.obs["sample"] == sample_name].copy()
        _plot_qc_panels(sample_subset, sample_dir, prefix="filtered")
        _write_obs_table(sample_subset, sample_dir / "cell_features_filtered.xls")

    _adata.write(out_path / "adata_counts.h5ad")
    return _adata


def _save_embedding_split_panels(
    adata: sc.AnnData,
    category: str,
    color: str,
    out_dir: PathLike,
    file_name: str,
) -> None:
    values = adata.obs[category].drop_duplicates().tolist()
    if len(values) <= 1:
        return

    fig, axes = plt.subplots(nrows=1, ncols=len(values), figsize=(len(values) * 5, 5))
    if len(values) == 2:
        axes = list(axes)

    for index, value in enumerate(values):
        axis = axes[index]
        sc.pl.umap(
            adata[adata.obs[category] == value],
            color=color,
            title=value,
            ax=axis,
            wspace=0.1,
            legend_loc="on data",
            legend_fontsize=12,
            legend_fontweight="normal",
            legend_fontoutline=2,
            show=False,
        )

    save_fig(fig, out_dir, file_name=file_name, fig_size={"width": len(values) * 5, "height": 5})


def dim_reduction_cluster(
    _adata: sc.AnnData,
    method: str,
    n_pcs: int,
    n_neighbors: int,
    resolution: float,
    batch_effect: Optional[str],
    out_dir: PathLike,
) -> sc.AnnData:
    out_path = _ensure_dir(out_dir)
    clustering_method = _normalize_method_name(method)
    batch_method = _normalize_method_name(batch_effect)

    _adata.layers["counts"] = _adata.X.copy()
    sc.pp.normalize_total(_adata)
    sc.pp.log1p(_adata)
    _adata.raw = _adata

    sc.pp.highly_variable_genes(
        _adata,
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_top_genes=2000,
    )
    sc.pl.highly_variable_genes(_adata, show=False)
    save_fig(plt.gcf(), out_path, file_name="highly_variable_genes", fig_size={"width": 16, "height": 8})

    _adata = _adata[:, _adata.var["highly_variable"]].copy()
    sc.pp.regress_out(_adata, ["total_counts", "pct_counts_mt"])
    sc.pp.scale(_adata, max_value=10)

    sc.tl.pca(_adata, svd_solver="arpack")
    sc.pl.pca(_adata, color="sample", show=False)
    save_fig(plt.gcf(), out_path, file_name="pca_sample")
    sc.pl.pca(_adata, color="group", show=False)
    save_fig(plt.gcf(), out_path, file_name="pca_group")
    sc.pl.pca_variance_ratio(_adata, log=True, show=False)
    save_fig(plt.gcf(), out_path, file_name="pca_variance_ratio")

    if batch_method == "harmony":
        sce.pp.harmony_integrate(_adata, "sample", adjusted_basis="X_pca")

    if batch_method == "bbknn":
        sce.pp.bbknn(_adata, batch_key="sample", n_pcs=n_pcs)
    else:
        sc.pp.neighbors(_adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    sc.tl.umap(_adata)
    if clustering_method == "leiden":
        sc.tl.leiden(_adata, resolution=resolution, key_added="cluster")
    elif clustering_method == "louvain":
        sc.tl.louvain(_adata, resolution=resolution, key_added="cluster")
    else:
        raise ValueError("method must be either leiden or louvain.")

    sc.pl.umap(
        _adata,
        color="cluster",
        title="Cluster",
        legend_loc="on data",
        legend_fontsize=10,
        legend_fontoutline=2,
        show=False,
    )
    save_fig(plt.gcf(), out_path, file_name="umap_dimplot")
    sc.pl.umap(_adata, color="sample", title="", show=False)
    save_fig(plt.gcf(), out_path, file_name="umap_sample_dimplot")
    sc.pl.umap(_adata, color="group", title="", show=False)
    save_fig(plt.gcf(), out_path, file_name="umap_group_dimplot")

    with rc_context({"figure.figsize": (4, 4)}):
        sc.pl.umap(
            _adata,
            color=["cluster", "total_counts", "n_genes_by_counts", "total_counts_mt"],
            title=["Cluster", "total_counts", "n_genes_by_counts", "total_counts_mt"],
            legend_loc="on data",
            frameon=True,
            legend_fontsize=9,
            legend_fontoutline=3,
            wspace=0.2,
            ncols=2,
            cmap="coolwarm",
            use_raw=True,
            show=False,
        )
        plt.savefig(out_path / "cluster_total_counts_n_genes_by_counts.png", facecolor="white", bbox_inches="tight", dpi=300)
        plt.close(plt.gcf())

    _adata.obs.to_csv(out_path / "all_cells_features.xls", sep="\t")

    result = (
        _adata.obs.groupby(["cluster", "sample"])
        .size()
        .reset_index(name="number")
        .pivot(index="cluster", columns="sample", values="number")
        .reset_index()
    )
    result.columns.name = None
    result.to_csv(out_path / "cluster_cells_statistics.xls", sep="\t", index=None)

    with rc_context({"figure.figsize": (4, 4)}):
        sc.pl.correlation_matrix(_adata, "cluster", show=False)
        plt.savefig(out_path / "cluster_correlation.png", facecolor="white", bbox_inches="tight", dpi=300)
        plt.close(plt.gcf())

    _save_embedding_split_panels(_adata, "group", "cluster", out_path, "umap_split_group_cluster")
    _save_embedding_split_panels(_adata, "sample", "cluster", out_path, "umap_split_sample_cluster")

    return _adata


def find_marker_gene(
    _adata: sc.AnnData,
    method: str,
    out_dir: PathLike,
    logfc: float = 1,
    fdr: float = 0.05,
    pct: float = 0.25,
) -> sc.AnnData:
    out_path = _as_path(out_dir)
    if out_path.exists():
        shutil.rmtree(out_path)
    out_path.mkdir(parents=True, exist_ok=True)

    sc.tl.rank_genes_groups(_adata, "cluster", method=method, pts=True)
    sc.pl.rank_genes_groups(_adata, n_genes=25, sharey=False, show=False)
    save_fig(plt.gcf(), out_path, file_name="rank_genes_groups", fig_size={"width": 20, "height": 12})

    result = _adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names

    all_diff = []
    for cluster_name in groups:
        cluster_dir = out_path / "cluster{}".format(cluster_name)
        cluster_dir.mkdir(parents=True, exist_ok=True)

        tmp = sc.get.rank_genes_groups_df(_adata, group=cluster_name)
        tmp = tmp[
            (tmp["logfoldchanges"].abs() > logfc)
            & (tmp["pvals_adj"] < fdr)
            & (tmp["pct_nz_group"] >= pct)
        ].copy()
        tmp["cluster"] = cluster_name
        tmp.to_csv(cluster_dir / "differentially_expressed_genes.xls", sep="\t", index=None)

        if not tmp.empty:
            all_diff.append(tmp)
            top_genes = list(tmp["names"].head(5))
            if top_genes:
                sc.pl.rank_genes_groups_violin(
                    _adata,
                    groups=cluster_name,
                    gene_names=top_genes,
                    jitter=0,
                    size=0,
                    show=False,
                )
                save_fig(plt.gcf(), cluster_dir, file_name="top_5_genes_violin")

    if all_diff:
        merged = pd.concat(all_diff, axis=0)
        merged.rename(columns={"names": "Symbol"}, inplace=True)
        merged.to_csv(out_path / "all_cluster_marker.xls", sep="\t", index=None)
        merged["cluster"].value_counts().sort_index().to_csv(out_path / "cluster_diff_statistics.xls", sep="\t")
    else:
        pd.DataFrame(columns=["Symbol", "cluster"]).to_csv(
            out_path / "all_cluster_marker.xls",
            sep="\t",
            index=None,
        )
        pd.Series(dtype=int).to_csv(out_path / "cluster_diff_statistics.xls", sep="\t")

    sc.pl.rank_genes_groups_dotplot(
        _adata,
        n_genes=3,
        values_to_plot="logfoldchanges",
        vmin=-3,
        vmax=3,
        min_logfoldchange=1,
        cmap="bwr",
        colorbar_title="log fold change",
        show=False,
    )
    save_fig(plt.gcf(), out_path, file_name="dotplot_top3_marker_gene", fig_size={"width": 20, "height": 8})

    sc.pl.rank_genes_groups_heatmap(
        _adata,
        n_genes=10,
        values_to_plot="logfoldchanges",
        min_logfoldchange=1,
        swap_axes=True,
        cmap="bwr",
        standard_scale="var",
        show_gene_labels=True,
        dendrogram=True,
        show=False,
    )
    save_fig(plt.gcf(), out_path, file_name="heatmap_top10_marker_gene", fig_size={"width": 26, "height": 20})

    return _adata


def _coerce_celltype_mapping(
    cell_type: Optional[Mapping[str, str]],
    clusters: Iterable[str],
) -> Optional[Dict[str, str]]:
    if cell_type is None:
        return None

    mapping = {str(key): str(value) for key, value in cell_type.items()}
    missing_clusters = sorted(set(map(str, clusters)).difference(mapping.keys()))
    if missing_clusters:
        raise ValueError(
            "cell_type mapping is missing clusters: {}".format(", ".join(missing_clusters))
        )
    return mapping


def plot_celltype_proportion(
    _adata: sc.AnnData,
    out_dir: PathLike,
    type: str = "sample",
    fig_size: Optional[Mapping[str, float]] = None,
) -> None:
    if not pd.api.types.is_categorical_dtype(_adata.obs["cell_type"]):
        _adata.obs["cell_type"] = pd.Categorical(_adata.obs["cell_type"])

    keys = _adata.obs["cell_type"].cat.categories.to_list()
    values = list(_adata.uns["cell_type_colors"])
    cell_type_colors = dict(zip(keys, values))
    proportion_table = pd.crosstab(
        _adata.obs[type],
        _adata.obs["cell_type"],
        normalize="index",
    )

    ax = proportion_table.plot(
        kind="bar",
        xlabel="",
        stacked=True,
        color=[cell_type_colors[cell] for cell in proportion_table.columns],
        edgecolor="black",
        width=0.8,
    )
    ax.set_xlabel(type.capitalize(), fontsize=14, family="sans-serif")
    ax.set_ylabel("Proportion", fontsize=14, family="sans-serif")

    legend_font = fm.FontProperties(family="sans-serif", size=10)
    ax.legend(
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        borderaxespad=0.0,
        fontsize=10,
        title="Cell Types",
        title_fontsize=12,
        prop=legend_font,
        frameon=False,
        handletextpad=1.5,
        labelspacing=1.0,
        columnspacing=1.0,
        markerscale=1,
        alignment="left",
        handlelength=2,
        handleheight=2,
    )
    ax.tick_params(axis="x", labelsize=12, labelrotation=90, labelcolor="black")
    ax.tick_params(axis="y", labelsize=12, labelcolor="black")
    ax.grid(False)
    ax.set_facecolor("white")

    for spine in ax.spines.values():
        spine.set_edgecolor("black")
        spine.set_linewidth(1.5)

    if fig_size is None:
        fig_size = {"width": max(len(set(_adata.obs[type])) * 2, 6), "height": 6}

    save_fig(plt.gcf(), out_dir, file_name="celltype_proportion_{}".format(type), fig_size=fig_size)


def celltype_annotation(
    _adata: sc.AnnData,
    database: str,
    organism: str,
    out_dir: PathLike,
    cell_type: Optional[Mapping[str, str]] = None,
    celltypist_model: Optional[str] = None,
) -> Tuple[sc.AnnData, Dict[str, str]]:
    out_path = _as_path(out_dir)
    if out_path.exists():
        shutil.rmtree(out_path)
    out_path.mkdir(parents=True, exist_ok=True)

    database_name = _normalize_method_name(database)
    organism_name = _normalize_method_name(organism)
    cluster_mapping = _coerce_celltype_mapping(cell_type, _adata.obs["cluster"].astype(str))

    if organism_name == "mouse":
        _adata.var["gene"] = _adata.var.index
        _adata.var.index = [gene.upper() for gene in _adata.var.index]

    if database_name == "celltypist":
        import celltypist

        if cluster_mapping:
            annotation_dict = cluster_mapping
            _adata.obs["cell_type"] = [annotation_dict[str(cluster)] for cluster in _adata.obs["cluster"]]
        else:
            if celltypist_model:
                model_name = celltypist_model
            elif organism_name == "human":
                model_name = "Immune_All_Low.pkl"
            elif organism_name == "mouse":
                model_name = "Adult_Mouse_Gut.pkl"
            else:
                raise ValueError("organism must be either human or mouse.")

            predictions = celltypist.annotate(
                _adata,
                model=model_name,
                majority_voting=True,
            )
            _adata = predictions.to_adata()
            _adata.obs["cell_type"] = _adata.obs["majority_voting"]
            grouped = (
                _adata.obs[["cluster", "cell_type"]]
                .groupby("cluster")["cell_type"]
                .apply(list)
                .to_dict()
            )
            annotation_dict = {
                str(cluster): "|".join(sorted(set(map(str, labels))))
                for cluster, labels in grouped.items()
            }
    else:
        raise ValueError("database must be CellTypist.")

    _adata.obs["cell_type"] = pd.Categorical(_adata.obs["cell_type"])

    with rc_context({"figure.figsize": (5, 5)}):
        sc.pl.umap(
            _adata,
            color="cell_type",
            show=False,
            title="",
            legend_loc="on data",
            legend_fontsize=9,
            legend_fontweight="normal",
            legend_fontoutline=2,
        )
        plt.savefig(out_path / "umap_dimplot_celltype.png", facecolor="white", bbox_inches="tight", dpi=300)
        plt.close(plt.gcf())

    _save_embedding_split_panels(_adata, "group", "cell_type", out_path, "umap_split_group_celltype")

    result = (
        _adata.obs.groupby(["cell_type", "sample"])
        .size()
        .reset_index(name="number")
        .pivot(index="cell_type", columns="sample", values="number")
        .reset_index()
    )
    result.columns.name = None
    result.to_csv(out_path / "celltype_statistics.xls", sep="\t", index=None)

    _adata.obs.to_csv(out_path / "all_cells_features.xls", sep="\t")
    plot_celltype_proportion(_adata, out_dir=out_path, type="group")
    plot_celltype_proportion(_adata, out_dir=out_path, type="sample")
    pd.DataFrame(
        _adata.obsm["X_umap"],
        index=_adata.obs_names,
        columns=["UMAP_1", "UMAP_2"],
    ).to_csv(out_path / "data_umap.xls", sep="\t")

    with rc_context({"figure.figsize": (5, 5)}):
        sc.pl.umap(
            _adata,
            color=["cluster", "cell_type"],
            wspace=0.1,
            legend_loc="on data",
            legend_fontsize=12,
            legend_fontweight="normal",
            legend_fontoutline=2,
            show=False,
        )
        plt.savefig(out_path / "umap_cluster_celltype.png", facecolor="white", bbox_inches="tight", dpi=300)
        plt.close(plt.gcf())

    if organism_name == "mouse":
        _adata.var.index = _adata.var["gene"]

    _adata.write(out_path / "result.h5ad")
    return _adata, annotation_dict


def visua_gene_express(
    _adata: sc.AnnData,
    genes: Sequence[str],
    out_dir: PathLike,
) -> Tuple[Image.Image, plt.Figure]:
    out_path = _ensure_dir(out_dir)

    with rc_context({"figure.figsize": (4, 4)}):
        sc.pl.umap(_adata, color=list(genes), s=15, ncols=3, cmap="coolwarm", use_raw=True, show=False)
        plt.savefig(out_path / "{}.png".format("_".join(genes)), facecolor="white", bbox_inches="tight", dpi=200)
        plt.close(plt.gcf())

    gene_expression = Image.open(out_path / "{}.png".format("_".join(genes)))

    if len(genes) > 1:
        fig, axes = plt.subplots(len(genes), 1, figsize=(10, len(genes) * 4), sharey=True)
        for index, gene_name in enumerate(genes):
            sc.pl.violin(
                _adata,
                gene_name,
                groupby="cluster",
                stripplot=False,
                xlabel="",
                inner="box",
                ax=axes[index],
                show=False,
            )
    else:
        fig, axis = plt.subplots(figsize=(4, 4))
        sc.pl.violin(
            _adata,
            genes[0],
            groupby="cluster",
            stripplot=False,
            xlabel="",
            inner="box",
            ax=axis,
            show=False,
        )

    return gene_expression, fig
