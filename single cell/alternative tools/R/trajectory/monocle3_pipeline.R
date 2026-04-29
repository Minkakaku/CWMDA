rm(list = ls())

get_script_dir_bootstrap <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)

  if (length(file_arg) > 0) {
    return(dirname(normalizePath(
      gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE),
      winslash = "/",
      mustWork = FALSE
    )))
  }

  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(
      sys.frames()[[1]]$ofile,
      winslash = "/",
      mustWork = FALSE
    )))
  }

  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_dir <- get_script_dir_bootstrap()
source(file.path(script_dir, "common_utils.R"), local = FALSE)

suppressPackageStartupMessages({
  library(monocle3)
  library(Seurat)
  library(ggplot2)
})

build_cds_from_seurat <- function(seu, assay_name, reduction_name, cluster_col = NULL) {
  counts_mat <- tryCatch(
    GetAssayData(seu, assay = assay_name, layer = "counts"),
    error = function(e) GetAssayData(seu, assay = assay_name, slot = "counts")
  )

  gene_metadata <- data.frame(
    gene_short_name = rownames(counts_mat),
    row.names = rownames(counts_mat),
    check.names = FALSE
  )

  cell_metadata <- seu@meta.data
  if (!is.null(cluster_col) && cluster_col %in% colnames(cell_metadata)) {
    cell_metadata$Clusters <- as.character(cell_metadata[[cluster_col]])
  } else {
    cell_metadata$Clusters <- as.character(Idents(seu))
  }

  cds <- new_cell_data_set(
    expression_data = counts_mat,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )
  cds <- estimate_size_factors(cds)

  reduction_label <- toupper(reduction_name)
  reducedDims(cds)[[reduction_label]] <- seu@reductions[[reduction_name]]@cell.embeddings

  cds <- cluster_cells(
    cds,
    reduction_method = reduction_label,
    cluster_method = "louvain",
    k = 30
  )

  names(cds@colData$Clusters) <- rownames(cds@colData)
  cds@colData$Clusters <- factor(cds@colData$Clusters)
  cds@clusters[[reduction_label]]$clusters <- cds@colData$Clusters
  cds
}

get_root_pr_node <- function(cds, reduction_name, root_col, root_value) {
  cell_ids <- which(colData(cds)[, root_col] == root_value)
  if (length(cell_ids) == 0) {
    stop("No cells found for root value ", root_value, " in column ", root_col)
  }

  reduction_label <- toupper(reduction_name)
  closest_vertex <- cds@principal_graph_aux[[reduction_label]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), , drop = FALSE])

  igraph::V(principal_graph(cds)[[reduction_label]])$name[
    as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
  ]
}

args <- parse_cli_args()
input_rds <- cli_string(args, "input-rds")
out_dir <- ensure_dir(cli_string(args, "out-dir", "monocle3_result"))
assay_name <- cli_string(args, "assay", "RNA")
reduction_name <- cli_string(args, "reduction", "umap")
cluster_col <- cli_string(args, "cluster-col")
root_col <- cli_string(args, "root-col")
root_value <- cli_string(args, "root-value")
minimal_branch_len <- cli_numeric(args, "minimal-branch-len", 25)
nn_k <- cli_numeric(args, "nn-k", 30)
use_partition <- cli_bool(args, "use-partition", TRUE)
cell_size <- cli_numeric(args, "cell-size", 1)

if (is.null(input_rds)) {
  stop("Please provide --input-rds.")
}

seu <- readRDS(input_rds)
if (!reduction_name %in% names(seu@reductions)) {
  stop("reduction not found: ", reduction_name)
}

cds <- build_cds_from_seurat(
  seu = seu,
  assay_name = assay_name,
  reduction_name = reduction_name,
  cluster_col = cluster_col
)
saveRDS(cds, file.path(out_dir, "monocle3_raw_cds.rds"))

reduction_label <- toupper(reduction_name)
cds <- learn_graph(
  cds,
  use_partition = use_partition,
  close_loop = FALSE,
  learn_graph_control = list(
    minimal_branch_len = minimal_branch_len,
    nn.k = nn_k
  ),
  verbose = TRUE
)
saveRDS(cds, file.path(out_dir, "monocle3_graph_cds.rds"))

p_cluster <- plot_cells(
  cds,
  reduction_method = reduction_label,
  color_cells_by = "Clusters",
  show_trajectory_graph = FALSE,
  label_cell_groups = TRUE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  cell_size = cell_size
)
save_plot(p_cluster, file.path(out_dir, "trajectory_cluster"), width = 14, height = 12)

p_branch <- plot_cells(
  cds,
  reduction_method = reduction_label,
  color_cells_by = "cluster",
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  trajectory_graph_segment_size = 2,
  graph_label_size = 2.5,
  cell_size = cell_size
)
save_plot(p_branch, file.path(out_dir, "trajectory_branch"), width = 14, height = 12)

if (!is.null(root_col) && !is.null(root_value)) {
  if (!root_col %in% colnames(colData(cds))) {
    stop("root column not found: ", root_col)
  }

  cds <- order_cells(
    cds,
    root_pr_nodes = get_root_pr_node(cds, reduction_name, root_col, root_value)
  )
  saveRDS(cds, file.path(out_dir, "monocle3_ordered_cds.rds"))

  p_pseudotime <- plot_cells(
    cds,
    reduction_method = reduction_label,
    color_cells_by = "pseudotime",
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    trajectory_graph_segment_size = 1.5,
    graph_label_size = 3,
    cell_size = cell_size
  )
  save_plot(p_pseudotime, file.path(out_dir, "trajectory_pseudotime"), width = 14, height = 12)
}

message("Monocle3 finished: ", out_dir)
