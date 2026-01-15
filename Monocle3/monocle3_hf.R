## ===============================
## 04_Monocle3 trajectory (CD8)
## ===============================
library(monocle3)
library(ggplot2)
library(dplyr)
library(Seurat)

## -------------------------------
## paths
## -------------------------------
mono_dir <- file.path(out_dir, "03_CD8_Subcluster/05_CD8_monocle3")
mono_rds <- file.path(mono_dir, "rds")
mono_fig <- file.path(mono_dir, "figures")

dir.create(mono_rds, recursive = TRUE, showWarnings = FALSE)
dir.create(mono_fig, recursive = TRUE, showWarnings = FALSE)

## -------------------------------
## load input object
## （明确来源：CD8 clean 对象）
## -------------------------------
scRNA <- readRDS(
  file.path(out_dir, "rds/cd8_clean_final.rds")
)

## -------------------------------
## Seurat -> Monocle3 CDS
## -------------------------------
rds2cds <- function(seuset, umap2dfile = NULL) {
  
  data <- GetAssayData(seuset[["RNA"]], layer = "counts")
  
  gene_metadata <- data.frame(
    gene_short_name = rownames(data),
    row.names = rownames(data),
    check.names = FALSE
  )
  
  cell_metadata <- seuset@meta.data
  cell_metadata$Clusters <- as.character(seuset@active.ident)
  
  cds <- new_cell_data_set(
    data,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )
  
  cds <- estimate_size_factors(cds)
  
  ## --- add UMAP coordinates ---
  if (is.null(umap2dfile)) {
    umap_df <- as.data.frame(seuset@reductions$umap@cell.embeddings)
  } else {
    umap_df <- read.table(umap2dfile, sep = "\t", header = TRUE, row.names = 1)
    umap_df <- umap_df[colnames(cds), , drop = FALSE]
  }
  
  reducedDims(cds)$UMAP <- umap_df
  
  ## --- generate monocle internal cluster structure ---
  cds <- cluster_cells(
    cds,
    reduction_method = "UMAP",
    cluster_method = "louvain",
    k = 30
  )
  
  ## --- overwrite clusters with Seurat clusters ---
  cds@colData$Clusters <- factor(cds@colData$Clusters)
  cds@clusters$UMAP$clusters <- cds@colData$Clusters
  
  cds
}

cds <- rds2cds(scRNA)
saveRDS(cds, file.path(mono_rds, "cd8_cds_raw.rds"))

## -------------------------------
## learn graph
## -------------------------------
cds <- learn_graph(
  cds,
  use_partition = TRUE,
  close_loop = FALSE,
  learn_graph_control = list(
    minimal_branch_len = 25,
    nn.k = 30
  ),
  verbose = TRUE
)

saveRDS(cds, file.path(mono_rds, "cd8_cds_graph.rds"))

## -------------------------------
## trajectory plots
## -------------------------------
p_cluster <- plot_cells(
  cds,
  reduction_method = "UMAP",
  color_cells_by = "Clusters",
  show_trajectory_graph = FALSE,
  label_cell_groups = TRUE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  cell_size = 1
)

ggsave(
  file.path(mono_fig, "Trajectory_cluster.png"),
  p_cluster, width = 15, height = 15
)
ggsave(
  file.path(mono_fig, "Trajectory_cluster.pdf"),
  p_cluster, width = 15, height = 15
)

p_branch <- plot_cells(
  cds,
  reduction_method = "UMAP",
  color_cells_by = "cluster",
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  trajectory_graph_segment_size = 2,
  graph_label_size = 2.5,
  cell_size = 1
)

ggsave(
  file.path(mono_fig, "Trajectory_branch.png"),
  p_branch, width = 15, height = 15
)
ggsave(
  file.path(mono_fig, "Trajectory_branch.pdf"),
  p_branch, width = 15, height = 15
)

## -------------------------------
## order cells (root = CD8 naive)
## -------------------------------
get_earliest_principal_node <- function(cds, time_bin = "Tpro CD8") {
  cell_ids <- which(colData(cds)[, "cd8_cell_type"] == time_bin)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
  ]
}

cds <- order_cells(
  cds,
  root_pr_nodes = get_earliest_principal_node(cds)
)

saveRDS(cds, file.path(mono_rds, "cd8_cds_ordered.rds"))

## -------------------------------
## pseudotime plot
## -------------------------------
p_pseudo <- plot_cells(
  cds,
  reduction_method = "UMAP",
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  trajectory_graph_segment_size = 1.5,
  graph_label_size = 3,
  cell_size = 1
)

ggsave(
  file.path(mono_fig, "Trajectory_pseudotime.png"),
  p_pseudo, width = 15, height = 15
)
ggsave(
  file.path(mono_fig, "Trajectory_pseudotime.pdf"),
  p_pseudo, width = 15, height = 15
)