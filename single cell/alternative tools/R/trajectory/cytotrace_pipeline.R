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
  library(Seurat)
  library(CytoTRACE)
  library(ggplot2)
  library(cowplot)
})

load_mapping_table <- function(path_value) {
  if (is.null(path_value) || path_value == "") {
    return(NULL)
  }

  read.table(
    path_value,
    sep = "\t",
    header = FALSE,
    comment.char = "",
    stringsAsFactors = FALSE
  )
}

metadata_plot <- function(object, feature_name, reduction_name = "umap", point_size = 1) {
  dim_names <- paste0(Key(object = object[[reduction_name]]), c(1, 2))
  plot_df <- FetchData(object = object, vars = c(dim_names, feature_name))

  ggplot(plot_df, aes_string(x = dim_names[1], y = dim_names[2], color = feature_name)) +
    geom_point(size = point_size) +
    scale_color_gradientn(colours = rainbow(8)[6:1]) +
    theme_cowplot() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.title = element_text(size = 12)
    ) +
    labs(color = "Predicted\norder", title = feature_name)
}

prepare_seurat_object <- function(seu, ident_col = NULL, celllist_file = NULL) {
  if (!is.null(celllist_file)) {
    cell_map <- load_mapping_table(celllist_file)
    if (ncol(cell_map) < 2) {
      stop("celllist file must contain at least two columns: cell and cluster.")
    }
    cell_map <- cell_map[cell_map[[1]] %in% colnames(seu), , drop = FALSE]
    seu <- subset(seu, cells = cell_map[[1]])
    Idents(seu) <- stats::setNames(cell_map[[2]], cell_map[[1]])[colnames(seu)]
    return(seu)
  }

  if (!is.null(ident_col) && ident_col %in% colnames(seu@meta.data)) {
    Idents(seu) <- seu[[ident_col, drop = TRUE]]
  }

  seu
}

random_sample_clusters <- function(cluster_vector, cells_percent, min_cells_per_cluster) {
  samples <- c()
  for (cluster_name in unique(cluster_vector)) {
    current_cells <- names(cluster_vector)[cluster_vector == cluster_name]
    target_size <- max(round(length(current_cells) * cells_percent), min_cells_per_cluster)
    target_size <- min(target_size, length(current_cells))
    samples <- c(samples, sample(current_cells, target_size, replace = FALSE))
  }
  unique(samples)
}

build_boxplot <- function(plot_df, group_name, color_map = NULL) {
  median_values <- tapply(plot_df$CytoTRACE, plot_df$Cluster, stats::median)
  order_levels <- names(sort(median_values, decreasing = TRUE))
  plot_df$Cluster <- factor(plot_df$Cluster, levels = order_levels)

  p <- ggplot(plot_df, aes(x = Cluster, y = CytoTRACE, color = Cluster, fill = Cluster)) +
    geom_boxplot(alpha = 0.5, outlier.size = 0) +
    geom_jitter(height = 0, size = 0.8, width = 0.1) +
    theme_cowplot() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(x = group_name, y = "CytoTRACE")

  if (!is.null(color_map)) {
    use_colors <- color_map[order_levels]
    p <- p +
      scale_color_manual(values = use_colors) +
      scale_fill_manual(values = use_colors)
  }

  p
}

args <- parse_cli_args()
input_rds <- cli_string(args, "input-rds")
out_dir <- ensure_dir(cli_string(args, "out-dir", "cytotrace_result"))
plot_gene_num <- cli_numeric(args, "plot-gene-num", 10)
color_tab <- cli_string(args, "color-tab")
celllist_file <- cli_string(args, "celllist")
ident_col <- cli_string(args, "ident-col")
group_col <- cli_string(args, "group-col", "group")
random_cells <- cli_bool(args, "random-cells", FALSE)
cells_percent <- cli_numeric(args, "cells-percent", 0.1)
min_cells_per_cluster <- cli_numeric(args, "min-cells-per-cluster", 1000)

if (is.null(input_rds)) {
  stop("Please provide --input-rds.")
}

seu <- readRDS(input_rds)
seu <- prepare_seurat_object(seu, ident_col = ident_col, celllist_file = celllist_file)

if (random_cells) {
  sampled_cells <- random_sample_clusters(Idents(seu), cells_percent, min_cells_per_cluster)
  seu <- subset(seu, cells = sampled_cells)
  write.table(
    data.frame(Cell = colnames(seu), Cluster = Idents(seu)),
    file.path(out_dir, "random_cells.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

if ("SCT" %in% names(seu@assays)) {
  DefaultAssay(seu) <- "SCT"
} else {
  DefaultAssay(seu) <- "RNA"
}

counts_mat <- GetAssayData(seu, slot = "counts")
counts_mat <- as.matrix(counts_mat[, colnames(seu), drop = FALSE])
counts_mat <- stats::na.omit(counts_mat)

cyto_results <- CytoTRACE(counts_mat, ncores = 1, enableFast = TRUE)
saveRDS(cyto_results, file.path(out_dir, "cytotrace_raw_result.rds"))

cyto_value <- cyto_results$CytoTRACE
cyto_meta <- data.frame(CytoTRACE = cyto_value, row.names = names(cyto_value), check.names = FALSE)
seu <- AddMetaData(seu, metadata = cyto_meta)

write.table(
  data.frame(Cell = names(cyto_value), CytoTRACE = cyto_value),
  file.path(out_dir, "cytotrace_value.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cyto_gene <- sort(cyto_results$cytoGenes, decreasing = TRUE)
gene_df <- data.frame(Gene = names(cyto_gene), GenesCorrelated = cyto_gene)
write.table(
  gene_df,
  file.path(out_dir, "cytotrace_gene.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if ("umap" %in% names(seu@reductions)) {
  p_umap <- metadata_plot(seu, "CytoTRACE", reduction_name = "umap")
  save_plot(p_umap, file.path(out_dir, "cytotrace_umap"), width = 9, height = 8)
}

if ("tsne" %in% names(seu@reductions)) {
  p_tsne <- metadata_plot(seu, "CytoTRACE", reduction_name = "tsne")
  save_plot(p_tsne, file.path(out_dir, "cytotrace_tsne"), width = 9, height = 8)
}

color_map <- NULL
color_table <- load_mapping_table(color_tab)
if (!is.null(color_table) && ncol(color_table) >= 2) {
  color_map <- stats::setNames(color_table[[2]], color_table[[1]])
}

cluster_df <- data.frame(
  Cluster = as.character(Idents(seu)),
  CytoTRACE = seu$CytoTRACE,
  stringsAsFactors = FALSE
)
write.table(
  cluster_df,
  file.path(out_dir, "cytotrace_by_cluster.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_cluster <- build_boxplot(cluster_df, group_name = ident_col %||% "active.ident", color_map = color_map)
save_plot(
  p_cluster,
  file.path(out_dir, "cytotrace_boxplot_cluster"),
  width = max(7, 0.4 * length(unique(cluster_df$Cluster)) + 6),
  height = 8
)

if (!is.null(group_col) && group_col %in% colnames(seu@meta.data)) {
  group_df <- data.frame(
    Cluster = as.character(seu@meta.data[[group_col]]),
    CytoTRACE = seu$CytoTRACE,
    stringsAsFactors = FALSE
  )
  p_group <- build_boxplot(group_df, group_name = group_col)
  save_plot(
    p_group,
    file.path(out_dir, "cytotrace_boxplot_group"),
    width = max(7, 0.4 * length(unique(group_df$Cluster)) + 6),
    height = 8
  )
}

gene_bar_df <- rbind(head(gene_df, plot_gene_num), tail(gene_df, plot_gene_num))
gene_bar_df$Negative <- gene_bar_df$GenesCorrelated <= 0
gene_bar_df$Gene <- factor(gene_bar_df$Gene, levels = rev(gene_bar_df$Gene))
limits_value <- ceiling(max(abs(gene_bar_df$GenesCorrelated)) * 100) / 100

p_gene <- ggplot(gene_bar_df, aes(x = Gene, y = GenesCorrelated, fill = Negative)) +
  geom_col() +
  coord_flip() +
  ylim(-limits_value, limits_value) +
  theme_cowplot() +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Correlation with CytoTRACE")

save_plot(
  p_gene,
  file.path(out_dir, "cytotrace_gene_barplot"),
  width = 7,
  height = max(5, 0.6 * plot_gene_num)
)

saveRDS(seu, file.path(out_dir, "cytotrace_seurat.rds"))
message("CytoTRACE finished: ", out_dir)
