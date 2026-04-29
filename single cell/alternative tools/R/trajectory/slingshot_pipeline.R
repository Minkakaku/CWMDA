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
  library(SingleCellExperiment)
  library(slingshot)
  library(ggplot2)
  library(RColorBrewer)
})

build_discrete_palette <- function(levels_value, palette_name) {
  n_value <- length(levels_value)
  color_values <- if (n_value <= 8) {
    brewer.pal(max(3, n_value), palette_name)[seq_len(n_value)]
  } else {
    colorRampPalette(brewer.pal(8, palette_name))(n_value)
  }
  stats::setNames(color_values, levels_value)
}

args <- parse_cli_args()
input_rds <- cli_string(args, "input-rds")
out_dir <- ensure_dir(cli_string(args, "out-dir", "slingshot_result"))
cluster_col <- cli_string(args, "cluster-col")
group_col <- cli_string(args, "group-col")
reduction_name <- cli_string(args, "reduction", "umap")
start_cluster <- cli_string(args, "start-cluster")

if (is.null(input_rds)) {
  stop("Please provide --input-rds.")
}

seu <- readRDS(input_rds)

if (!reduction_name %in% names(seu@reductions)) {
  stop("reduction not found: ", reduction_name)
}

cluster_values <- if (!is.null(cluster_col) && cluster_col %in% colnames(seu@meta.data)) {
  as.character(seu@meta.data[[cluster_col]])
} else {
  as.character(Idents(seu))
}

sce <- as.SingleCellExperiment(seu)
reducedDims(sce)$CUSTOM <- seu@reductions[[reduction_name]]@cell.embeddings
sce$cluster <- factor(cluster_values)

slingshot_args <- list(
  data = sce,
  clusterLabels = "cluster",
  reducedDim = "CUSTOM"
)
if (!is.null(start_cluster) && start_cluster != "") {
  sce <- slingshot::slingshot(
    data = sce,
    clusterLabels = "cluster",
    reducedDim = "CUSTOM",
    start.clus = start_cluster
  )
} else {
  sce <- slingshot::slingshot(
    data = sce,
    clusterLabels = "cluster",
    reducedDim = "CUSTOM"
  )
}

pt <- slingPseudotime(sce)
sce$pseudotime_slingshot <- pt[, 1]
seu$pseudotime_slingshot <- sce$pseudotime_slingshot

saveRDS(sce, file.path(out_dir, "sce_slingshot.rds"))
saveRDS(seu, file.path(out_dir, "seurat_with_slingshot_pseudotime.rds"))

embedding_df <- as.data.frame(seu@reductions[[reduction_name]]@cell.embeddings)
colnames(embedding_df)[1:2] <- c("Dim1", "Dim2")
embedding_df$cluster <- cluster_values
embedding_df$pseudotime <- seu$pseudotime_slingshot

if (!is.null(group_col) && group_col %in% colnames(seu@meta.data)) {
  embedding_df$group <- as.character(seu@meta.data[[group_col]])
}

curve_list <- slingCurves(sce)
curve_rows <- lapply(seq_along(curve_list), function(index_value) {
    df <- as.data.frame(curve_list[[index_value]]$s)
    colnames(df) <- c("Dim1", "Dim2")
    df$lineage <- paste0("Lineage_", index_value)
    df
})
curve_df <- if (length(curve_rows) == 1L) curve_rows[[1L]] else Reduce(rbind, curve_rows)

celltype_levels <- unique(embedding_df$cluster)
cluster_palette <- build_discrete_palette(celltype_levels, "Set2")

base_cluster_plot <- ggplot() +
  geom_point(
    data = embedding_df,
    aes(x = Dim1, y = Dim2, fill = cluster),
    shape = 21,
    size = 3,
    stroke = 0.2,
    color = "black"
  ) +
  geom_path(
    data = curve_df,
    aes(x = Dim1, y = Dim2, group = lineage),
    color = "black",
    linewidth = 1
  ) +
  scale_fill_manual(values = cluster_palette, name = "Cluster") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(title = "Slingshot trajectory by cluster")

base_pseudotime_plot <- ggplot() +
  geom_point(
    data = embedding_df,
    aes(x = Dim1, y = Dim2, fill = pseudotime),
    shape = 21,
    size = 3,
    stroke = 0.2,
    color = "black"
  ) +
  geom_path(
    data = curve_df,
    aes(x = Dim1, y = Dim2, group = lineage),
    color = "black",
    linewidth = 1
  ) +
  scale_fill_viridis_c(option = "plasma", na.value = "grey80", name = "Pseudotime") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(title = "Slingshot pseudotime")

if ("group" %in% colnames(embedding_df)) {
  p_cluster <- base_cluster_plot + facet_wrap(~group, nrow = 1)
  p_pseudotime <- base_pseudotime_plot + facet_wrap(~group, nrow = 1)
} else {
  p_cluster <- base_cluster_plot
  p_pseudotime <- base_pseudotime_plot
}

save_plot(p_cluster, file.path(out_dir, "slingshot_cluster"), width = 16, height = 7)
save_plot(p_pseudotime, file.path(out_dir, "slingshot_pseudotime"), width = 16, height = 7)

message("Slingshot finished: ", out_dir)
