if (!exists("script_dir")) {
  stop("请先在步骤脚本中定义 script_dir，再 source 这个公共脚本。")
}

project_dir <- normalizePath(
  file.path(script_dir, ".."),
  winslash = "/",
  mustWork = FALSE
)
setwd(project_dir)

required_packages <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "Startrac",
  "pheatmap",
  "patchwork",
  "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "缺少以下 R 包，请先安装后再运行：",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  invisible(lapply(required_packages, library, character.only = TRUE))
})

PARAM <- list(
  assay = "RNA",
  celltype_col = "cell_type",
  target_celltype = "Epithelial", # 改这里，选择要做亚群细分的主群
  sample_col = "orig.ident",
  use_harmony = FALSE,
  harmony_group = "orig.ident",
  regress_vars = c("nCount_RNA", "percent.mt"),
  nfeatures = 2000,
  npcs = 30,
  dims = 1:30,
  resolution = 0.6,
  remove_clusters = character(0),
  group_map = NULL,
  input = "04_rds/01.sc_annoted.rds"
)

SUBCLUSTER_NAME <- PARAM$target_celltype
SUB_ROOT <- "03_Subcluster"
SUB_DIR <- file.path(SUB_ROOT, SUBCLUSTER_NAME)
dir.create(SUB_DIR, recursive = TRUE, showWarnings = FALSE)

SEEDS <- list(
  subcluster = 2001L,
  markers = 2002L,
  annotation = 2003L,
  proportion = 2004L
)

RDS_CLUSTER <- file.path(SUB_DIR, "subcluster_clustered.rds")
RDS_CLUSTER_HARMONY <- file.path(SUB_DIR, "subcluster_harmony.rds")
RDS_CLUSTER_NO_HARMONY <- file.path(SUB_DIR, "subcluster_no_harmony.rds")
RDS_MARKERS <- file.path(SUB_DIR, "markers_all_clusters.rds")
RDS_ANNOTATED <- file.path(SUB_DIR, "subcluster_annotated.rds")
RDS_ROE <- file.path(SUB_DIR, "ROE_matrix.rds")

set_pipeline_seed <- function(seed) {
  set.seed(seed)
  invisible(seed)
}

save_plot <- function(p, filename, width = 8, height = 6) {
  ggplot2::ggsave(
    paste0(filename, ".png"),
    p,
    width = width,
    height = height,
    dpi = 300
  )
  ggplot2::ggsave(
    paste0(filename, ".pdf"),
    p,
    width = width,
    height = height
  )
}

save_pheatmap <- function(p, filename, width = 8, height = 6, dpi = 300) {
  grDevices::png(
    filename = paste0(filename, ".png"),
    width = width * dpi,
    height = height * dpi,
    res = dpi
  )
  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  grDevices::dev.off()

  grDevices::pdf(
    file = paste0(filename, ".pdf"),
    width = width,
    height = height
  )
  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  grDevices::dev.off()
}

plot_celltype_distribution <- function(
  data,
  sample_col = "sample",
  celltype_col = "celltype",
  palette_name = "Set1",
  output_dir = ".",
  output_file = "plotC.pdf",
  width = 7,
  height = 5,
  text_size = 12,
  x_angle = 45,
  legend_position = "right"
) {
  if (!inherits(data, "Seurat")) {
    stop("data 必须是 Seurat 对象")
  }
  if (!sample_col %in% colnames(data@meta.data)) {
    stop(paste0("meta.data 中未找到样本列：", sample_col))
  }
  if (!celltype_col %in% colnames(data@meta.data)) {
    stop(paste0("meta.data 中未找到细胞类型列：", celltype_col))
  }

  meta_data <- data@meta.data[, c(sample_col, celltype_col), drop = FALSE]
  colnames(meta_data) <- c("Sample", "CellType")

  if (anyNA(meta_data$CellType)) {
    stop("注释列中仍有 NA，请先完成亚群注释。")
  }

  meta_data$CellType <- factor(meta_data$CellType)
  colour_count <- length(levels(meta_data$CellType))

  if (!palette_name %in% rownames(RColorBrewer::brewer.pal.info)) {
    warning(paste0("颜色板 ", palette_name, " 不存在，已切换为 Set1"))
    palette_name <- "Set1"
  }

  max_colors <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
  base_colors <- RColorBrewer::brewer.pal(min(max_colors, colour_count), palette_name)
  celltype_colors <- grDevices::colorRampPalette(base_colors)(colour_count)
  names(celltype_colors) <- levels(meta_data$CellType)

  plot_data <- as.data.frame(table(meta_data$Sample, meta_data$CellType))
  colnames(plot_data) <- c("Sample", "CellType", "Number")

  p_count <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = Sample, y = Number, fill = CellType)
  ) +
    ggplot2::geom_bar(stat = "identity", width = 0.8) +
    ggplot2::scale_fill_manual(values = celltype_colors) +
    ggplot2::labs(x = NULL, y = "Cell number") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = text_size, colour = "black"),
      axis.title.y = ggplot2::element_text(size = text_size),
      axis.text.x = ggplot2::element_text(angle = x_angle, hjust = 0.5, vjust = 0.8),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA),
      legend.position = legend_position,
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = text_size - 1)
    )

  p_prop <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = Sample, y = Number, fill = CellType)
  ) +
    ggplot2::geom_bar(stat = "identity", width = 0.8, position = "fill") +
    ggplot2::scale_fill_manual(values = celltype_colors) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(x = NULL, y = "Cell proportion") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = text_size, colour = "black"),
      axis.title.y = ggplot2::element_text(size = text_size),
      axis.text.x = ggplot2::element_text(angle = x_angle, hjust = 0.5, vjust = 0.8),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA),
      legend.position = legend_position,
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = text_size - 1)
    )

  p_combined <- p_count + p_prop +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = legend_position)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  ggplot2::ggsave(
    filename = file.path(output_dir, output_file),
    plot = p_combined,
    width = width,
    height = height,
    device = "pdf",
    dpi = 600
  )

  invisible(p_combined)
}

join_layers_if_needed <- function(obj) {
  if (!"Seurat" %in% loadedNamespaces()) {
    return(obj)
  }
  if (!"JoinLayers" %in% getNamespaceExports("Seurat")) {
    return(obj)
  }

  tryCatch(
    Seurat::JoinLayers(obj),
    error = function(e) {
      message("JoinLayers 跳过：", conditionMessage(e))
      obj
    }
  )
}

prepare_subcluster_input <- function(obj) {
  if (!PARAM$sample_col %in% colnames(obj@meta.data)) {
    stop("输入对象缺少样本列：", PARAM$sample_col)
  }

  if (!"samples" %in% colnames(obj@meta.data)) {
    obj$samples <- obj[[PARAM$sample_col, drop = TRUE]]
  }

  if (!is.null(PARAM$group_map)) {
    obj$group <- unname(PARAM$group_map[obj[[PARAM$sample_col, drop = TRUE]]])
  } else if (!"group" %in% colnames(obj@meta.data)) {
    sample_values <- obj[[PARAM$sample_col, drop = TRUE]]
    obj$group <- dplyr::case_when(
      grepl("^Con", sample_values) ~ "CON",
      grepl("^Exp", sample_values) ~ "Exp",
      TRUE ~ NA_character_
    )
  }

  obj
}

read_required_rds <- function(path) {
  if (!file.exists(path)) {
    stop("未找到 RDS：", path)
  }
  readRDS(path)
}
