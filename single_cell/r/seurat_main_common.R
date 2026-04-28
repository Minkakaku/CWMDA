script_dir = "/home/hf/DataFile/split_pipeline"
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
  "SeuratWrappers",
  "SingleCellExperiment",
  "scDblFinder",
  "celda",
  "dplyr",
  "ggplot2",
  "patchwork",
  "scales",
  "harmony",
  "pheatmap"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "缺少以下 R 包，请先在服务器安装后再运行：",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  invisible(lapply(required_packages, library, character.only = TRUE))
})

PIPELINE_SEEDS <- list(
  qc = 1001L,
  clustering = 1002L,
  annotation = 1003L,
  proportion = 1004L
)

DATA_DIR <- "data2analysis"
QC_DIR <- "01_QC"
GLOBAL_DIR <- "02_GlobalClustering"
MARKER_DIR <- "03_FindMarkers"
ANNOTATION_DIR <- "04_CellAnnotation"
PROPORTION_DIR <- "05_CellProportion"
RDS_DIR <- "04_rds"

RDS_QC <- file.path(RDS_DIR, "01.sc_after_QC_merged.rds")
RDS_CLUSTER <- file.path(RDS_DIR, "02.sc_global_clustered.rds")
RDS_MARKERS <- file.path(RDS_DIR, "03.all_markers.rds")
RDS_ANNOTATED <- file.path(RDS_DIR, "04.sc_annotated.rds")
RDS_ANNOTATED_LEGACY <- file.path(RDS_DIR, "01.sc_annoted.rds")
RDS_FINAL <- file.path(RDS_DIR, "05.sc_final.rds")
RDS_ROE <- file.path(RDS_DIR, "04.roe_matrix.rds")

init_output_dirs <- function() {
  dirs <- c(QC_DIR, GLOBAL_DIR, MARKER_DIR, ANNOTATION_DIR, PROPORTION_DIR, RDS_DIR)
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

set_pipeline_seed <- function(seed) {
  set.seed(seed)
  invisible(seed)
}

save_plot <- function(p, filename, width = 8, height = 6) {
  ggplot2::ggsave(
    filename = paste0(filename, ".png"),
    plot = p,
    width = width,
    height = height,
    dpi = 300
  )
  ggplot2::ggsave(
    filename = paste0(filename, ".pdf"),
    plot = p,
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

make_qc_plots <- function(obj) {
  list(
    violin = Seurat::VlnPlot(
      obj,
      features = c(
        "nFeature_RNA",
        "nCount_RNA",
        "percent.mt",
        "percent.hb",
        "log10GenesPerUMI"
      ),
      ncol = 5,
      pt.size = 0.1
    ),
    scatter_gene_umi = Seurat::FeatureScatter(
      obj,
      feature1 = "nCount_RNA",
      feature2 = "nFeature_RNA"
    ),
    scatter_mt = Seurat::FeatureScatter(
      obj,
      feature1 = "nCount_RNA",
      feature2 = "percent.mt"
    ),
    cells_per_sample = ggplot2::ggplot(obj@meta.data, ggplot2::aes(orig.ident)) +
      ggplot2::geom_bar() +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      ),
    doublet_ratio = ggplot2::ggplot(
      obj@meta.data,
      ggplot2::aes(orig.ident, fill = doublet_class)
    ) +
      ggplot2::geom_bar(position = "fill") +
      ggplot2::scale_y_continuous(labels = scales::percent) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      ),
    contamination = Seurat::VlnPlot(
      obj,
      features = "Contamination",
      group.by = "orig.ident",
      pt.size = 0
    ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
  )
}

save_qc <- function(obj, stage) {
  plots <- make_qc_plots(obj)

  for (nm in names(plots)) {
    save_plot(
      p = plots[[nm]],
      filename = file.path(QC_DIR, paste0("QC_", nm, "_", stage)),
      width = ifelse(nm == "violin", 16, 10),
      height = 5
    )
  }
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
    stop("cell_type 中仍有 NA，请先完成所有 cluster 的注释。")
  }

  meta_data$CellType <- factor(meta_data$CellType)

  colour_count <- length(levels(meta_data$CellType))

  if (!palette_name %in% rownames(RColorBrewer::brewer.pal.info)) {
    warning(paste0("颜色板 ", palette_name, " 不存在，已切换为 Set1"))
    palette_name <- "Set1"
  }

  max_colors <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
  base_colors <- RColorBrewer::brewer.pal(
    min(max_colors, colour_count),
    palette_name
  )
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
      axis.text.x = ggplot2::element_text(
        angle = x_angle,
        hjust = 0.5,
        vjust = 0.8
      ),
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
      axis.text.x = ggplot2::element_text(
        angle = x_angle,
        hjust = 0.5,
        vjust = 0.8
      ),
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

merge_seurat_objects <- function(obj_list) {
  if (length(obj_list) == 0) {
    stop("没有任何样本通过初步检查，无法合并。")
  }

  if (length(obj_list) == 1) {
    return(obj_list[[1]])
  }

  merge(
    x = obj_list[[1]],
    y = obj_list[-1],
    add.cell.ids = names(obj_list)
  )
}

read_step_rds <- function(path) {
  if (!file.exists(path)) {
    stop("未找到上一步输出的 RDS：", path)
  }
  readRDS(path)
}

save_step_rds <- function(obj, path) {
  saveRDS(obj, path)
  message("已保存 RDS: ", path)
  invisible(path)
}

assign_group_metadata <- function(obj) {
  obj$samples <- obj$orig.ident
  obj$group <- dplyr::case_when(
    obj$orig.ident %in% c("MIA_LYH_N", "MIA_FXX_N", "MIA_HLQ_N") ~ "Normal",
    obj$orig.ident %in% c("MIA_HXL_T", "MIA_LYP_T", "MIA_HLQ_T") ~ "Microinvasive",
    obj$orig.ident %in% c("MIA_HXL_P", "MIA_LYP_P", "MIA_HLQ_P") ~ "Microinvasive_Para",
    obj$orig.ident %in% c("MIA_ZGM_T", "MIA_FXX_T", "MIA_LYH_T") ~ "Invasive",
    obj$orig.ident %in% c("MIA_ZGM_P", "MIA_FXX_P", "MIA_LYH_P") ~ "Invasive_Para",
    TRUE ~ NA_character_
  )

  if (anyNA(obj$group)) {
    warning("有样本名没有匹配到 CON/Exp 分组，请检查 orig.ident 的命名规则。")
  }

  obj
}

join_layers_if_needed <- function(obj) {
  if (!"Seurat" %in% loadedNamespaces()) {
    return(obj)
  }

  if (!inherits(obj, "Seurat")) {
    return(obj)
  }

  # 1) 如果是 SCT 流程，优先走 PrepSCTFindMarkers
  if ("SCT" %in% Seurat::Assays(obj)) {
    sct_class <- class(obj[["SCT"]])[1]

    if ("PrepSCTFindMarkers" %in% getNamespaceExports("Seurat")) {
      return(
        tryCatch(
          {
            message("检测到 SCT assay (", sct_class, ")，执行 PrepSCTFindMarkers")
            Seurat::PrepSCTFindMarkers(obj, assay = "SCT", verbose = TRUE)
          },
          error = function(e) {
            message("PrepSCTFindMarkers 跳过：", conditionMessage(e))
            obj
          }
        )
      )
    } else {
      message("PrepSCTFindMarkers 不可用，跳过 SCT 预处理")
      return(obj)
    }
  }

  # 2) 非 SCT 情况下，再尝试 JoinLayers（主要针对 v5 RNA assay）
  if (!"JoinLayers" %in% getNamespaceExports("Seurat")) {
    return(obj)
  }

  if ("RNA" %in% Seurat::Assays(obj)) {
    rna_class <- class(obj[["RNA"]])[1]

    if (inherits(obj[["RNA"]], "Assay5") || inherits(obj[["RNA"]], "StdAssay")) {
      return(
        tryCatch(
          {
            message("检测到 RNA assay (", rna_class, ")，执行 JoinLayers")
            Seurat::JoinLayers(obj, assay = "RNA")
          },
          error = function(e) {
            message("JoinLayers 跳过：", conditionMessage(e))
            obj
          }
        )
      )
    } else {
      message("JoinLayers 跳过：RNA assay 不是 Assay5/StdAssay，而是 ", rna_class)
      return(obj)
    }
  }

  return(obj)
}

init_output_dirs()
