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
  ## -------------------------------
  ## 0. 基本检查
  ## -------------------------------
  if (!inherits(data, "Seurat")) {
    stop("data 必须是 Seurat 对象")
  }
  if (!sample_col %in% colnames(data@meta.data)) {
    stop(paste0("meta.data 中未找到样本列：", sample_col))
  }
  if (!celltype_col %in% colnames(data@meta.data)) {
    stop(paste0("meta.data 中未找到细胞类型列：", celltype_col))
  }
  
  ## -------------------------------
  ## 1. 提取 meta.data
  ## -------------------------------
  meta_data <- data@meta.data[, c(sample_col, celltype_col)]
  colnames(meta_data) <- c("Sample", "CellType")
  meta_data$CellType <- factor(meta_data$CellType)
  
  ## -------------------------------
  ## 2. 颜色映射（安全版）
  ## -------------------------------
  colourCount <- length(levels(meta_data$CellType))
  
  if (!palette_name %in% rownames(RColorBrewer::brewer.pal.info)) {
    warning(paste0("颜色板 ", palette_name, " 不存在，已切换为 Set1"))
    palette_name <- "Set1"
  }
  
  max_colors <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
  base_colors <- RColorBrewer::brewer.pal(min(max_colors, colourCount), palette_name)
  celltype_colors <- colorRampPalette(base_colors)(colourCount)
  names(celltype_colors) <- levels(meta_data$CellType)
  
  ## -------------------------------
  ## 3. 样本 × 细胞类型计数
  ## -------------------------------
  plot_data <- as.data.frame(
    table(meta_data$Sample, meta_data$CellType)
  )
  colnames(plot_data) <- c("Sample", "CellType", "Number")
  
  ## -------------------------------
  ## 4. 细胞数堆叠柱状图
  ## -------------------------------
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
        angle = x_angle, hjust = 0.5, vjust = 0.8
      ),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA),
      legend.position = legend_position,
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = text_size - 1)
    )
  
  ## -------------------------------
  ## 5. 百分比堆叠柱状图
  ## -------------------------------
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
        angle = x_angle, hjust = 0.5, vjust = 0.8
      ),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA),
      legend.position = legend_position,
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = text_size - 1)
    )
  
  ## -------------------------------
  ## 6. 组合图（共享图例）
  ## -------------------------------
  p_combined <- p_count + p_prop +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = legend_position)
  
  ## -------------------------------
  ## 7. 保存
  ## -------------------------------
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
