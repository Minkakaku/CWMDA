rm(list=ls())
set.seed(123456)
## -------------------------------
## 0. Global parameters
## -------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(Startrac)
  library(pheatmap)
})

PARAM <- list(
  assay = "RNA",

  celltype_col = "cell_type", # 主注释列
  target_celltype = "Epithelial", # 当前亚群名称（用于 subset）

  sample_col = "orig.ident",

  use_harmony = FALSE, # <<< 是否做 Harmony
  harmony_group = "orig.ident", # Harmony 的 batch 变量

  regress_vars = c("nCount_RNA", "percent.mt"),
  nfeatures = 2000,
  npcs = 30,
  dims = 1:30,
  resolution = 0.6,

  remove_clusters = character(0), # 需要移除的 cluster id
  group_map = NULL, # sample -> group 映射（可选）
  input = "04_rds/01.sc_annoted.rds"
)

SUBCLUSTER_NAME <- PARAM$target_celltype   # 目录名 = 亚群名

save_plot <- function(p, filename, width = 8, height = 6) {
  ggsave(paste0(filename, ".png"), p, width = width, height = height, dpi = 300)
  ggsave(paste0(filename, ".pdf"), p, width = width, height = height)
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
  
  meta_data <- data@meta.data[, c(sample_col, celltype_col)]
  colnames(meta_data) <- c("Sample", "CellType")
  meta_data$CellType <- factor(meta_data$CellType)
  
  colourCount <- length(levels(meta_data$CellType))
  
  if (!palette_name %in% rownames(RColorBrewer::brewer.pal.info)) {
    warning(paste0("颜色板 ", palette_name, " 不存在，已切换为 Set1"))
    palette_name <- "Set1"
  }
  
  max_colors <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
  base_colors <- RColorBrewer::brewer.pal(min(max_colors, colourCount), palette_name)
  celltype_colors <- colorRampPalette(base_colors)(colourCount)
  names(celltype_colors) <- levels(meta_data$CellType)
  
  plot_data <- as.data.frame(
    table(meta_data$Sample, meta_data$CellType)
  )
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
        angle = x_angle, hjust = 0.5, vjust = 0.8
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
        angle = x_angle, hjust = 0.5, vjust = 0.8
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

## =========================================================
## 1. Directory structure (flat)
## =========================================================
sub_root <- "03_Subcluster"
sub_dir <- file.path(sub_root, SUBCLUSTER_NAME)
dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)

## =========================================================
## 2. Input selection
## =========================================================
sc = readRDS(PARAM$input)
DefaultAssay(sc) <- PARAM$assay

if (!is.null(PARAM$target_celltype)) {
  stopifnot(PARAM$celltype_col %in% colnames(sc@meta.data))
  obj <- subset(sc, subset = cell_type == PARAM$target_celltype)
} else {
  obj <- sc
}

## =========================================================
## 3. Preprocessing (PCA only)
## =========================================================
obj <- obj |>
  NormalizeData(verbose = FALSE) |>
  FindVariableFeatures(
    nfeatures = PARAM$nfeatures,
    verbose = FALSE
  ) |>
  ScaleData(
    vars.to.regress = intersect(
      PARAM$regress_vars,
      colnames(obj@meta.data)
    ),
    verbose = FALSE
  ) |>
  RunPCA(
    npcs = PARAM$npcs,
    verbose = FALSE
  )

## =========================================================
## 4. Batch effect check (decision point)
## =========================================================
obj <- RunUMAP(
  obj,
  reduction = "pca",
  dims = PARAM$dims,
  verbose = FALSE
)

save_plot(
  DimPlot(obj, group.by = PARAM$sample_col, pt.size = 0.2),
  file.path(sub_dir, "UMAP_raw_by_sample"),
  9, 7
)

## =========================================================
## 5. Optional Harmony integration
## =========================================================
if (PARAM$use_harmony) {

  suppressPackageStartupMessages(library(harmony))
  stopifnot(PARAM$harmony_group %in% colnames(obj@meta.data))

  obj <- RunHarmony(
    obj,
    group.by.vars = PARAM$harmony_group,
    reduction = "pca",
    assay.use = PARAM$assay,
    verbose = FALSE
  )

  obj <- obj |>
    FindNeighbors(
      reduction = "harmony",
      dims = PARAM$dims,
      verbose = FALSE
    ) |>
    FindClusters(
      resolution = PARAM$resolution,
      verbose = FALSE
    ) |>
    RunUMAP(
      reduction = "harmony",
      dims = PARAM$dims,
      verbose = FALSE
    )

  saveRDS(obj, file.path(sub_dir, "subcluster_harmony.rds"))

} else {

  obj <- obj |>
    FindNeighbors(
      reduction = "pca",
      dims = PARAM$dims,
      verbose = FALSE
    ) |>
    FindClusters(
      resolution = PARAM$resolution,
      verbose = FALSE
    ) |>
    RunUMAP(
      reduction = "pca",
      dims = PARAM$dims,
      verbose = FALSE
    )

  saveRDS(obj, file.path(sub_dir, "subcluster_no_harmony.rds"))
}

## =========================================================
## 6. UMAP visualization
## =========================================================
save_plot(
  DimPlot(obj, label = TRUE, repel = TRUE),
  file.path(sub_dir, "UMAP_clusters"),
  9,
  7
)

save_plot(
  DimPlot(obj, group.by = PARAM$sample_col),
  file.path(sub_dir, "UMAP_by_sample"),
  9,
  7
)

## =========================================================
## 7. Marker detection
## =========================================================
obj <- JoinLayers(obj)
obj$cluster_id <- as.character(Idents(obj))

markers_all <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.25
) |>
  arrange(cluster, desc(avg_log2FC))

write.csv(
  markers_all,
  file.path(sub_dir, "markers_all_clusters.csv"),
  row.names = FALSE
)


ANNOTATION <- c(
  "0" = "malignant_epithelial",
  "1" = "malignant_epithelial",
  "2" = "malignant_epithelial",
  "3" = "malignant_epithelial",
  "4" = "malignant_epithelial",
  "5" = "malignant_epithelial",
  "6" = "normal_epithelial",
  "7" = "malignant_epithelial",
  "8" = "normal_epithelial",
  "9" = "malignant_epithelial"
)

if (!is.null(ANNOTATION)) {
  stopifnot(all(unique(obj$cluster_id) %in% names(ANNOTATION)))

  sub_cell_type <- ANNOTATION[as.character(obj$cluster_id)]
  names(sub_cell_type) <- colnames(obj)

  obj$sub_cell_type <- sub_cell_type
} else {
  obj$sub_cell_type <- obj$cluster_id
}
table(obj$sub_cell_type)
saveRDS(obj, file.path(sub_dir, "subcluster_annotated.rds"))

p <- DimPlot(
  obj,
  group.by = "sub_cell_type",
  label = TRUE,
  repel = TRUE
)

save_plot(
  p,
  file.path(sub_dir, "UMAP_by_sub_cell_type"),
  9,
  7
)

p <- DimPlot(
  obj,
  group.by = "sub_cell_type",
  split.by = "group",
  ncol = 3
)

save_plot(
  p,
  file.path(sub_dir, "UMAP_split_by_group"),
  15,
  7
)

p <- DimPlot(
  obj,
  group.by = "sub_cell_type",
  split.by = "samples",
  ncol = 3
)

save_plot(
  p,
  file.path(sub_dir, "UMAP_split_by_samples"),
  9,
  7
)

plot_celltype_distribution(
  data = obj,
  sample_col = "samples",
  celltype_col = "sub_cell_type",
  palette_name = "Paired",
  output_dir = sub_dir,
  output_file = "Epithelial_subtype_by_sample.pdf",
  width = 6,
  height = 6,
  text_size = 7,
  x_angle = 45,
  legend_position = "bottom"
)

plot_celltype_distribution(
  data = obj,
  sample_col = "group",
  celltype_col = "sub_cell_type",
  palette_name = "Paired",
  output_dir = sub_dir,
  output_file = "Epithelial_subtype_by_group.pdf",
  width = 6,
  height = 6,
  text_size = 7,
  x_angle = 45,
  legend_position = "bottom"
)

seu_obj <- readRDS("03_Subcluster\\Epithelial\\subcluster_annotated.rds")

data_use <- seu_obj@meta.data
colnames(data_use)
# 计算Ro/e
Roe <- calTissueDist(
  data_use,
  byPatient = F, # FALSE: 以组别（tissue/Group）为单位计算总体富集度,TRUE: 为每个患者单独计算富集度（适合患者异质性分析）
  colname.patient = "orig.ident", # 样本名称
  colname.cluster = "sub_cell_type", # 细胞类型
  colname.tissue = "group", # 样本分组
  method = "chisq", # 计算方法，可选："chisq", "fisher", and "freq"
  min.rowSum = 0 # 行总和过滤阈值，行和（某个细胞类型在所有组中的总细胞数）小于等于此值的细胞类型会被过滤掉
)
Roe

# 绘制热图
p = pheatmap(
  Roe,
  scale = "none",
  display_numbers = TRUE,
  number_color = "black",
  number_format = "%.2f",
  color = colorRampPalette(c("#2166AC", "#b55c67ff", "#B2182B"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = "grey80",
  cellwidth = 30,
  cellheight = 25
)

### 这里需要改
save_plot(p, filename = "03_Subcluster\\Epithelial\\ROE")

