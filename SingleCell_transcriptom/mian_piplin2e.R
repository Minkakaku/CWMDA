# 清除系统环境的干扰，开始前跑一次就可
rm(list = ls())
# 加载需要用到包
# suppressPackageStartupMessages()是为了不看到包加载时输出很多不需要的输出，可以不放
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(celda)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(harmony)
  library(Startrac)
  library(pheatmap)
  library(Seurat)
  library(ggplot2)
})
# 这个是数据放的目录，格式是XX样本/filtered_feature_bc_matrix/然后三个文件
data_dir <- "00_matrix"
# 创建一下要保存的目录，我这里设置了01 02 03 04，也就是先做到亚分群这一步，后面有需要再新建目录
dir.create("01_QC", recursive = TRUE, showWarnings = FALSE)
dir.create("02_GlobalClustering", recursive = TRUE, showWarnings = FALSE)
dir.create("03_Subcluster", recursive = TRUE, showWarnings = FALSE)
dir.create("04_rds", recursive = TRUE, showWarnings = FALSE)
# 接下来是function，函数就是把需要重复运行的代码整合到一起，可以多次调用
# 这个是保存图片的函数，同时保存png和pdf
save_plot <- function(p, filename, width = 8, height = 6) {
  ggsave(paste0(filename, ".png"), p, width = width, height = height, dpi = 300)
  ggsave(paste0(filename, ".pdf"), p, width = width, height = height)
}
# 这个是批量绘制QC的图，可以自己添加
make_qc_plots <- function(obj) {
  list(
    violin = VlnPlot(
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

    scatter_gene_umi = FeatureScatter(obj, "nCount_RNA", "nFeature_RNA"),

    scatter_mt = FeatureScatter(obj, "nCount_RNA", "percent.mt"),

    cells_per_sample = ggplot(obj@meta.data, aes(orig.ident)) +
      geom_bar() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)),

    doublet_ratio = ggplot(
      obj@meta.data,
      aes(orig.ident, fill = doublet_class)
    ) +
      geom_bar(position = "fill") +
      scale_y_continuous(labels = scales::percent) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)),

    contamination = VlnPlot(
      obj,
      "Contamination",
      group.by = "orig.ident",
      pt.size = 0
    ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
}
# 这个是保存QC图片，先调用前面的make_qc_plots，再save
save_qc <- function(obj, stage) {
  plots <- make_qc_plots(obj)

  for (nm in names(plots)) {
    save_plot(
      plots[[nm]],
      paste0("01_QC/QC_", nm, "_", stage),
      ifelse(nm == "violin", 16, 10),
      ifelse(nm == "violin", 5, 5)
    )
  }
}

# 这个是保存比例图的函数
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


# 检查样本在不在
samples <- list.dirs(data_dir, recursive = FALSE, full.names = FALSE)
if (!dir.exists(data_dir)) {
  stop(paste("错误：数据目录不存在！目录路径：", data_dir))
}

print(paste0("process data in sample:", samples))

obj_list <- list()
# 批量读取然后处理，这里用了一个for循环，计算了线粒体比例，污染数等等
for (s in samples) {
  message("Processing sample: ", s)

  counts <- Read10X(file.path(
    data_dir,
    s,
    "Outs",
    "filtered_feature_bc_matrix"
  ))

  obj <- CreateSeuratObject(
    counts = counts,
    project = s,
    min.cells = 3,
    min.features = 200
  )

  ## -------------------------------
  ## QC metrics
  ## -------------------------------
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, "^mt-")
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, "^Hb[^(P)]")
  obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)

  if (ncol(obj) < 1000) {
    next
  }

  ## -------------------------------
  ## normalize + HVG
  ## -------------------------------
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000)

  ## -------------------------------
  ## doublet detection
  ## -------------------------------
  sce <- as.SingleCellExperiment(obj)
  sce <- scDblFinder(sce)

  obj$doublet_class <- sce$scDblFinder.class
  obj$doublet_score <- sce$scDblFinder.score

  ## -------------------------------
  ## ambient RNA estimation
  ## -------------------------------
  decont <- decontX(counts(sce))
  obj$Contamination <- decont$contamination

  obj_list[[s]] <- obj
}
# 处理完后合并成一个
sc <- merge(
  x = obj_list[[1]],
  y = obj_list[-1],
  add.cell.ids = names(obj_list)
)
# 在过滤之前先看看样本的原始状态
save_qc(sc, stage = "before")
# 这里设定了了过滤的阈值，以及去除多胞
sc <- subset(
  sc,
  subset = nFeature_RNA > 500 &
    nFeature_RNA < 7500 &
    percent.mt < 8 &
    percent.hb < 1 &
    nCount_RNA < 60000 &
    log10GenesPerUMI > 0.8
)

sc <- subset(
  sc,
  subset = doublet_class == "singlet"
)
# 这里去除了环境RNA污染，只有在分群比较杂的情况下用
sc = sc[, sc$Contamination < 0.1]
save_qc(sc, stage = "after")
saveRDS(sc, "04_rds/01.sc_after_QC_merged.rds")

sc = readRDS("04_rds/01.sc_after_QC_merged.rds")
options(future.globals.maxSize = 20 * 1024^3)
sc <- SCTransform(
  sc,
  vars.to.regress = c("percent.mt"),
  verbose = FALSE
)

sc <- RunPCA(sc, npcs = 50, verbose = FALSE)

# 运行 Harmony 整合，由于前面用了 SCT，此时的 DefaultAssay 已经是 "SCT"
sc <- RunHarmony(
  sc,
  group.by.vars = "orig.ident",
  reduction.use = "pca",
  assay.use = "SCT" # 显式指定使用 SCT assay 以保证后续分析调用正确的矩阵
)

## PCA 结构
p_elbow <- ElbowPlot(sc, ndims = 50)
save_plot(p_elbow, "02_GlobalClustering/PCA_elbow")

# 保持基于 harmony 的邻接图和聚类降维不变
sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:25)
sc <- FindClusters(sc, resolution = 0.35)
sc <- RunUMAP(sc, reduction = "harmony", dims = 1:25)

# sc = readRDS("04_rds/01.sc_after_QC_merged.rds")
# # 开始标准流程
# sc <- NormalizeData(sc, verbose = FALSE)
# sc <- FindVariableFeatures(sc, nfeatures = 2000)
# s.genes <- stringr::str_to_title(cc.genes$s.genes)
# g2m.genes <- stringr::str_to_title(cc.genes$g2m.genes)
# # sc <- CellCycleScoring(
# #   sc,
# #   s.features = s.genes,
# #   g2m.features = g2m.genes,
# #   set.ident = TRUE
# # )
# # sc <- ScaleData(
# #   sc,
# #   vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
# #   verbose = FALSE
# # )

# sc <- ScaleData(
#   sc,
#   vars.to.regress = c("percent.mt"),
#   verbose = FALSE
# )
# sc <- RunPCA(sc, npcs = 50, verbose = FALSE)
# sc <- RunHarmony(
#   sc,
#   group.by.vars = "orig.ident",
#   reduction.use = "pca"
# )
# ## PCA 结构
# p_elbow <- ElbowPlot(sc, ndims = 50)
# save_plot(p_elbow, "02_GlobalClustering/PCA_elbow")

# sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:24)
# sc <- FindClusters(sc, resolution = 0.3)
# sc <- RunUMAP(sc, reduction = "harmony", dims = 1:24)

sc$samples = sc$orig.ident
table(sc$samples)
sc$group <- dplyr::case_when(
  grepl("^Con", sc$orig.ident) ~ "CON",
  grepl("^Exp", sc$orig.ident) ~ "Exp",
  TRUE ~ NA_character_
)
table(sc$orig.ident, sc$group)


## 1️⃣ 主分群（cluster）
p_umap_cluster <- DimPlot(
  sc,
  label = TRUE,
  repel = TRUE
)
save_plot(p_umap_cluster, "02_GlobalClustering/UMAP_clusters")

## 2️⃣ 样本分布（orig.ident）
p_umap_sample <- DimPlot(
  sc,
  group.by = "samples"
)
save_plot(p_umap_sample, "02_GlobalClustering/UMAP_by_sample")

## 4️⃣ contamination#这个在骨关节炎这种不干净的数据里用，也可以用来判断肿瘤
p_feat2 <- FeaturePlot(sc, features = "Contamination")
save_plot(
  p_feat2,
  "02_GlobalClustering/Feature_Contamination",
  width = 5,
  height = 5
)

# 找marker基因
DefaultAssay(sc) <- "RNA"

## 如果你用的是 Seurat v5 layered assay
sc <- JoinLayers(sc)

markers_all <- FindAllMarkers(
  sc,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.25,
  return.thresh = 0.05
)

markers_all <- markers_all %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE)

write.csv(
  markers_all,
  file = "02_GlobalClustering/AllClusters_markers.csv",
  row.names = FALSE
)

pdac_markers_broad <- list(
  # 1. 上皮/肿瘤细胞
  "Epithelial" = c("Krt18", "Sox9", "Krt8"),
  # 2. 泛成纤维细胞
  "Fibroblast" = c("Col1a1", "Pdgfra", "Col3a1", "Dcn"),
  # 3. 髓系免疫细胞
  "Myeloid" = c("Ms4a7", "Cd68", "Lyz2"),
  # 4. 淋巴系免疫细胞
  "Lymphocyte" = c("Ptprc", "Cd3d", "Cd79a", "Nkg7"),
  # 5. 内皮细胞
  "Endothelial" = c("Pecam1", "Cdh5", "Kdr")
)

# FeaturePlot(sc,c("Cd4","Cd8a"))

p = DotPlot(sc, features = pdac_markers_broad) + RotatedAxis()
p
save_plot(
  p,
  "02_GlobalClustering/cell_markers",
  10
)
# 在这里替换
sc$cell_type <- NA
sc$cell_type[sc$seurat_clusters %in% c("0")] <- "CAF"
sc$cell_type[sc$seurat_clusters %in% c("1")] <- "Macrophage"
sc$cell_type[sc$seurat_clusters %in% c("2")] <- "CAF"
sc$cell_type[sc$seurat_clusters %in% c("3")] <- "Epithelial"
sc$cell_type[sc$seurat_clusters %in% c("4")] <- "pEpithelial"
sc$cell_type[sc$seurat_clusters %in% c("5")] <- "T"
sc$cell_type[sc$seurat_clusters %in% c("6")] <- "ISG+ Epithelial"
sc$cell_type[sc$seurat_clusters %in% c("7")] <- "Macrophage"
sc$cell_type[sc$seurat_clusters %in% c("8")] <- "VSMCs"
sc$cell_type[sc$seurat_clusters %in% c("9")] <- "Neutrophil"
sc$cell_type[sc$seurat_clusters %in% c("10")] <- "T"
sc$cell_type[sc$seurat_clusters %in% c("11")] <- "DC"
sc$cell_type[sc$seurat_clusters %in% c("12")] <- "B"
sc$cell_type[sc$seurat_clusters %in% c("13")] <- "Plasma"
sc$cell_type[sc$seurat_clusters %in% c("14")] <- "Mesothelial"
sc$cell_type[sc$seurat_clusters %in% c("15")] <- "iNeutrophil"
sc$cell_type[sc$seurat_clusters %in% c("16")] <- "pEpithelial"
sc$cell_type[sc$seurat_clusters %in% c("17")] <- "Endothelial"
DimPlot(
  sc,
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
)

# 画图+计算细胞比例
p_umap_celltype <- DimPlot(
  sc,
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
)
save_plot(
  p_umap_celltype,
  "02_GlobalClustering/UMAP_by_celltype"
)

p_umap_by_sample <- DimPlot(
  sc,
  group.by = "cell_type",
  split.by = "samples",
  ncol = 3
)
save_plot(
  p_umap_by_sample,
  "02_GlobalClustering/UMAP_split_by_sample",
  width = 15,
  height = 8
)

p_umap_by_group <- DimPlot(
  sc,
  group.by = "cell_type",
  split.by = "group",
  ncol = 3
)
save_plot(
  p_umap_by_group,
  "02_GlobalClustering/UMAP_split_by_group",
  width = 15,
  height = 6
)

plot_celltype_distribution(
  data = sc,
  sample_col = "samples",
  celltype_col = "cell_type",
  palette_name = "Paired",
  output_dir = "02_GlobalClustering",
  output_file = "sc_proportion_samples.pdf",
  width = 6,
  height = 6,
  text_size = 7,
  x_angle = 45,
  legend_position = "bottom"
)

plot_celltype_distribution(
  data = sc,
  sample_col = "group",
  celltype_col = "cell_type",
  palette_name = "Paired",
  output_dir = "02_GlobalClustering",
  output_file = "sc_proportion_group.pdf",
  width = 6,
  height = 6,
  text_size = 7,
  x_angle = 45,
  legend_position = "bottom"
)

sample_col <- "orig.ident" # 样本列
celltype_col <- "cell_type" # 细胞类型列（可换成 cell_type）

meta <- sc@meta.data
count_mat <- table(
  meta[[sample_col]],
  meta[[celltype_col]]
)

count_mat <- as.data.frame.matrix(count_mat)


write.csv(
  count_mat,
  file = "02_GlobalClustering/sample_by_celltype_count.csv",
  quote = FALSE
)

saveRDS(sc, "04_rds/01.sc_annoted.rds")
seu_obj <- sc

data_use <- seu_obj@meta.data

# 计算Ro/e
Roe <- calTissueDist(
  data_use,
  byPatient = F, # FALSE: 以组别（tissue/Group）为单位计算总体富集度,TRUE: 为每个患者单独计算富集度（适合患者异质性分析）
  colname.patient = "orig.ident", # 样本名称
  colname.cluster = "cell_type", # 细胞类型
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

save_plot(p, filename = "02_GlobalClustering/ROE")
