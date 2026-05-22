############################################
## Stable scRNA-seq pipeline (final version)
############################################

rm(list = ls())
set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(celda)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(harmony)
  library(scales)
})

## ===============================
## paths & helpers
## ===============================
data_dir <- "matrix"
out_dir  <- "output"

dir.create(out_dir, showWarnings = FALSE)
dir.create(file.path(out_dir, "00_QC"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "01_GlobalClustering"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "02_FilteredClustering"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "03_CD8_Subcluster"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "rds"), recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, filename, width = 8, height = 6) {
  ggsave(paste0(filename, ".png"), p, width = width, height = height, dpi = 300)
  ggsave(paste0(filename, ".pdf"), p, width = width, height = height)
}

## ===============================
## sample-level QC
## ===============================
samples <- list.dirs(data_dir, recursive = FALSE, full.names = FALSE)
obj_list <- list()

for (s in samples) {
  
  message("Processing sample: ", s)
  
  counts <- Read10X(file.path(data_dir, s))
  obj <- CreateSeuratObject(
    counts = counts,
    project = s,
    min.cells = 3,
    min.features = 200
  )
  
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, "^mt-")
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, "^Hb[^(P)]")
  obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
  
  obj <- subset(
    obj,
    subset =
      nFeature_RNA > 500 &
      nFeature_RNA < 6000 &
      percent.mt < 15 &
      percent.hb < 1 &
      log10GenesPerUMI > 0.8
  )
  
  if (ncol(obj) < 1000) next
  
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000)
  
  sce <- as.SingleCellExperiment(obj)
  sce <- scDblFinder(sce, dbr = 0.1)
  
  obj$doublet_class <- sce$scDblFinder.class
  obj$doublet_score <- sce$scDblFinder.score
  
  decont <- decontX(counts(sce))
  obj$Contamination <- decont$contamination
  
  obj <- subset(
    obj,
    subset =
      doublet_class == "singlet"
  )
  
  obj_list[[s]] <- obj
}

## ===============================
## merge
## ===============================
sc <- merge(
  x = obj_list[[1]],
  y = obj_list[-1],
  add.cell.ids = names(obj_list)
)

DefaultAssay(sc) <- "RNA"
saveRDS(sc, file.path(out_dir, "rds/sc_after_QC_merged.rds"))

## ===============================
## QC plots (merged)
## ===============================
p1 <- VlnPlot(
  sc,
  features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.hb","log10GenesPerUMI"),
  ncol = 5,
  pt.size = 0.1
)
save_plot(p1, file.path(out_dir, "00_QC/QC_violin"), 16, 5)

p2 <- FeatureScatter(sc, "nCount_RNA", "nFeature_RNA")
save_plot(p2, file.path(out_dir, "00_QC/QC_scatter_nCount_nFeature"))

p3 <- FeatureScatter(sc, "nCount_RNA", "percent.mt")
save_plot(p3, file.path(out_dir, "00_QC/QC_scatter_nCount_mt"))

p4 <- ggplot(sc@meta.data, aes(orig.ident)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p4, file.path(out_dir, "00_QC/QC_cells_per_sample"), 10, 5)

p5 <- ggplot(sc@meta.data, aes(orig.ident, fill = doublet_class)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p5, file.path(out_dir, "00_QC/QC_doublet_ratio"), 10, 5)

p6 <- VlnPlot(sc, "Contamination", group.by = "orig.ident", pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p6, file.path(out_dir, "00_QC/QC_contamination"), 10, 5)

## ===============================
## global clustering + harmony
## ===============================
sc <- NormalizeData(sc, verbose = FALSE)
sc <- FindVariableFeatures(sc, nfeatures = 2000)

sc <- ScaleData(
  sc,
  vars.to.regress = "percent.mt",
  verbose = FALSE
)

sc <- RunPCA(sc, npcs = 50, verbose = FALSE)

sc <- RunHarmony(
  sc,
  group.by.vars = "orig.ident",
  reduction.use = "pca"
)

sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:30)
sc <- FindClusters(sc, resolution = 0.7)
sc <- RunUMAP(sc, reduction = "harmony", dims = 1:30)

saveRDS(sc, file.path(out_dir, "rds/sc_global.rds"))

## global plots
p_elbow <- ElbowPlot(sc, ndims = 50)
save_plot(p_elbow, file.path(out_dir, "01_GlobalClustering/PCA_elbow"))

p_umap1 <- DimPlot(sc, label = TRUE)
save_plot(p_umap1, file.path(out_dir, "01_GlobalClustering/UMAP_clusters"))

p_umap2 <- DimPlot(sc, group.by = "orig.ident")
save_plot(p_umap2, file.path(out_dir, "01_GlobalClustering/UMAP_by_sample"))

p_feat <- FeaturePlot(sc, features = c("Cd3d","Cd3e"), ncol = 2)
save_plot(p_feat, file.path(out_dir, "01_GlobalClustering/Feature_CD3D_CD8A"), 10, 5)

p_feat2 <- FeaturePlot(sc, features = c("Contamination"))
save_plot(p_feat2, file.path(out_dir, "01_GlobalClustering/Feature_Contamination"), 5, 5)

DefaultAssay(sc) <- "RNA"
sc = JoinLayers(sc)
markers_all <- FindAllMarkers(
  sc,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.25,
  return.thresh = 0.05
)

write.csv(
  markers_all,
  file = file.path(out_dir, "01_GlobalClustering/AllClusters_markers.csv"),
  row.names = FALSE
)

markers_top <- markers_all %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 20) %>%
  ungroup()

write.csv(
  markers_top,
  file = file.path(out_dir, "01_GlobalClustering/AllClusters_top20_markers.csv"),
  row.names = FALSE
)

markers_top <-markers_all %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 5) %>%
  ungroup()

top_genes <- unique(markers_top$gene)

p_dot_global <- DotPlot(
  sc,
  features = top_genes
) + RotatedAxis()

save_plot(
  p_dot_global,
  file.path(out_dir, "01_GlobalClustering/DotPlot_topMarkers"),
  width = 40,
  height = 8
)

## mouse-style (your current naming)
markers_prostate <- list(
  Prostate_Luminal = c("Ar","Krt8","Krt18","Epcam","Tacstd2","Foxa1"),
  Prostate_Basal   = c("Krt5","Krt14","Trp63"),
  NEPC             = c("Chga","Chgb","Syp","Ascl1","Insm1","Tubb3"),
  
  T_Pan   = c("Cd3d","Cd3e","Trac","Lck"),
  T_CD8   = c("Cd8a","Cd8b1","Nkg7","Gzmb","Prf1","Ifng"),
  T_CD4   = c("Cd4","Il7r","Ltb"),
  T_Treg  = c("Foxp3","Il2ra","Ctla4","Ikzf2"),
  T_Exhaust = c("Pdcd1","Lag3","Havcr2","Tigit","Tox"),
  
  NK      = c("Ncr1","Klrk1","Klrd1","Nkg7","Prf1","Gzmb"),
  Myeloid = c("Lyz2","Tyrobp","Fcerg1","Aif1"),
  Mono_Inflam = c("S100a8","S100a9","Ly6c2","Lcn2"),
  TAM_C1QC = c("Apoe","C1qa","C1qb","C1qc","Lpl","Mrc1"),
  DC      = c("Itgax","H2-Ab1","Fcer1a","Clec10a"),
  
  Fibro_CAF = c("Col1a1","Col1a2","Dcn","Lum","Col3a1"),
  Myofibro_Pericyte = c("Acta2","Tagln","Myh11","Rgs5"),
  Endothelial = c("Pecam1","Kdr","Vwf","Emcn"),
  
  Cycling = c("Mki67","Top2a","Cenpf","Ube2c"),
  ISG     = c("Isg15","Ifit1","Ifit3","Irf7","Oasl1"),
  Stress  = c("Fos","Jun","Dusp1","Hspa1a","Hspb1"),
  RBC     = c("Hbb-bs","Hbb-bt","Hba-a1")
)

## Example:
p = DotPlot(sc, features = unique(unlist(markers_prostate))) + RotatedAxis()

save_plot(
  p,
  file.path(out_dir, "01_GlobalClustering/DotPlot_cellMarkers"),
  width = 30,
  height = 8
)

####################
sc$cell_type <- NA 
sc$cell_type[sc$seurat_clusters %in% c("0","3","4","8","9")] <- "Prostate_Luminal" 
sc$cell_type[sc$seurat_clusters %in% c("11","13","16","21")] <- "T_Pan"
sc$cell_type[sc$seurat_clusters %in% c("14")] <- "NK" 
sc$cell_type[sc$seurat_clusters %in% c("1","2","7","10","12")] <- "TAM"
sc$cell_type[sc$seurat_clusters %in% c("18")] <- "Cxcr2_positive TANs" 
sc$cell_type[sc$seurat_clusters %in% c("5")] <- "Trem1_positive TANs" 
sc$cell_type[sc$seurat_clusters %in% c("19")] <- "cDC1" 
sc$cell_type[sc$seurat_clusters %in% c("15","17")] <- "cDC2" 
sc$cell_type[sc$seurat_clusters %in% c("25")] <- "Mast" 
sc$cell_type[sc$seurat_clusters %in% c("22")] <- "Fibro_CAF" 
sc$cell_type[sc$seurat_clusters %in% c("24")] <- "Endothelial"
sc$cell_type[sc$seurat_clusters %in% c("6")] <- "Activated myeloid APCs" 
sc$cell_type[sc$seurat_clusters %in% c("20")] <- "B" 
sc$cell_type[sc$seurat_clusters %in% c("23")] <- "Osteoclast-like macrophages"

p_umap_celltype <- DimPlot(
  sc,
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
)

save_plot(p_umap_celltype, file.path(out_dir, "01_GlobalClustering/UMAP_by_annotations") ,13, 10)
saveRDS(sc, file.path(out_dir, "rds/sc_global_annotations.rds"))
## ===============================
## 03_CD8_Subcluster
##  - T extraction -> CD8 extraction -> CD8 subclustering
## ===============================

dir.create(file.path(out_dir, "03_CD8_Subcluster/00_T_cells"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "03_CD8_Subcluster/01_CD8_extract"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "03_CD8_Subcluster/02_CD8_subcluster"), recursive = TRUE, showWarnings = FALSE)

DefaultAssay(sc) <- "RNA"

## -------------------------------
## Step A: extract T cells (Pan-T)
## -------------------------------
tcell <- subset(sc, subset = cell_type == "T_Pan")

tcell <- NormalizeData(tcell, verbose = FALSE)
tcell <- FindVariableFeatures(tcell, nfeatures = 2000, verbose = FALSE)

tcell <- ScaleData(
  tcell,
  vars.to.regress = c("nCount_RNA", "percent.mt"),
  verbose = FALSE
)

tcell <- RunPCA(tcell, npcs = 30, verbose = FALSE)
tcell <- RunUMAP(tcell, dims = 1:20)

tcell <- AddModuleScore(
  tcell,
  features = list(
    CD4 = c("Cd4","Il7r","Ccr7","Ltb"),
    CD8 = c("Cd8a","Cd8b1","Gzmb","Nkg7")
  ),
  name = c("CD4","CD8")
)

tcell$T_type <- "Unassigned"

tcell$T_type[tcell$CD41 > tcell$CD82] <- "CD4"
tcell$T_type[tcell$CD82 > tcell$CD41] <- "CD8"

tcell$T_type <- factor(
  tcell$T_type,
  levels = c("CD4","CD8")
)

p_umap_type <- DimPlot(
  tcell,
  group.by = "T_type",
  label = FALSE
)

save_plot(
  p_umap_type,
  file.path(out_dir, "03_CD8_Subcluster/00_T_cells/UMAP_T_CD4_CD8"),
  8, 6
)

p_score <- FeaturePlot(
  tcell,
  features = c("CD41","CD82"),
  ncol = 2
)

save_plot(
  p_score,
  file.path(out_dir, "03_CD8_Subcluster/00_T_cells/Feature_CD4_CD8_scores"),
  10, 4
)

############################
## Robust CD8 extraction
## （已在上游做过 Harmony，这里不再重复）
############################

## 只保留已经判定为 CD8 的 T 细胞
cd8 <- subset(
  tcell,
  subset = T_type == "CD8"
)

## sanity check
table(cd8$T_type)
ncol(cd8)

dir.create(file.path(out_dir, "rds"), showWarnings = FALSE, recursive = TRUE)
saveRDS(cd8, file.path(out_dir, "rds/cd8_raw_from_tcell.rds"))

## CD8 marker QC
p_cd8_qc <- VlnPlot(
  cd8,
  features = c("Cd8a","Cd8b1","Cd4","Il7r","Nkg7","Gzmb","Trac"),
  ncol = 7,
  pt.size = 0
)

dir.create(
  file.path(out_dir, "03_CD8_Subcluster/01_CD8_extract"),
  recursive = TRUE,
  showWarnings = FALSE
)

save_plot(
  p_cd8_qc,
  file.path(out_dir, "03_CD8_Subcluster/01_CD8_extract/Vln_CD8_markerQC"),
  18, 4
)

############################
## CD8 re-subclustering
## （不再使用 Harmony，保留真实生物学差异）
############################

dir.create(
  file.path(out_dir, "03_CD8_Subcluster/02_CD8_subcluster"),
  recursive = TRUE,
  showWarnings = FALSE
)

DefaultAssay(cd8) <- "RNA"

## 统一预处理（防止对象状态不一致）
cd8 <- NormalizeData(cd8, verbose = FALSE)

cd8 <- FindVariableFeatures(
  cd8,
  nfeatures = 2000,
  verbose = FALSE
)

## 回归测序深度与线粒体比例
cd8 <- ScaleData(
  cd8,
  vars.to.regress = c("nCount_RNA", "percent.mt"),
  verbose = FALSE
)

## PCA（不做 Harmony）
cd8 <- RunPCA(
  cd8,
  npcs = 30,
  verbose = FALSE
)

## 邻接图 / 聚类 / UMAP
cd8 <- FindNeighbors(
  cd8,
  reduction = "pca",
  dims = 1:30,
  verbose = FALSE
)

## 分辨率先用中等偏高，后续根据 marker 再收敛
cd8 <- FindClusters(
  cd8,
  resolution = 0.7,
  verbose = FALSE
)

cd8 <- RunUMAP(
  cd8,
  reduction = "pca",
  dims = 1:30,
  verbose = FALSE
)

saveRDS(
  cd8,
  file.path(out_dir, "rds/cd8_resubcluster_noHarmony.rds")
)

## -------------------------------
## Plots (PNG+PDF)
## -------------------------------
p_cd8_umap <- DimPlot(cd8, label = TRUE, repel = TRUE)
save_plot(p_cd8_umap, file.path(out_dir, "03_CD8_Subcluster/02_CD8_subcluster/UMAP_CD8_clusters_resub"), 10, 7)

p_cd8_by_sample <- DimPlot(cd8, group.by = "orig.ident")
save_plot(p_cd8_by_sample, file.path(out_dir, "03_CD8_Subcluster/02_CD8_subcluster/UMAP_CD8_by_sample_resub"), 9, 6)

# 关键状态marker（用于快速判断：naive-like/effector/exhausted/prolif）
p_cd8_states <- FeaturePlot(
  cd8,
  features = c("Tcf7","Il7r","Ltb","Gzmk","Gzmb","Prf1","Nkg7","Ccl5","Pdcd1","Lag3","Havcr2","Tox","Mki67"),
  ncol = 5
)
save_plot(p_cd8_states, file.path(out_dir, "03_CD8_Subcluster/02_CD8_subcluster/Feature_CD8_states_resub"), 20, 10)


## -------------------------------
## FindAllMarkers (CD8) + remove TAM contamination cluster
## -------------------------------

DefaultAssay(cd8) <- "RNA"
cd8 <- JoinLayers(cd8)

## make a stable cluster id (avoid numeric/character mismatch)
if ("seurat_clusters" %in% colnames(cd8@meta.data)) {
  cd8$cluster_id <- as.character(cd8$seurat_clusters)
} else {
  cd8$cluster_id <- as.character(Idents(cd8))
}

## --- markers ---
cd8_markers_all <- FindAllMarkers(
  cd8,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

dir.create(file.path(out_dir, "03_CD8_Subcluster/02_CD8_subcluster"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "rds"), showWarnings = FALSE, recursive = TRUE)

write.csv(
  cd8_markers_all,
  file = file.path(out_dir, "03_CD8_Subcluster/02_CD8_subcluster/CD8Clusters_markers_resub.csv"),
  row.names = FALSE
)

cd8_markers_top20 <- cd8_markers_all %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::ungroup()

write.csv(
  cd8_markers_top20,
  file = file.path(out_dir, "03_CD8_Subcluster/02_CD8_subcluster/CD8Clusters_top20_markers_resub.csv"),
  row.names = FALSE
)

cd8_markers_top5 <- cd8_markers_all %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup()

cd8_top_genes <- unique(cd8_markers_top5$gene)

p_cd8_dot_top <- DotPlot(cd8, features = cd8_top_genes) + RotatedAxis()
save_plot(
  p_cd8_dot_top,
  file.path(out_dir, "03_CD8_Subcluster/02_CD8_subcluster/DotPlot_CD8_topMarkers_resub"),
  24, 6
)

top10_each <- cd8_markers_all %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()

p_cd8_heat <- DoHeatmap(cd8, features = unique(top10_each$gene), size = 3) + NoLegend()
save_plot(
  p_cd8_heat,
  file.path(out_dir, "03_CD8_Subcluster/02_CD8_subcluster/Heatmap_CD8_top10_resub"),
  14, 10
)



## --- annotate subtypes (robust to cluster_id type) ---
cd8$cd8_cell_type <- NA_character_

cd8$cd8_cell_type[cd8$cluster_id %in% c("0")] <- "Teff CD8"
cd8$cd8_cell_type[cd8$cluster_id %in% c("1")] <- "ISG high CD8"
cd8$cd8_cell_type[cd8$cluster_id %in% c("2")] <- "Xcl1 high Teff CD8"
cd8$cd8_cell_type[cd8$cluster_id %in% c("3")] <- "TAM"
cd8$cd8_cell_type[cd8$cluster_id %in% c("4")] <- "Tpro CD8"
cd8$cd8_cell_type[cd8$cluster_id %in% c("5")] <- "Tex CD8"
cd8$cd8_cell_type[cd8$cluster_id %in% c("6")] <- "Treg CD8"
## --- remove TAM contamination cluster ---
tam_clusters <- c("3")
cd8_clean <- subset(cd8, subset = !cluster_id %in% tam_clusters)

## re-level after removal
cd8_clean$cd8_cell_type <- factor(
  cd8_clean$cd8_cell_type,
  levels = c("Tex CD8", "Tpro CD8", "Teff CD8", "ISG high CD8", "Xcl1 high Teff CD8", "Treg CD8")
)

saveRDS(cd8_clean, file.path(out_dir, "rds/cd8_clean_noTAM.rds"))

## -------------------------------
## Step 0: 准备 CD8 的 sample / celltype 列
## -------------------------------
cd8_clean$sample   <- cd8_clean$orig.ident
cd8_clean$celltype <- cd8_clean$cd8_cell_type

## -------------------------------
## Step 1: 明确分组
## -------------------------------
cd8_clean$group <- NA_character_
cd8_clean$group[cd8_clean$sample %in% c("Control1", "Control2", "Control3")] <- "Control"
cd8_clean$group[cd8_clean$sample %in% c("MC_1", "MC_2", "MC_3")] <- "MC"
cd8_clean$group <- factor(cd8_clean$group, levels = c("Control", "MC"))
stopifnot(!any(is.na(cd8_clean$group)))
library(scplotter)
p_cd8_umap <- CellDimPlot(
  cd8_clean,
  group_by = "cd8_cell_type",pt_size=3,highlight_size = 3,
  highlight = TRUE, theme = "theme_blank"
)
save_plot(
  p_cd8_umap,
  file.path(out_dir, "03_CD8_Subcluster/02_CD8_subcluster/UMAP_CD8_clean_celltypes"),
  8, 8
)

p_cd8_umap <- CellDimPlot(cd8_clean,
            group_by = "cd8_cell_type", facet_by = "group", pt_size=3,highlight_size = 3,
            highlight = TRUE, theme = "theme_blank"
)
save_plot(
  p_cd8_umap,
  file.path(out_dir, "03_CD8_Subcluster/02_CD8_subcluster/UMAP_CD8_clean_celltypes_by_group"),
  10, 10
)
### 堆叠

library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)



dir.create(file.path(out_dir, "03_CD8_Subcluster/03_CD8_proportion"), showWarnings = FALSE, recursive = TRUE)

## -------------------------------
## Step 2: 按样本（用你的 function）
## -------------------------------
plot_celltype_distribution(
  data = cd8_clean,
  sample_col = "sample",
  celltype_col = "celltype",
  palette_name = "Paired",
  output_dir = file.path(out_dir, "03_CD8_Subcluster/03_CD8_proportion"),
  output_file = "CD8_subtype_count_and_percentage_by_sample.pdf",
  width = 14,
  height = 6,
  text_size = 13,
  x_angle = 0,
  legend_position = "bottom"
)

## -------------------------------
## Step 3: 按组汇总（不依赖 uncount）
## 3.1 组内：细胞数（直接汇总 meta.data）
## -------------------------------
group_count <- cd8_clean@meta.data %>%
  dplyr::count(group, celltype, name = "n_cells")

## 3.2 组内：百分比
group_prop <- group_count %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(
    proportion = n_cells / sum(n_cells)
  ) %>%
  dplyr::ungroup()

write.csv(
  group_count,
  file = file.path(out_dir, "03_CD8_Subcluster/03_CD8_proportion/CD8_subtype_counts_by_group.csv"),
  row.names = FALSE
)
write.csv(
  group_prop,
  file = file.path(out_dir, "03_CD8_Subcluster/03_CD8_proportion/CD8_subtype_percentage_by_group.csv"),
  row.names = FALSE
)

## -------------------------------
## Step 4: 组水平组合图（数量 + 百分比）
## -------------------------------
p_group_count <- ggplot(
  group_count,
  aes(x = group, y = n_cells, fill = celltype)
) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(x = "", y = "Cell count", title = "CD8 subtypes (group counts)") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

p_group_percent <- ggplot(
  group_prop,
  aes(x = group, y = proportion, fill = celltype)
) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(x = "", y = "Cell proportion (%)", title = "CD8 subtypes (group percentage)") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

p_group_combined <- p_group_count / p_group_percent +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

save_plot(
  p_group_combined,
  file.path(out_dir, "03_CD8_Subcluster/03_CD8_proportion/Combined_Barplot_CD8_Count_and_Percentage_by_group"),
  8, 8
)
saveRDS(cd8_clean, file.path(out_dir, "rds/cd8_clean_final.rds"))


