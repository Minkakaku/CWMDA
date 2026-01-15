library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

## ===============================
## paths
## ===============================
mono_dir <- file.path(out_dir, "03_CD8_Subcluster/06_CD8_slingshot")
mono_rds <- file.path(mono_dir, "rds")
mono_fig <- file.path(mono_dir, "figures")

dir.create(mono_rds, recursive = TRUE, showWarnings = FALSE)
dir.create(mono_fig, recursive = TRUE, showWarnings = FALSE)

## ===============================
## load input (CD8 clean)
## ===============================
scRNA <- readRDS(file.path(out_dir, "rds/cd8_clean_final.rds"))

group_col <- "group"

## ===============================
## Slingshot (global fit)
## ===============================
sce <- as.SingleCellExperiment(scRNA)
reducedDims(sce)$UMAP <- scRNA@reductions$umap@cell.embeddings
sce$cluster <- factor(scRNA$cd8_cell_type)

sce <- slingshot(
  sce,
  clusterLabels = "cluster",
  reducedDim    = "UMAP",
  start.clus    = "Tpro CD8"
)

saveRDS(sce, file.path(mono_rds, "cd8_sce_slingshot.rds"))

## pseudotime
pt <- slingPseudotime(sce)
sce$pseudotime_slingshot <- pt[, 1]
scRNA$pseudotime_slingshot <- sce$pseudotime_slingshot
saveRDS(scRNA, file.path(mono_rds, "cd8_with_slingshot_pseudotime.rds"))

## ===============================
## prepare UMAP dataframe
## ===============================
umap_df <- as.data.frame(scRNA@reductions$umap@cell.embeddings)
colnames(umap_df) <- c("UMAP1", "UMAP2")

umap_df$cluster    <- as.character(scRNA$cd8_cell_type)
umap_df$group      <- as.character(scRNA@meta.data[[group_col]])
umap_df$pseudotime <- scRNA$pseudotime_slingshot
umap_df$in_lineage <- !is.na(umap_df$pseudotime)

## ===============================
## slingshot curves
## ===============================
curve_list <- slingCurves(sce)

curve_df <- do.call(
  rbind,
  lapply(seq_along(curve_list), function(i) {
    df <- as.data.frame(curve_list[[i]]$s)
    colnames(df) <- c("UMAP1", "UMAP2")
    df$lineage <- paste0("Lineage_", i)
    df
  })
)

lineage_cols <- setNames(
  brewer.pal(max(3, length(unique(curve_df$lineage))), "Set1"),
  unique(curve_df$lineage)
)

## 颜色（固定顺序，避免跨图变化）
celltype_levels <- unique(umap_df$cluster)
celltype_cols <- setNames(
  brewer.pal(max(3, length(celltype_levels)), "Set2")[seq_along(celltype_levels)],
  celltype_levels
)

p_combined <- ggplot() +
  
  ## 所有细胞：CD8 subtype
  geom_point(
    data = umap_df,
    aes(
      x = UMAP1,
      y = UMAP2,
      fill = cluster
    ),
    shape  = 21,
    size   = 3.5,      # 类似 pt_size
    stroke = 0.3,      # 轻描边，避免太“糊”
    color  = "black"
  ) +
  
  ## Slingshot trajectories（全局一致）
  geom_path(
    data = curve_df,
    aes(
      x = UMAP1,
      y = UMAP2,
      group = lineage
    ),
    color = "black",
    linewidth = 1.1
  ) +
  
  scale_fill_manual(
    values = celltype_cols,
    name = "CD8 subtype"
  ) +
  
  facet_wrap(~ group, nrow = 1) +
  
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "CD8 subtypes with Slingshot trajectories"
  )

ggsave(
  file.path(mono_fig, "CD8_Subtype_Slingshot_CombinedGroup.png"),
  p_combined,
  width = 16,
  height = 7
)

ggsave(
  file.path(mono_fig, "CD8_Subtype_Slingshot_CombinedGroup.pdf"),
  p_combined,
  width = 16,
  height = 7
)
p_pseudotime <- ggplot() +
  
  ## 所有细胞：pseudotime 填充
  geom_point(
    data = umap_df,
    aes(
      x = UMAP1,
      y = UMAP2,
      fill = pseudotime
    ),
    shape  = 21,
    size   = 3.5,      # 和 subtype 图一致
    stroke = 0.3,
    color  = "black"
  ) +
  
  ## Slingshot trajectories
  geom_path(
    data = curve_df,
    aes(
      x = UMAP1,
      y = UMAP2,
      group = lineage
    ),
    color = "black",
    linewidth = 1.1
  ) +
  
  scale_fill_viridis_c(
    option = "plasma",
    na.value = "grey80",
    name = "Slingshot pseudotime"
  ) +
  
  facet_wrap(~ group, nrow = 1) +
  
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "CD8 Slingshot pseudotime"
  )

ggsave(
  file.path(mono_fig, "CD8_Pseudotime_Slingshot_CombinedGroup.png"),
  p_pseudotime,
  width = 16,
  height = 7
)

ggsave(
  file.path(mono_fig, "CD8_Pseudotime_Slingshot_CombinedGroup.pdf"),
  p_pseudotime,
  width = 16,
  height = 7
)