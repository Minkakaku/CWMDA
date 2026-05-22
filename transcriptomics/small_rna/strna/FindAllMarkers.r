rm(list = ls())
# hongfan
suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(patchwork)
    library(dplyr)
    library(SeuratDisk)
    library(purrr)
})
#####################################################
dr <- "/data/gpfs02/wmo/work/HF/06.ST-seq/project"
setwd(dr)
wk <- readRDS("wk2/finally_wk2_te.rds")
DefaultAssay(wk) <- "RNA"
Idents(wk) <- wk$seurat_clusters
wk.markers <- FindAllMarkers(wk,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.5,
    test.use = "wilcox"
)
write.csv(wk.markers, "tmpmarkers.csv")
# wk.markers <- read.csv("wk.markers.csv")
top6 <- wk.markers %>%
    group_by(cluster) %>%
    top_n(
        n = 6,
        wt = avg_log2FC
    )
wkp1 <- DoHeatmap(wk,
    group.by = "seurat_clusters", assay = "SCT",
    features = top6$gene, size = 5.5
)
pdf("wk4/wk4.markers.pdf", width = 15, height = 13)
wkp1
dev.off()
