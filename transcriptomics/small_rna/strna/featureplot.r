rm(list = ls())
# hongfan
suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(patchwork)
    library(dplyr)
    library(SeuratDisk)
})
#####################################################
dr <- "/data/gpfs02/wmo/work/HF/06.ST-seq/project"
setwd(dr)
wk <- readRDS("wk2/finally_wk2_te.rds")
Idents(wk) <- wk$seurat_clusters
p1 <- DimPlot(wk, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(wk, label = TRUE, label.size = 3)
p3 <- p1 + p2
ggsave("wk2.FeaturePlot.pdf", plot = p3, width = 25, height = 15)





pairs <- c(
    "SPP1", "CD44", "VEGFA", "NRP1", "NRP2",
    "PSAP", "GPR37L1", "EPHB2", "FLT1", "KDR"
)

setwd("/data/gpfs02/wmo/work/HF/06.ST-seq/project/wk2")
for (i in pairs) {
    p <- SpatialFeaturePlot(wk, features = i, alpha = c(0.1, 1))
    ggsave(paste0(i, "wk2.FeaturePlot.pdf"), plot = p)
}
