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
wk$seurat_clusters <- Idents(wk)
sublabelcelltype <- c(
    "0" = "Microglial cell", "1" = "Astrocyte",
    "2" = "Astrocyte", "3" = "Neuron",
    "4" = "Astrocyte", "5" = "Oligodendrocyte", "6" = "Neuron",
    "7" = "Red blood cell/Astrocytes", "8" = "Astrocyte/Microglial cell",
    "9" = "Oligodendrocyte precursor cells", "10" = "Oligodendrocyte precursor cells",
    "11" = "Neuron"
)
wk[["seurat_clusters"]] <- unname(sublabelcelltype[wk@meta.data$seurat_clusters])
saveRDS(wk, "wk2/finally_wk2_te.rds")
