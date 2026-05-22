rm(list = ls())
# hongfan
suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(patchwork)
    library(dplyr)
    library(scCATCH)
})
#####################################################
dr <- "/data/gpfs02/wmo/work/HF/06.ST-seq/project"
setwd(dr)
wk <- readRDS("wk3.te.2.rds")
wk.data.cluster = as.character(wk$seurat_clusters)
DefaultAssay(wk) <- "RNA"
wk = NormalizeData(wk)
wk.data <- GetAssayData(wk, slot = "data", assay = "RNA")
wk.data <- rev_gene(data = wk.data, data_type = "data", species = "Human", geneinfo = geneinfo)
obj <- createscCATCH(data = wk.data, cluster = wk.data.cluster)
cellmatch_new <- cellmatch[cellmatch$species == "Human", ]
tissue = c("Brain", "Dorsolateral prefrontal cortex", "Embryonic brain", "Embryonic prefrontal cortex", "Fetal brain", "Hippocampus", "Inferior colliculus", "Midbrain", "Sympathetic ganglion")
cancer = c("Astrocytoma", "CNS Primitive Neuroectodermal Tumor", "Ependymoma", "Glioma", "Glioblastoma", "Head and Neck Cancer", "High-grade glioma", "Intracranial Aneurysm", "Malignant Peripheral Nerve Sheath Tumor", "Medulloblastoma")
cellmatch_new <- cellmatch_new[cellmatch_new$cancer %in% cancer | cellmatch_new$tissue %in% tissue, ]
obj <- findmarkergene(object = obj, if_use_custom_marker = TRUE, marker = cellmatch_new)
obj <- findcelltype(obj)
write.csv(obj@celltype, "wk3.celltype.csv")
