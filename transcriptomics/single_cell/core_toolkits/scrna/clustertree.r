## written by hongfan
## R4.1.3 HPC
## 20220616
## for GBM scRNAseq
rm(list = ls())
## packages

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  # library(harmony)
  library(SingleR)
  library(celldex)
  library(limma)
  library(BiocParallel)
  library(clustree)
})
## path
path <- "/data/gpfs02/wmo/work/HF/scPractice/GBMrecurrent"
setwd(path)
list.files()
###############################################################################
#################################  1.prepare  #################################
## 1 prepare counts
# d1
d1 <- read.csv("prim/Richards_NatureCancer_GBM_scRNAseq_counts.csv",
  row.names = 1
)
d1 <- d1[c(grep("G967|G983", colnames(d1)))]
rownames(d1) <- gsub("_", "-", rownames(d1))
colnames(d1) <- gsub("_", "-", colnames(d1))

# d2 needs QC
d2 <- Read10X("prim/gnmjk2")
f <- read.csv("prim/GSE173278_scRNAseq_all_cells_metadata.csv", header = T)
f1 <- f[f$source == "tis" & (f$patient == "JK136" | f$patient == "JK202" | f$patient == "JK142" | f$patient == "JK196"), ]
f1_cn <- f1[, 1]
d2 <- d2[, f1_cn]
d2 <- avereps(d2)
rownames(d2) <- gsub("_", "-", rownames(d2))
colnames(d2) <- gsub("_", "-", colnames(d2))
# m1
m1 <- read.csv("metasis/m1.csv",
  header = T,
  row.names = 1
)
rownames(m1) <- gsub("_", "-", rownames(m1))
colnames(m1) <- gsub("_", "-", colnames(m1))
# m2
m2 <- read.csv("metasis/GSM3516668_MSK_LX255B_METASTASIS_dense.csv",
  header = T, row.names = 1
)
m2 = as.data.frame(t(m2))
rownames(m2) <- gsub("_", "-", rownames(m2))
colnames(m2) <- gsub("_", "-", colnames(m2))
# m3
m3 <- read.csv("metasis/GSM3516671_MSK_LX681_METASTASIS_dense.csv",
  header = T, row.names = 1
)
m3 = as.data.frame(t(m3))
rownames(m3) <- gsub("_", "-", rownames(m3))
colnames(m3) <- gsub("_", "-", colnames(m3))
# m4
m4 <- read.csv("metasis/GSM3516678_MSK_LX701_METASTASIS_dense.csv",
  header = T, row.names = 1
)
m4 = as.data.frame(t(m4))
rownames(m4) <- gsub("_", "-", rownames(m4))
colnames(m4) <- gsub("_", "-", colnames(m4))
########################################################################
## build resuratobject
d1s <- CreateSeuratObject(d1, min.cells = 3, min.features = 200)
d2s <- CreateSeuratObject(d2, min.cells = 3, min.features = 200)
m1s <- CreateSeuratObject(m1, min.cells = 3, min.features = 200)
m2s <- CreateSeuratObject(m2, min.cells = 3, min.features = 200)
m3s <- CreateSeuratObject(m3, min.cells = 3, min.features = 200)
m4s <- CreateSeuratObject(m4, min.cells = 3, min.features = 200)

##########  2. Calculate QC metrics ##########
tmpname <- c(d1s, d2s, m1s, m2s, m3s, m4s)
saveRDS(tmpname, file = "tmpname.rds")
symbol <- c("prim", "prim2", "metas", "metas2", "metas3", "metas4")
for (i in seq_len(length(tmpname))) {
  tmpname[[i]][["percent.mt"]] <- PercentageFeatureSet(
    object = tmpname[[i]],
    assay = "RNA",
    pattern = "^MT-"
  )
  tmpname[[i]]$log10GenesPerUMI <- log10(tmpname[[i]]$nFeature_RNA) / log10(tmpname[[i]]$nCount_RNA)
  # pdf(
  #     file = paste(symbol[i], "featureViolin.pdf", sep = "_"),
  #     width = 10, height = 6
  # )
  # VlnPlot(
  #     object = tmpname[[i]],
  #     features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  #     ncol = 3
  # )
  # dev.off()
  tmpname[[i]] <- subset(
    x = tmpname[[i]],
    subset = nFeature_RNA > 200 &
      percent.mt < 15 &
      nCount_RNA < 30000 &
      log10GenesPerUMI > 0.8
  )
}
for (i in seq_len(length(tmpname))) {
  counts <- GetAssayData(object = tmpname[[i]], slot = "counts")
  nonzero <- counts > 0
  keep_genes <- Matrix::rowSums(nonzero) >= 10
  filtered_counts <- counts[keep_genes, ]
  tmpname[[i]] <- CreateSeuratObject(
    filtered_counts,
    meta.data = tmpname[[i]]@meta.data
  )
  tmpname[[i]]$group <- symbol[i]
}
# for (i in seq_len(length(tmpname))) {
#     tmpname[[i]] <- NormalizeData(tmpname[[i]])
#     tmpname[[i]] <- FindVariableFeatures(
#         tmpname[[i]],
#         selection.method = "vst",
#         nfeatures = 2000
#     )
# }
## SCT
d.list <- list()
for (i in seq_len(length(tmpname))) {
  d.list[[i]] <- SCTransform(tmpname[[i]], verbose = FALSE)
}
###############################################################################
############################  3.standard workflow  ############################
d.features <- SelectIntegrationFeatures(
  object.list = d.list,
  nfeatures = 3000
)
d.list <- PrepSCTIntegration(
  object.list = d.list,
  anchor.features = d.features
)
d.anchors <- FindIntegrationAnchors(
  object.list = d.list,
  normalization.method = "SCT",
  anchor.features = d.features
)
d.integrated <- IntegrateData(
  anchorset = d.anchors,
  normalization.method = "SCT"
)

DefaultAssay(d.integrated) <- "integrated"
saveRDS(d.integrated, file = "20220628results/20220626d.integrated.rds")

d.integrated <- RunPCA(
  object = d.integrated,
  verbose = FALSE
)
## 20PCAnumber
d.integrated <- RunUMAP(
  object = d.integrated,
  dims = 1:20,
  reduction = "pca"
)
d.integrated <- FindNeighbors(
  d.integrated,
  reduction = "pca",
  dims = 1:20
)
# Check clustering stability at given resolution
# Set different resolutions
res.used <- seq(0.1, 1, by = 0.1)
# Loop over and perform clustering of different resolutions
for (i in res.used) {
  sce <- FindClusters(object = d.integrated, verbose = T, resolution = res.used)
}
# Make plot
pdf("clustertree.pdf")
clustree(sce) +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
dev.off()
