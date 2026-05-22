## written by hongfan
## R4.1.3 HPC
## 20220616
## for GBM scRNAseq
rm(list = ls())
## packages
library(Seurat)
library(tidyverse)
library(patchwork)
# library(harmony)
library(SingleR)
library(celldex)
library(limma)
library(BiocParallel)
## path
path <- "/data/gpfs02/wmo/work/HF/scPractice/GBMrecurrent"
setwd(path)
list.files()
###############################################################################
#################################  1.prepare  #################################
## 1 prepare counts
# d1
d1 <- read.csv("Richards_NatureCancer_GBM_scRNAseq_counts.csv",
    row.names = 1
)
d1 <- d1[c(grep("G967|G983", colnames(d1)))]

# d2 needs QC
d2 <- Read10X("/data/gpfs02/wmo/work/HF/scPractice/GBMrecurrent/gnmjk2")
f <- read.csv("GSE173278_scRNAseq_all_cells_metadata.csv", header = T)
f1 <- f[f$source == "tis" & (f$patient == "JK136" | f$patient == "JK202" | f$patient == "JK142" | f$patient == "JK196"), ]
f1_cn <- f1[, 1]
d2 <- d2[, f1_cn]
d2 <- avereps(d2)

# d3
d3 <- read.delim("GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv",
    sep = "\t",
    row.names = 1
)
d3 <- d3[c(grep("MGH85|BT786|BT1187", colnames(d3)))]
########################################################################
## build resuratobject
d1s <- CreateSeuratObject(d1, min.cells = 3, min.features = 200, names.delim = "-")
d2s <- CreateSeuratObject(d2, min.cells = 3, min.features = 200, names.delim = "-")
d3s <- CreateSeuratObject(d3, min.cells = 3, min.features = 200, names.delim = "-")

##########  2. Calculate QC metrics ##########
tmpname <- c(d1s, d2s, d3s)
symbol <- c("group1", "group2", "group3")
for (i in seq_len(length(tmpname))) {
    tmpname[[i]][["percent.mt"]] <- PercentageFeatureSet(
        object = tmpname[[i]],
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
d.integrated <- RunPCA(
    object = d.integrated,
    verbose = FALSE
)
d.integrated <- RunUMAP(
    object = d.integrated,
    dims = 1:30
)
saveRDS(d.integrated, file = "0617integrated.RDS")
plots <- DimPlot(d.integrated,
    group.by = "group",
    combine = FALSE
)
pdf("combine.pdf")
plots <- lapply(X = plots, FUN = function(x) {
    x + theme(legend.position = "top") + guides(color = guide_legend(
        nrow = 4,
        byrow = TRUE, override.aes = list(size = 2.5)
    ))
})
CombinePlots(plots)
dev.off()
###############################################################################
##############################   4.DE analysis   ##############################
d.integrated.markers <- FindAllMarkers(
    object = d.integrated,
    min.pct = 0.25,
    logfc.threshold = 0.5
)
d.filtered.markers <- d.integrated.markers[
    (abs(as.numeric(as.vector(d.integrated.markers$avg_log2FC))) > 0.5 &
        as.numeric(as.vector(d.integrated.markers$p_val_adj)) < 0.05)
]
write.csv(d.filtered.markers, file = "markersgenes.csv")
top10 <- d.filtered.markers %>% top_n(n = 10, wt = avg_log2FC)

pdf(file = "umapHeatmap.pdf", width = 48, height = 36)
DoHeatmap(object = d.integrated, features = top10$gene) + NoLegend()
dev.off()

###############################################################################
################################   5.singleR   ################################
immuneref1 <- readRDS("immuneref1.rds")
immuneref2 <- readRDS("immuneref2.rds")
immuneref3 <- readRDS("immuneref3.rds")
immuneref4 <- readRDS("immuneref4.rds")
singler_labels <- SingleR(
    test = d.integrated,
    ref = list(f1 = immuneref1, f2 = immuneref2, f3 = immuneref3, f4 = immuneref4),
    labels = list(f1$label.main, f2$label.main, f3$label.main, f4$label.main)
)

pdf(file = "celllabel.pdf", width = 10, height = 10)
plotScoreHeatmap(singler_labels,
    clusters = singler_labels@rownames,
    fontsize.row = 9,
    show_colnames = T
)
dev.off()

# new.cluster.ids <- singler_labels$pruned.labels
# names(new.cluster.ids) <- levels(d.integrated)
# d.integrated <- RenameIdents(d.integrated, new.cluster.ids)
# pdf(file = "uMAP_label.pdf", width = 6.5, height = 6)
# DimPlot(d.integrated, reduction = "umap", pt.size = 0.5, label = TRUE)
# dev.off()

d.integrated@meta.data$labels <- singler_labels$labels
pdf(file = "celllabel.pdf")
DimPlot(pbmc3, group.by = c("seurat_clusters", "labels"), reduction = "umap")
dev.off()
