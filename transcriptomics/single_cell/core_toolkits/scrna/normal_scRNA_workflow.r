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
wirte.csv(d1, file = "G967|G983_counts.csv")


rownames(d1) <- gsub("_", "-", rownames(d1))
colnames(d1) <- gsub("_", "-", colnames(d1))

# d2 needs QC
d2 <- Read10X("prim/gnmjk2")
f <- read.csv("prim/GSE173278_scRNAseq_all_cells_metadata.csv", header = T)
f1 <- f[f$source == "tis" & (f$patient == "JK136" | f$patient == "JK202" | f$patient == "JK142" | f$patient == "JK196"), ]
f1_cn <- f1[, 1]
d2 <- d2[, f1_cn]
d2 <- avereps(d2)
write.csv(d2, file = "GSE173278_scRNAseq_all_cells_cpunts.csv")

rownames(d2) <- gsub("_", "-", rownames(d2))
colnames(d2) <- gsub("_", "-", colnames(d2))
# m1
m1 <- read.csv("metasis/m1.csv",
    header = T,
    row.names = 1
)
write.csv(m1, file = "test.csv")
rownames(m1) <- gsub("_", "-", rownames(m1))
colnames(m1) <- gsub("_", "-", colnames(m1))
# m2
m2 <- read.csv("metasis/GSM3516668_MSK_LX255B_METASTASIS_dense.csv",
    header = T, row.names = 1
)
m2 = as.data.frame(t(m2))

write.csv(m2, file = "m2.csv")
rownames(m2) <- gsub("_", "-", rownames(m2))
colnames(m2) <- gsub("_", "-", colnames(m2))
# m3
m3 <- read.csv("metasis/GSM3516671_MSK_LX681_METASTASIS_dense.csv",
    header = T, row.names = 1
)
m3 = as.data.frame(t(m3))
write.csv(m3, file = "m3.csv")
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
pdf("elbowplot.pdf")
ElbowPlot(d.integrated)
dev.off()
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
d.integrated <- FindClusters(
    d.integrated,
    resolution = 0.4
)
saveRDS(d.integrated, file = "20220628results/20220626d.alldone.rds")

plots <- DimPlot(d.integrated,
    group.by = "group",
    combine = FALSE
)

pdf("20220628results/0626combine1.pdf")
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
# DefaultAssay(d.integrated) = "RNA"
# d.integrated.markers <- FindAllMarkers(
#     object = d.integrated,
#     min.pct = 0.25,
#     logfc.threshold = 0.5
# )
# d.filtered.markers <- d.integrated.markers[
#     (abs(as.numeric(as.vector(d.integrated.markers$avg_log2FC))) > 0.5 &
#         as.numeric(as.vector(d.integrated.markers$p_val_adj)) < 0.05),
# ]
# write.csv(d.filtered.markers, file = "markersgenes.csv")
# top10 <- d.filtered.markers %>% top_n(n = 10, wt = avg_log2FC)

# pdf(file = "umapHeatmap.pdf", width = 48, height = 36)
# DoHeatmap(object = d.integrated, features = top10$gene) + NoLegend()
# dev.off()

###############################################################################
################################   5.singleR   ################################

immuneref2 = readRDS("/data/gpfs02/wmo/work/HF/scPractice/GBMrecurrent/immuneref2.rds")
testdata <- GetAssayData(d.integrated, slot = "data")
singler_labels <- SingleR(
    test = testdata,
    ref = immuneref2,
    labels = immuneref2$label,
    method = "cluster",
    clusters = d.integrated@meta.data$seurat_clusters
)

pdf(file = "20220628results/0626celllabel1.pdf", width = 30, height = 30)
plotScoreHeatmap(singler_labels,
    clusters = singler_labels@rownames,
    fontsize.row = 9,
    show_colnames = T
)
dev.off()

new.cluster.ids <- singler_labels$pruned.labels
names(new.cluster.ids) <- levels(d.integrated)
d.integrated <- RenameIdents(d.integrated, new.cluster.ids)
pdf(file = "20220628results/0626uMAP_label.pdf", width = 6.5, height = 6)
DimPlot(d.integrated,
    reduction = "umap", pt.size = 0.5, label = TRUE
)
dev.off()
