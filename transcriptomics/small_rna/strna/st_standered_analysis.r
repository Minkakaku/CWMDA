rm(list = ls())
# hongfan
suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(patchwork)
    library(dplyr)
})
#####################################################
dr <- "/data/gpfs02/wmo/work/HF/06.ST-seq/project"
setwd(dr)
##################################################### read 4 datasets
# dir <- c(
#     "/data/gpfs02/wmo/work/HF/06.ST-seq/project/WK1",
#     "/data/gpfs02/wmo/work/HF/06.ST-seq/project/WK2",
#     "/data/gpfs02/wmo/work/HF/06.ST-seq/project/WK3",
#     "/data/gpfs02/wmo/work/HF/06.ST-seq/project/WK4"
# )
# names(dir) <- c("wk1", "wk2", "wk3", "wk4")
# ##################################################### edit
# wk <- list()
# for (i in seq_len(length(dir))) {
#     wk[[i]] <- Load10X_Spatial(data.dir = dir[i])
#     wk[[i]]@meta.data$orig.ident <- names(dir)[i]
# }
# for (i in seq_len(length(wk))) {
#     wk[[i]] <- SCTransform(wk[[i]],
#         assay = "Spatial",
#         return.only.var.genes = FALSE,
#         verbose = FALSE
#     )
# }
# wk.merge <- merge(wk[[1]], y = c(wk[[2]], wk[[3]], wk[[4]]))
# DefaultAssay(wk.merge) <- "SCT"
# VariableFeatures(wk.merge) <- c(
#     VariableFeatures(wk[[1]]), VariableFeatures(wk[[2]]),
#     VariableFeatures(wk[[3]]), VariableFeatures(wk[[4]])
# )
# wk.merge <- RunPCA(wk.merge, verbose = FALSE)
# wk.merge <- FindNeighbors(wk.merge, dims = 1:20)
# wk.merge <- FindClusters(wk.merge, verbose = FALSE)
# wk.merge <- RunUMAP(wk.merge, dims = 1:30)
# saveRDS(wk.merge, "wk.merge.rds")
# wk.merge <- readRDS("wk.merge.rds")
# p1 <- DimPlot(wk.merge, reduction = "umap", label = TRUE)
# p2 <- SpatialDimPlot(wk.merge, label = TRUE, label.size = 3)
# p1 + p2
# ggsave("wk.plot.pdf", height = 40, width = 200, limitsize = FALSE)


d <- "/data/gpfs02/wmo/work/HF/06.ST-seq/project/WK1"
wk1 <- Load10X_Spatial(data.dir = d)
wk1 <- SCTransform(wk1,
        assay = "Spatial",
        return.only.var.genes = FALSE,
        verbose = FALSE
    )
DefaultAssay(wk1) <- "SCT"
VariableFeatures(wk1) <- VariableFeatures(wk1)
wk1 <- RunPCA(wk1, verbose = FALSE)
wk1 <- FindNeighbors(wk1, dims = 1:20)
wk1 <- FindClusters(wk1, verbose = FALSE)
wk1 <- RunUMAP(wk1, dims = 1:30)

library(STdeconvolve)
pos <- wk1$slice1@coordinates ## x and y positions of each pixel
cd <- wk1$Spatial@counts ## matrix of gene counts in each pixel
# annot <- wk1$annot ## annotated tissue layers assigned to each pixel
## remove pixels with too few genes
counts <- cleanCounts(counts = cd,
                      min.lib.size = 100,
                      min.reads = 1,
                      min.detected = 1,
                      verbose = TRUE)
## feature select for genes
corpus <- restrictCorpus(counts,
                         removeAbove = 1.0,
                         removeBelow = 0.05,
                         alpha = 0.05,
                         plot = TRUE,
                         verbose = TRUE)
## Note: the input corpus needs to be an integer count matrix of pixels x genes
ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(2, 9, by = 1),
               perc.rare.thresh = 0.05,
               plot=TRUE,
               verbose=TRUE)
## select model with minimum perplexity
optLDA <- optimalModel(models = ldas, opt = "min")

## extract pixel cell-type proportions (theta) and cell-type gene expression profiles (beta) for the given dataset
## we can also remove cell-types from pixels that contribute less than 5% of the pixel proportion
## and scale the deconvolved transcriptional profiles by 1000 
results <- getBetaTheta(optLDA,
                        perc.filt = 0.05,
                        betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta
colnames(pos) = c("tissue" ,  "x"   ,   "y"     , "imagerow", "imagecol")
pos2=pos[,c(2,3)]
pdf("vizAllTopics.pdf")
vizAllTopics(deconProp, pos2, 
            #  groups = annot, 
            #  group_cols = rainbow(length(levels(annot))),
             r=0.4)
dev.off()