# hongfan
library(dplyr)
library(Seurat)
rm(list = ls())

dir <- "/data/gpfs02/wmo/work/HF/scRNA-seq/220817ST_hep/SC/RDS"
setwd(dir)
fs <- list.files(path = dir, pattern = ".rds")
df <- list()
# intergrate all data
df <- lapply(fs, function(x) {
    readRDS(x)
})
for (i in seq_len(length(fs))) {
    df[[i]] <- SCTransform(df[[i]], vst.flavor = "v2", verbose = FALSE)
}
features <- SelectIntegrationFeatures(object.list = df, nfeatures = 3000)
df <- PrepSCTIntegration(object.list = df, anchor.features = features)

df_anchors <- FindIntegrationAnchors(
    object.list = df, normalization.method = "SCT",
    anchor.features = features
)
df_combined_sct <- IntegrateData(
    anchorset = df_anchors,
    normalization.method = "SCT"
)
df_combined_sct <- RunPCA(df_combined_sct, verbose = FALSE)
df_combined_sct <- RunUMAP(df_combined_sct,
    reduction = "pca",
    dims = 1:30, verbose = FALSE
)
df_combined_sct <- FindNeighbors(df_combined_sct,
    reduction = "pca", dims = 1:30
)
df_combined_sct <- FindClusters(df_combined_sct, resolution = 0.3)
saveRDS(df_combined_sct, "df_combined_sct.rds")

DimPlot(df_combined_sct, reduction = "umap")
