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
wk <- readRDS("wk1.te.2.rds")
p1 <- DimPlot(wk, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(wk, label = TRUE, label.size = 3)
p3 <- p1 + p2
ggsave("wk1_feature_plots.pdf", plot = p3, width = 20, height = 15)
wk <- FindSpatiallyVariableFeatures(wk,
    assay = "SCT",
    features = VariableFeatures(wk)[1:2000],
    selection.method = "markvariogram"
)
top.features <- head(
    SpatiallyVariableFeatures(wk,
        selection.method = "markvariogram"
    ),
    6
)
wkp1 <- DoHeatmap(wk,
    group.by = "seurat_clusters",
    features = top.features, size = 5.5
)
pdf("wkp1.pdf", width = 20, height = 15)
wkp1
dev.off()
SpatialFeaturePlot(wk, features = top.features, ncol = 3, alpha = c(0.1, 1))
ggsave("wkSpatialFeaturePlot.pdf", width = 30, height = 15)
wk.markers <- FindAllMarkers(wk,
    only.pos = FALSE,
    min.pct = 0.25,
    logfc.threshold = 0.5,
    test.use = "wilcox"
)
# write.csv(wk.markers, "wk.markers.csv")
wkp1 <- DoHeatmap(wk,
    group.by = "seurat_clusters",
    features = wk.markers, size = 5.5
)
pdf("wkFindAllMarkersp2.pdf", width = 20, height = 15)
wkp1
dev.off()
get_conserved <- function(seurat_clusters) {
    FindConservedMarkers(wk,
        ident.1 = seurat_clusters,
        grouping.var = "sample",
        only.pos = TRUE
    ) %>%
        rownames_to_column(var = "gene") %>%
        left_join(
            y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name")
        ) %>%
        cbind(cluster_id = cluster, .)
}
conserved_markers <- map_dfr(
    c(0:0:(length(levels(wk$seurat_clusters)) - 1)),
    get_conserved
)
wkp1 <- DoHeatmap(wk,
    group.by = "seurat_clusters",
    features = top.features, size = 5.5
)
pdf("wkconservedp1.pdf", width = 20, height = 15)
wkp1
dev.off()
# write.csv(conserved_markers, "wk.conserved.markers.csv")
