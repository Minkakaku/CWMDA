rm(list = ls())
# hongfan
suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(patchwork)
    library(dplyr)
    library(SeuratDisk)
    library(purrr)
    library(tibble)
})
#####################################################
dr <- "/data/gpfs02/wmo/work/HF/06.ST-seq/project"
setwd(dr)
wk <- readRDS("wk4.te.2.rds")
get_conserved <- function(cluster) {
    FindConservedMarkers(wk,
        ident.1 = cluster,
        grouping.var = "orig.ident",
        only.pos = TRUE
    ) %>%
        rownames_to_column(var = "gene") %>%
        cbind(cluster_id = cluster, .)
}

conserved_markers <- map_dfr(
    c(0:(length(levels(wk$seurat_clusters)) - 1)),
    get_conserved
)
write.csv(conserved_markers, "wk4/wk4.conserved.markers.csv")
# conserved_markers <- read.csv("wk4/wk4.conserved.markers.csv")
top6 <- conserved_markers %>%
    group_by(cluster_id) %>%
    top_n(
        n = 6,
        wt = wk.te.2_avg_log2FC
    )
wkp1 <- DoHeatmap(wk,
    group.by = "seurat_clusters",
    features = top6$gene, size = 5.5
)
pdf("wk4/wk4.conserved.marker.pdf", width = 15, height = 13)
wkp1
dev.off()
