rm(list = ls())
suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(Matrix)
    library(SeuratDisk)
})
path = "/data/gpfs02/wmo/work/HF/06.ST-seq/project"
setwd(path)
wk <- readRDS("wk2/finally_wk2_te.rds")
table(wk$seurat_clusters)
# so <- LoadH5Seurat("final.h5seurat")
# B_cells <- subset(so, idents = "B cells")
# NK_cells <- subset(so, idents = "NK cells")
# B_cells <- so[, so$cell_type %in% "B cells"]
# NK_cells <- so[, so$cell_type %in% "NK cells"]

# # Save normalised counts - NOT scaled!
# writeMM(wk@assays$RNA@data, file = "adult_S/matrix.mtx")
# save gene and cell names
# write(x = rownames(adult_S@assays$RNA@data), file = "adult_S/features.tsv")
# write(x = colnames(adult_S@assays$RNA@data), file = "adult_S/barcodes.tsv")
wk@meta.data$Cell <- rownames(wk@meta.data)
df <- wk@meta.data[, c("Cell", "seurat_clusters")]
write.table(df,
    file = "wk2_meta.txt",
    sep = "\t",
    quote = F,
    row.names = F
)

# # Save normalised counts - NOT scaled!
# writeMM(aged_S@assays$RNA@data, file = "aged_S/matrix.mtx")
# # save gene and cell names
# write(x = rownames(aged_S@assays$RNA@data), file = "aged_S/features.tsv")
# write(x = colnames(aged_S@assays$RNA@data), file = "aged_S/barcodes.tsv")
# aged_S@meta.data$Cell <- rownames(aged_S@meta.data)
# df <- aged_S@meta.data[, c("Cell", "cell_type")]
# write.table(df,
#     file = "aged_meta.txt",
#     sep = "\t",
#     quote = F,
#     row.names = F
# )

SaveH5Seurat(
    wk,
    filename = "wk2/wk2.h5Seurat",
    overwrite = FALSE,
    verbose = TRUE,
)
Convert("wk2/wk2.h5Seurat", "h5ad",
    assay = "RNA",
    overwrite = FALSE,
    verbose = TRUE,
)

# SaveH5Seurat(
#     NK_cells,
#     filename = "NK_cells.h5Seurat",
#     overwrite = FALSE,
#     verbose = TRUE,
# )
# Convert("aged_S.h5seurat", "h5ad",
#     assay = "RNA",
#     overwrite = FALSE,
#     verbose = TRUE,
# )
