rm(list=ls())
Sys.which("python")
options(python_cmd = "E:\\miniconda3\\python.exe")
BiocManager::install("infercnv")
library(infercnv)
library(Seurat)
system("python --version")
obj <- readRDS("04_rds\\01.sc_annoted.rds")
DefaultAssay(obj) <- "RNA"
mat <- GetAssayData(obj, layer = "counts")
mat <- as.matrix(mat)
mat <- mat[rowSums(mat) > 0, ]

# 假设你有 celltype_major
reference_types <- c("T", "B", "Myeloid", "Endothelial", "Fibroblast")
ref_cells <- colnames(obj)[obj$celltype_major %in% reference_types]
obs_cells <- setdiff(colnames(obj), ref_cells)

# 如果你只想在上皮里判肿瘤 vs 正常：obs 用上皮，ref 用免疫/基质
obs_cells <- colnames(obj)[obj$celltype_major == "Epithelial"]

# 这里给一个兜底：没有分型时，先全部当 observation（不推荐用于最终结论）
ref_cells <- colnames(obj)[obj$celltype_major %in% c("T", "B", "Myeloid")]
obs_cells <- setdiff(colnames(obj), ref_cells)

use_cells <- c(ref_cells, obs_cells)
mat_use <- mat[, use_cells, drop = FALSE]

anno <- data.frame(
    cell = use_cells,
    group = ifelse(use_cells %in% ref_cells, "reference", "observation")
)

# infercnv 的 annotation 文件要求：两列，无表头也行；我建议写成文件
outdir <- "infercnv_out"
dir.create(outdir, showWarnings = FALSE)

anno_file <- file.path(outdir, "cell_annotations.txt")
write.table(
    anno,
    file = anno_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
)