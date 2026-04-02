library(Seurat)
library(anndata)
library(Matrix)

# 假设你的原始对象名为 seurat_obj
sc = readRDS("04_rds\\01.sc_annoted.rds")
# ---------------------------------------------------------
# 第一步：使用 DietSeurat 提取纯净的 RNA Assay
# 注意：如果你想把降维结果（如 umap, pca）一起导过去，必须在 dimreducs 中指定！
# ---------------------------------------------------------
seurat_rna <- DietSeurat(
  sc,
  assays = "RNA",
  dimreducs = c("pca", "umap") # 请根据你需要导出的降维结果名称进行修改
)

# 确保默认 Assay 是 RNA
DefaultAssay(seurat_rna) <- "RNA"

# ---------------------------------------------------------
# 第二步：合并图层 (针对 Seurat v5)
# 将可能存在的 counts.1, counts.2 等图层合并
# ---------------------------------------------------------
seurat_rna <- JoinLayers(seurat_rna)

# 提取表达矩阵 (通常提取 counts 层)
# 注意：Python (AnnData) 要求【细胞为行，基因为列】，所以必须用 t() 进行转置
X_mat <- t(LayerData(seurat_rna, assay = "RNA", layer = "counts"))

# 组装 AnnData 对象
ad <- AnnData(
  X = X_mat,
  obs = seurat_rna@meta.data # 放入 metadata
)

# 提取并装入降维坐标 (前提是你之前跑过)
if ("pca" %in% names(seurat_rna@reductions)) {
  ad$obsm$X_pca <- Embeddings(seurat_rna, "pca")
}
if ("umap" %in% names(seurat_rna@reductions)) {
  ad$obsm$X_umap <- Embeddings(seurat_rna, "umap")
}

# 直接保存为 .h5ad
ad$write_h5ad("04_rds\\01.sc_annoted.h5ad")
