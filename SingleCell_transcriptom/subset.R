library(Seurat)

seu_obj = readRDS("04_rds\\01.sc_annoted.rds")
seu_target_subset <- subset(
  x = seu_obj, 
  subset = cell_type %in% c("Epithelial", "ISG + Epithelial", "CAF", "pEpithelial", "DC")
)

# 检查提取后的细胞类型分布，确认是否提取成功
table(seu_target_subset$cell_type)

saveRDS(seu_target_subset, "04_rds//copykat.RDS")
