rm(list = ls())

get_current_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)

  if (length(file_arg) > 0) {
    return(dirname(normalizePath(
      gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE),
      winslash = "/",
      mustWork = FALSE
    )))
  }

  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(
      sys.frames()[[1]]$ofile,
      winslash = "/",
      mustWork = FALSE
    )))
  }

  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_dir <- get_current_script_dir()
source(file.path(script_dir, "seurat_main_common.R"), local = FALSE)

set_pipeline_seed(PIPELINE_SEEDS$annotation)

sc <- read_step_rds(RDS_CLUSTER)
sc <- assign_group_metadata(sc)

if (!file.exists(RDS_MARKERS)) {
  warning("未找到 marker 结果：", RDS_MARKERS, "。建议先运行 03_find_markers.R 再做注释。")
}

# 手工注释在这里改最方便，后续如果 cluster 编号变化，只需要更新这个向量。
cluster_to_celltype <- c(
  "0" = "CD4 T",
  "1" = "Extra_alveolar capillary EC",
  "2" = "NK",
  "3" = "COL14A1 matrix FB",
  "4" = "CD8 T",
  "5" = "Neutrophil",
  "6" = "Alveolar capillary EC",
  "7" = "Club Epi",
  "8" = "Resident Macrophage",
  "9" = "Monocyte",
  "10" = "MAST",
  "11" = "Alveolar Macrophage",
  "12" = "Extra_alveolar capillary EC",
  "13" = "Basal Epi",
  "14" = "PTGDS matrix FB",
  "15" = "AT1 Epi",
  "16" = "SMC",
  "17" = "Venous EC",
  "18" = "Pericyte",
  "19" = "B",
  "20" = "Artery EC",
  "21" = "Plasma",
  "22" = "Lymphatic endothelial",
  "23" = "Ciliated",
  "24" = "pDC"
)

sc$cell_type <- unname(cluster_to_celltype[as.character(sc$seurat_clusters)])

if (anyNA(sc$cell_type)) {
  warning("有部分 cluster 尚未被注释，请检查 cluster_to_celltype。")
}

mapping_df <- data.frame(
  cluster = names(cluster_to_celltype),
  cell_type = unname(cluster_to_celltype),
  row.names = NULL
)
write.csv(
  mapping_df,
  file = file.path(ANNOTATION_DIR, "cluster_celltype_mapping.csv"),
  row.names = FALSE,
  quote = FALSE
)

p_umap_celltype <- Seurat::DimPlot(
  sc,
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
)
save_plot(p_umap_celltype, file.path(ANNOTATION_DIR, "UMAP_by_celltype"),  width = 13,
  height = 8)

p_umap_by_sample <- Seurat::DimPlot(
  sc,
  group.by = "cell_type",
  split.by = "samples",
  ncol = 3
)
save_plot(
  p_umap_by_sample,
  file.path(ANNOTATION_DIR, "UMAP_split_by_sample"),
  width = 12,
  height = 10
)

p_umap_by_group <- Seurat::DimPlot(
  sc,
  group.by = "cell_type",
  split.by = "group",
  ncol = 3
)
save_plot(
  p_umap_by_group,
  file.path(ANNOTATION_DIR, "UMAP_split_by_group"),
  width = 12,
  height = 6
)

save_step_rds(sc, RDS_ANNOTATED)
save_step_rds(sc, RDS_ANNOTATED_LEGACY)

message("细胞注释完成。兼容旧流程的 RDS 也已同步更新：", RDS_ANNOTATED_LEGACY)
