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
source(file.path(script_dir, "seurat_subcluster_common.R"), local = FALSE)

set_pipeline_seed(SEEDS$annotation)

obj <- read_required_rds(RDS_CLUSTER)
obj <- prepare_subcluster_input(obj)
obj$cluster_id <- as.character(Seurat::Idents(obj))

if (!file.exists(RDS_MARKERS)) {
  warning("尚未发现 marker 结果，建议先运行 02_subcluster_find_markers.R 再注释。")
}

# 在这里填写 cluster 到亚群类型的映射。
# 先跑完 02_subcluster_find_markers.R，看 markers_all_clusters.csv，再回来补这里。
ANNOTATION <- c(
  # "0" = "Subtype_A",
  # "1" = "Subtype_B"
)

if (length(ANNOTATION) == 0) {
  stop("ANNOTATION 还是空的。请先根据 marker 结果填写 cluster 注释。")
}

if (!all(unique(obj$cluster_id) %in% names(ANNOTATION))) {
  stop("ANNOTATION 没有覆盖所有 cluster，请先补全映射。")
}

sub_cell_type <- ANNOTATION[as.character(obj$cluster_id)]
names(sub_cell_type) <- colnames(obj)
obj$sub_cell_type <- sub_cell_type

write.csv(
  data.frame(
    cluster_id = names(ANNOTATION),
    sub_cell_type = unname(ANNOTATION),
    row.names = NULL
  ),
  file.path(SUB_DIR, "subcluster_annotation_mapping.csv"),
  row.names = FALSE,
  quote = FALSE
)

saveRDS(obj, RDS_ANNOTATED)

save_plot(
  Seurat::DimPlot(
    obj,
    group.by = "sub_cell_type",
    label = TRUE,
    repel = TRUE
  ),
  file.path(SUB_DIR, "UMAP_by_sub_cell_type"),
  9,
  7
)

save_plot(
  Seurat::DimPlot(
    obj,
    group.by = "sub_cell_type",
    split.by = "group",
    ncol = 3
  ),
  file.path(SUB_DIR, "UMAP_split_by_group"),
  15,
  7
)

save_plot(
  Seurat::DimPlot(
    obj,
    group.by = "sub_cell_type",
    split.by = "samples",
    ncol = 3
  ),
  file.path(SUB_DIR, "UMAP_split_by_samples"),
  9,
  7
)

message("亚群注释已保存到 ", RDS_ANNOTATED, "。后续改注释只需要改这个脚本里的 ANNOTATION。")
