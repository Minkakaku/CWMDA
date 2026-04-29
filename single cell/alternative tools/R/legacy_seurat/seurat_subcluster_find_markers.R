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

set_pipeline_seed(SEEDS$markers)

obj <- read_required_rds(RDS_CLUSTER)
Seurat::DefaultAssay(obj) <- PARAM$assay
obj <- join_layers_if_needed(obj)
obj$cluster_id <- as.character(Seurat::Idents(obj))

markers_all <- Seurat::FindAllMarkers(
  obj,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.25
) |>
  dplyr::arrange(cluster, desc(avg_log2FC))

write.csv(
  markers_all,
  file.path(SUB_DIR, "markers_all_clusters.csv"),
  row.names = FALSE
)
saveRDS(markers_all, RDS_MARKERS)

message("marker 已输出：", file.path(SUB_DIR, "markers_all_clusters.csv"))
