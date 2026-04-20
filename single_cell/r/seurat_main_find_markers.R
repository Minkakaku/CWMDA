rm(list = ls())

get_current_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)

  if (length(file_arg) > 0) {
    return(dirname(normalizePath(
      sub("^--file=", "", file_arg[1]),
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

Seurat::Idents(sc) <- "seurat_clusters"
sc <- join_layers_if_needed(sc)
DefaultAssay(sc) <- "SCT"
markers_all <- Seurat::FindAllMarkers(
  sc,
  assay = "SCT",
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.25,
  return.thresh = 0.05
)

markers_all <- markers_all %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(desc(avg_log2FC), .by_group = TRUE)

write.csv(
  markers_all,
  file = file.path(MARKER_DIR, "AllClusters_markers.csv"),
  row.names = FALSE
)
saveRDS(markers_all, RDS_MARKERS)

pdac_markers_broad <- list(
  Epithelial = c("Krt18", "Sox9", "Krt8"),
  Fibroblast = c("Col1a1", "Pdgfra", "Col3a1", "Dcn"),
  Myeloid = c("Ms4a7", "Cd68", "Lyz2"),
  Lymphocyte = c("Ptprc", "Cd3d", "Cd79a", "Nkg7"),
  Endothelial = c("Pecam1", "Cdh5", "Kdr")
)

p_markers <- Seurat::DotPlot(sc, features = pdac_markers_broad) +
  Seurat::RotatedAxis()
save_plot(
  p_markers,
  file.path(MARKER_DIR, "cell_markers"),
  width = 10,
  height = 6
)

message("FindAllMarkers 完成。请先查看 ", file.path(MARKER_DIR, "AllClusters_markers.csv"))
