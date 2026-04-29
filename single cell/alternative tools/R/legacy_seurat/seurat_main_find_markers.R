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

human_lung_markers_broad_core <- list(
  Epithelial   = c("EPCAM", "KRT18", "KRT19"),
  Endothelial  = c("PECAM1", "CLDN5", "FLT1"),
  Fibroblast   = c("COL1A1", "DCN", "PDGFRA"),
  Myeloid      = c("TYROBP", "AIF1", "LST1"),
  Macrophage   = c("CD68", "C1QA", "APOE"),
  Monocyte     = c("FCN1", "CCR2", "CTSS"),
  Neutrophil   = c("S100A8", "S100A9", "CSF3R"),
  Dendritic    = c("FCER1A", "ITGAX", "HLA-DRA"),
  T_cell       = c("CD3D", "CD3E", "TRAC"),
  NK           = c("NKG7", "KLRD1", "GNLY"),
  B_cell       = c("CD79A", "MS4A1", "CD74"),
  Plasma       = c("JCHAIN", "MZB1", "SDC1"),
  Mast         = c("MS4A2", "KIT", "CPA3"),
  Pericyte     = c("RGS5", "PDGFRB", "CSPG4"),
  SMC          = c("ACTA2", "TAGLN", "MYH11")
)

human_lung_markers_broad_core <- lapply(
  human_lung_markers_broad_core,
  function(x) x[x %in% rownames(sc)]
)
human_lung_markers_broad_core <- human_lung_markers_broad_core[
  lengths(human_lung_markers_broad_core) > 0
]

p_markers <- Seurat::DotPlot(
  sc,
  features = human_lung_markers_broad_core,
  assay = "SCT"
) + Seurat::RotatedAxis()

save_plot( p_markers, file.path(MARKER_DIR, "cell_markers"), width = 15, height = 6 )

message("FindAllMarkers 完成。请先查看 ", file.path(MARKER_DIR, "AllClusters_markers.csv"))
