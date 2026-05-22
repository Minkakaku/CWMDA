`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) y else x
}

if (!exists("script_dir")) {
  script_dir <- dirname(normalizePath(sys.frame(1)$ofile %||% ".", mustWork = FALSE))
}

env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

# Default project layout differs by platform:
# - Linux: pipeline scripts live inside the project working directory.
# - Windows: data2analysis is resolved relative to the current working directory.
SCRIPT_DIR <- normalizePath(script_dir, winslash = "/", mustWork = FALSE)
DEFAULT_PROJECT_DIR <- if (identical(.Platform$OS.type, "windows")) {
  getwd()
} else {
  file.path(SCRIPT_DIR, "..")
}
PROJECT_DIR <- normalizePath(env("PROJECT_DIR", DEFAULT_PROJECT_DIR), winslash = "/", mustWork = FALSE)
setwd(PROJECT_DIR)

env_num <- function(name, default) {
  value <- env(name, "")
  if (!nzchar(value)) return(default)
  out <- suppressWarnings(as.numeric(value))
  if (is.na(out)) {
    stop("Invalid numeric environment variable ", name, ": ", value, call. = FALSE)
  }
  out
}

env_int <- function(name, default) {
  value <- env(name, "")
  if (!nzchar(value)) return(as.integer(default))
  out <- suppressWarnings(as.numeric(value))
  if (is.na(out) || out != as.integer(out)) {
    stop("Invalid integer environment variable ", name, ": ", value, call. = FALSE)
  }
  as.integer(out)
}

env_bool <- function(name, default) {
  value <- tolower(env(name, ""))
  if (!nzchar(value)) return(default)
  value %in% c("1", "true", "t", "yes", "y")
}

env_csv <- function(name, default) {
  value <- env(name, "")
  if (!nzchar(value)) return(default)
  value <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  value[nzchar(value)]
}

normalization_method <- function(value) {
  method <- gsub("[^A-Z0-9]+", "", toupper(value))

  if (method %in% c("SCT", "SCTRANSFORM")) {
    return("SCT")
  }

  if (method %in% c("NORMALIZEDATA", "LOGNORMALIZE", "LOGNORM", "RNA")) {
    return("NormalizeData")
  }

  stop("Unsupported normalization method: ", value, ". Use NormalizeData or SCT.", call. = FALSE)
}

batch_correction_method <- function(value) {
  method <- gsub("[^A-Z0-9]+", "", toupper(value))

  if (method %in% c("NONE", "NO", "FALSE", "0", "PCA", "RAW")) {
    return("none")
  }

  if (method %in% c("HARMONY")) {
    return("harmony")
  }

  if (method %in% c("CCA", "CANONICALCORRELATIONANALYSIS")) {
    return("cca")
  }

  stop("Unsupported batch correction method: ", value, ". Use none, harmony, or cca.", call. = FALSE)
}

PACKAGES <- c(
  "Seurat",
  "SingleCellExperiment",
  "scDblFinder",
  "celda",
  "dplyr",
  "ggplot2",
  "scales",
  "RColorBrewer",
  "qs2"
)

OPTIONAL_PACKAGES <- c(
  "ggalluvial",
  "harmony",
  "anndata",
  "Matrix",
  "infercnv",
  "AnnoProbe"
)

DATA_DIR <- env("DATA_DIR", "data2analysis")
QC_DIR <- "01_QC"
GLOBAL_DIR <- "02_GlobalClustering"
MARKER_DIR <- "03_FindMarkers"
ANNOTATION_DIR <- "04_CellAnnotation"
PROPORTION_DIR <- "05_CellProportion"
QS2_DIR <- "04_qs2"
SAMPLE_SHEET_PATH <- env("SAMPLE_SHEET", file.path(QC_DIR, "sample_groups.csv"))
QC_FILTER_SHEET_PATH <- env("QC_FILTER_SHEET", file.path(QC_DIR, "sample_qc_thresholds.csv"))

QS2_PRE_QC <- env("QS2_PRE_QC", file.path(QS2_DIR, "00.sc_before_QC_merged.qs2"))
QS2_QC <- env("QS2_QC", file.path(QS2_DIR, "01.sc_after_QC_merged.qs2"))
QS2_CLUSTER <- env("QS2_CLUSTER", file.path(QS2_DIR, "02.sc_global_clustered.qs2"))
QS2_MARKERS <- env("QS2_MARKERS", file.path(QS2_DIR, "03.all_markers.qs2"))
QS2_ANNOTATED <- env("QS2_ANNOTATED", file.path(QS2_DIR, "04.sc_annotated.qs2"))
QS2_FINAL <- env("QS2_FINAL", file.path(QS2_DIR, "05.sc_final.qs2"))

SEED <- 1001L

QC_READ <- list(
  create_min_cells = 3,
  create_min_features = 200,
  skip_sample_cells_below = 1000,
  run_doublet = TRUE,
  run_decontx = FALSE
)
options(future.globals.maxSize = 24 * 1024^3)  # 4 GiB
# NormalizeData
PIPELINE_NORMALIZATION <- normalization_method(env(
  "PIPELINE_NORMALIZATION",
  "SCT"
))
PIPELINE_ASSAY <- if (identical(PIPELINE_NORMALIZATION, "SCT")) "SCT" else "RNA"

NORMALIZE <- list(
  method = PIPELINE_NORMALIZATION,
  assay = env("PIPELINE_ASSAY", PIPELINE_ASSAY),
  source_assay = env("PIPELINE_SOURCE_ASSAY", "RNA"),
  nfeatures = env_int("PIPELINE_NFEATURES", 3000L),
  regress_vars = env_csv("PIPELINE_REGRESS_VARS", "percent.mt"),
  scale_features = env("PIPELINE_SCALE_FEATURES", "variable"),
  sct_method = env("PIPELINE_SCT_METHOD", ""),
  sct_vst_flavor = env("PIPELINE_SCT_VST_FLAVOR", "v2"),
  sct_return_only_var_genes = env_bool("PIPELINE_SCT_RETURN_ONLY_VAR_GENES", TRUE),
  sct_conserve_memory = env_bool("PIPELINE_SCT_CONSERVE_MEMORY", TRUE),
  verbose = TRUE
)

CLUSTER <- list(
  npcs = env_int("PIPELINE_NPCS", 50L),
  dims = seq_len(env_int("PIPELINE_NDIMS", 35L)),
  resolution = env_num("PIPELINE_RESOLUTION", 0.6),
  batch_method = batch_correction_method(env("PIPELINE_BATCH_METHOD", "harmony")),
  batch_group = env("PIPELINE_BATCH_GROUP", env("PIPELINE_HARMONY_GROUP", "orig.ident")),
  cca_reduction = env("PIPELINE_CCA_REDUCTION", "integrated.cca")
)

GROUP_LEVELS <- c(
  "PBS",
  "TPF",
  "BiTEs@TPF",
  "aCD20@TPF"
)

MARKER_PARAMS <- list(
  assay = NORMALIZE$assay,
  only_pos = env_bool("MARKER_ONLY_POS", TRUE),
  test_use = env("MARKER_TEST_USE", "wilcox"),
  min_pct = env_num("MARKER_MIN_PCT", 0.25),
  logfc_threshold = env_num("MARKER_LOGFC_THRESHOLD", 0.25),
  return_thresh = env_num("MARKER_RETURN_THRESH", 0.05),
  top_n_for_template = env_int("MARKER_TOP_N_FOR_TEMPLATE", 10L)
)

BROAD_MARKER_GENES <- list(
  Epithelial = c("EPCAM", "KRT18", "KRT19"),
  Endothelial = c("PECAM1", "CLDN5", "FLT1"),
  Fibroblast = c("COL1A1", "DCN", "PDGFRA"),
  Myeloid = c("TYROBP", "AIF1", "LST1"),
  T_cell = c("CD3D", "CD3E", "TRAC"),
  NK = c("NKG7", "KLRD1", "GNLY"),
  B_cell = c("CD79A", "MS4A1", "CD74"),
  Plasma = c("JCHAIN", "MZB1", "SDC1"),
  Mast = c("MS4A2", "KIT", "CPA3"),
  Pericyte = c("RGS5", "PDGFRB", "CSPG4"),
  SMC = c("ACTA2", "TAGLN", "MYH11")
)

ANNOTATION <- list(
  template_path = env("ANNOTATION_TEMPLATE", file.path(MARKER_DIR, "cluster_celltype_template.csv")),
  mapping_path = env("ANNOTATION_MAPPING", file.path(ANNOTATION_DIR, "cluster_celltype_mapping.csv")),
  remove_clusters = env_csv("ANNOTATION_REMOVE_CLUSTERS", character(0)),
  rerun_neighbors_clusters = env_bool("ANNOTATION_RERUN_NEIGHBORS_CLUSTERS", FALSE),
  rerun_umap = env_bool("ANNOTATION_RERUN_UMAP", FALSE),
  rerun_resolution = env_num("ANNOTATION_RERUN_RESOLUTION", CLUSTER$resolution),
  rerun_dims = seq_len(env_int("ANNOTATION_RERUN_NDIMS", max(CLUSTER$dims))),
  cluster_backup_col = env("ANNOTATION_CLUSTER_BACKUP_COL", "global_cluster_before_annotation_cleanup")
)

PROPORTION_COMPARISONS <- list(
  c("PBS", "TPF"),
  c("PBS", "BiTEs@TPF"),
  c("PBS", "aCD20@TPF"),
  c("TPF", "BiTEs@TPF"),
  c("TPF", "aCD20@TPF"),
  c("BiTEs@TPF", "aCD20@TPF")
)

PROPORTION <- list(
  comparisons = PROPORTION_COMPARISONS,
  make_sankey = env_bool("PROPORTION_MAKE_SANKEY", TRUE),
  sankey_min_count = env_int("PROPORTION_SANKEY_MIN_COUNT", 1L)
)

SUBCLUSTER_NORMALIZATION <- normalization_method(env("SUBCLUSTER_NORMALIZATION", NORMALIZE$method))
SUBCLUSTER_ASSAY <- if (identical(SUBCLUSTER_NORMALIZATION, "SCT")) "SCT" else "RNA"

SUBCLUSTER <- list(
  input_qs2 = env("SUBCLUSTER_INPUT_QS2", QS2_ANNOTATED),
  parent_col = env("SUBCLUSTER_PARENT_COL", "cell_type"),
  target_celltype = env_csv("SUBCLUSTER_TARGET", c("Treg", "Cycling Treg")),
  sample_col = env("SUBCLUSTER_SAMPLE_COL", "samples"),
  group_col = env("SUBCLUSTER_GROUP_COL", "group"),
  method = SUBCLUSTER_NORMALIZATION,
  assay = env("SUBCLUSTER_ASSAY", SUBCLUSTER_ASSAY),
  source_assay = env("SUBCLUSTER_SOURCE_ASSAY", "RNA"),
  nfeatures = env_int("SUBCLUSTER_NFEATURES", 3000L),
  regress_vars = env_csv("SUBCLUSTER_REGRESS_VARS", "percent.mt"),
  scale_features = env("SUBCLUSTER_SCALE_FEATURES", "variable"),
  sct_method = env("SUBCLUSTER_SCT_METHOD", NORMALIZE$sct_method),
  sct_vst_flavor = env("SUBCLUSTER_SCT_VST_FLAVOR", NORMALIZE$sct_vst_flavor),
  sct_return_only_var_genes = env_bool(
    "SUBCLUSTER_SCT_RETURN_ONLY_VAR_GENES",
    NORMALIZE$sct_return_only_var_genes
  ),
  sct_conserve_memory = env_bool("SUBCLUSTER_SCT_CONSERVE_MEMORY", NORMALIZE$sct_conserve_memory),
  verbose = TRUE,
  npcs = env_int("SUBCLUSTER_NPCS", 15L),
  dims = seq_len(env_int("SUBCLUSTER_NDIMS", 15L)),
  resolution = env_num("SUBCLUSTER_RESOLUTION", 0.2),
  batch_method = batch_correction_method(env("SUBCLUSTER_BATCH_METHOD", "harmony")),
  batch_group = env("SUBCLUSTER_BATCH_GROUP", env("SUBCLUSTER_HARMONY_GROUP", "orig.ident")),
  cca_reduction = env("SUBCLUSTER_CCA_REDUCTION", "subcluster.integrated.cca"),
  remove_clusters = env_csv("SUBCLUSTER_REMOVE_CLUSTERS", character(0)),
  rerun_neighbors_clusters = env_bool("SUBCLUSTER_RERUN_NEIGHBORS_CLUSTERS", FALSE),
  rerun_umap = env_bool("SUBCLUSTER_RERUN_UMAP", FALSE),
  rerun_resolution = env_num("SUBCLUSTER_RERUN_RESOLUTION", env_num("SUBCLUSTER_RESOLUTION", 0.2)),
  rerun_dims = seq_len(env_int("SUBCLUSTER_RERUN_NDIMS", env_int("SUBCLUSTER_NDIMS", 15L))),
  cluster_backup_col = env("SUBCLUSTER_CLUSTER_BACKUP_COL", "subcluster_before_annotation_cleanup"),
  marker_params = list(
    only_pos = env_bool("SUBCLUSTER_MARKER_ONLY_POS", TRUE),
    test_use = env("SUBCLUSTER_MARKER_TEST_USE", "wilcox"),
    min_pct = env_num("SUBCLUSTER_MARKER_MIN_PCT", 0.25),
    logfc_threshold = env_num("SUBCLUSTER_MARKER_LOGFC_THRESHOLD", 0.25),
    return_thresh = env_num("SUBCLUSTER_MARKER_RETURN_THRESH", 0.05),
    top_n_for_template = env_int("SUBCLUSTER_MARKER_TOP_N_FOR_TEMPLATE", 10L)
  ),
  template_path = env("SUBCLUSTER_TEMPLATE", file.path("06_Subcluster", "subcluster_annotation_template.csv")),
  comparisons = PROPORTION_COMPARISONS,
  make_sankey = env_bool("SUBCLUSTER_MAKE_SANKEY", TRUE)
)

SUBCLUSTER_NAME <- gsub("[^A-Za-z0-9_]+", "_", paste(SUBCLUSTER$target_celltype, collapse = "_"))
if (!nzchar(SUBCLUSTER_NAME)) SUBCLUSTER_NAME <- "AllCells"
SUBCLUSTER_DIR <- file.path("06_Subcluster", SUBCLUSTER_NAME)

QS2_SUBCLUSTER_CLUSTERED <- env(
  "QS2_SUBCLUSTER_CLUSTERED",
  file.path(QS2_DIR, paste0("06.", SUBCLUSTER_NAME, ".subcluster_clustered.qs2"))
)
QS2_SUBCLUSTER_MARKERS <- env(
  "QS2_SUBCLUSTER_MARKERS",
  file.path(QS2_DIR, paste0("07.", SUBCLUSTER_NAME, ".subcluster_markers.qs2"))
)
QS2_SUBCLUSTER_ANNOTATED <- env(
  "QS2_SUBCLUSTER_ANNOTATED",
  file.path(QS2_DIR, paste0("08.", SUBCLUSTER_NAME, ".subcluster_annotated.qs2"))
)

SUBCLUSTER$template_path <- env(
  "SUBCLUSTER_TEMPLATE",
  file.path(SUBCLUSTER_DIR, "subcluster_annotation_template.csv")
)

dir.create(QC_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(GLOBAL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MARKER_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(ANNOTATION_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PROPORTION_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(SUBCLUSTER_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(QS2_DIR, showWarnings = FALSE, recursive = TRUE)
