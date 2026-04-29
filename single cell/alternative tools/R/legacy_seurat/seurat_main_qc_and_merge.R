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

set_pipeline_seed(PIPELINE_SEEDS$qc)

if (!dir.exists(DATA_DIR)) {
  stop("数据目录不存在：", DATA_DIR)
}

samples <- list.dirs(DATA_DIR, recursive = FALSE, full.names = FALSE)

if (length(samples) == 0) {
  stop("在 ", DATA_DIR, " 下没有找到任何样本目录。")
}

message("待处理样本：", paste(samples, collapse = ", "))

obj_list <- list()

for (sample_idx in seq_along(samples)) {
  s <- samples[[sample_idx]]
  matrix_dir <- file.path(DATA_DIR, s, "outs", "filtered_feature_bc_matrix")

  if (!dir.exists(matrix_dir)) {
    warning("跳过样本 ", s, "：未找到目录 ", matrix_dir)
    next
  }

  message("Processing sample: ", s)

  counts <- Seurat::Read10X(matrix_dir)

  obj <- Seurat::CreateSeuratObject(
    counts = counts,
    project = s,
    min.cells = 3,
    min.features = 200
  )

  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, "^MT-")
  obj[["percent.hb"]] <- Seurat::PercentageFeatureSet(obj, "^HB[^(P)]")
  obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)

  if (ncol(obj) < 1000) {
    message("跳过样本 ", s, "：细胞数少于 1000。")
    next
  }

  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = 2000)

  set_pipeline_seed(PIPELINE_SEEDS$qc + sample_idx)
  sce <- Seurat::as.SingleCellExperiment(obj)
  sce <- scDblFinder::scDblFinder(sce)

  obj$doublet_class <- sce$scDblFinder.class
  obj$doublet_score <- sce$scDblFinder.score

  set_pipeline_seed(PIPELINE_SEEDS$qc + sample_idx)
  decont <- celda::decontX(counts(sce))
  obj$Contamination <- decont$contamination

  obj_list[[s]] <- obj
}

sc <- merge_seurat_objects(obj_list)
save_qc(sc, stage = "before")

sc <- subset(
  sc,
  subset = nFeature_RNA > 500 &
    nFeature_RNA < 7500 &
    percent.mt < 8 &
    percent.hb < 1 &
    nCount_RNA < 60000 &
    log10GenesPerUMI > 0.8
)

sc <- subset(sc, subset = doublet_class == "singlet")
sc <- sc[, !is.na(sc$Contamination) & sc$Contamination < 0.2]

save_qc(sc, stage = "after")
save_step_rds(sc, RDS_QC)

message("QC 完成，保留细胞数：", ncol(sc))
