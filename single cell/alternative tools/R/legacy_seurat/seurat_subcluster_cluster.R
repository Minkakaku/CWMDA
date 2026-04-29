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

set_pipeline_seed(SEEDS$subcluster)

sc <- read_required_rds(PARAM$input)
sc <- prepare_subcluster_input(sc)
Seurat::DefaultAssay(sc) <- PARAM$assay

if (!is.null(PARAM$target_celltype)) {
  if (!PARAM$celltype_col %in% colnames(sc@meta.data)) {
    stop("输入对象缺少主注释列：", PARAM$celltype_col)
  }

  obj <- sc[, sc[[PARAM$celltype_col, drop = TRUE]] == PARAM$target_celltype]
} else {
  obj <- sc
}

if (ncol(obj) == 0) {
  stop("subset 后没有细胞，请检查 target_celltype 设置。")
}

obj <- obj |>
  Seurat::NormalizeData(verbose = FALSE) |>
  Seurat::FindVariableFeatures(
    nfeatures = PARAM$nfeatures,
    verbose = FALSE
  ) |>
  Seurat::ScaleData(
    vars.to.regress = intersect(PARAM$regress_vars, colnames(obj@meta.data)),
    verbose = FALSE
  ) |>
  Seurat::RunPCA(
    npcs = PARAM$npcs,
    verbose = FALSE
  )

obj <- Seurat::RunUMAP(
  obj,
  reduction = "pca",
  dims = PARAM$dims,
  seed.use = SEEDS$subcluster,
  verbose = FALSE
)

save_plot(
  Seurat::DimPlot(obj, group.by = PARAM$sample_col, pt.size = 0.2),
  file.path(SUB_DIR, "UMAP_raw_by_sample"),
  9,
  7
)

if (PARAM$use_harmony) {
  suppressPackageStartupMessages(library(harmony))

  if (!PARAM$harmony_group %in% colnames(obj@meta.data)) {
    stop("Harmony 分组列不存在：", PARAM$harmony_group)
  }

  set_pipeline_seed(SEEDS$subcluster)
  obj <- RunHarmony(
    obj,
    group.by.vars = PARAM$harmony_group,
    reduction = "pca",
    assay.use = PARAM$assay,
    verbose = FALSE
  )

  obj <- obj |>
    Seurat::FindNeighbors(
      reduction = "harmony",
      dims = PARAM$dims,
      verbose = FALSE
    ) |>
    Seurat::FindClusters(
      resolution = PARAM$resolution,
      random.seed = SEEDS$subcluster,
      verbose = FALSE
    ) |>
    Seurat::RunUMAP(
      reduction = "harmony",
      dims = PARAM$dims,
      seed.use = SEEDS$subcluster,
      verbose = FALSE
    )

  saveRDS(obj, RDS_CLUSTER_HARMONY)
} else {
  obj <- obj |>
    Seurat::FindNeighbors(
      reduction = "pca",
      dims = PARAM$dims,
      verbose = FALSE
    ) |>
    Seurat::FindClusters(
      resolution = PARAM$resolution,
      random.seed = SEEDS$subcluster,
      verbose = FALSE
    ) |>
    Seurat::RunUMAP(
      reduction = "pca",
      dims = PARAM$dims,
      seed.use = SEEDS$subcluster,
      verbose = FALSE
    )

  saveRDS(obj, RDS_CLUSTER_NO_HARMONY)
}

if (length(PARAM$remove_clusters) > 0) {
  obj <- subset(obj, idents = setdiff(levels(Seurat::Idents(obj)), PARAM$remove_clusters))
}

saveRDS(obj, RDS_CLUSTER)

save_plot(
  Seurat::DimPlot(obj, label = TRUE, repel = TRUE),
  file.path(SUB_DIR, "UMAP_clusters"),
  9,
  7
)

save_plot(
  Seurat::DimPlot(obj, group.by = PARAM$sample_col),
  file.path(SUB_DIR, "UMAP_by_sample"),
  9,
  7
)

message("亚群细分完成。先看 ", file.path(SUB_DIR, "UMAP_clusters.pdf"), "，确认 cluster 再跑 marker。")
