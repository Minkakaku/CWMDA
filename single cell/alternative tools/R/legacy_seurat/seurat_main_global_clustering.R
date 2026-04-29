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

set_pipeline_seed(PIPELINE_SEEDS$clustering)
options(future.globals.maxSize = 20 * 1024^3)

sc <- read_step_rds(RDS_QC)

sc <- Seurat::SCTransform(
  sc,
  vars.to.regress = "percent.mt",
  verbose = FALSE
)

sc <- Seurat::RunPCA(sc, npcs = 50, verbose = FALSE)

set_pipeline_seed(PIPELINE_SEEDS$clustering)
sc <- RunHarmony(
  sc,
  group.by.vars = "orig.ident",
  reduction.use = "pca",
  assay.use = "SCT"
)

p_elbow <- Seurat::ElbowPlot(sc, ndims = 50)
save_plot(p_elbow, file.path(GLOBAL_DIR, "PCA_elbow"))

sc <- Seurat::FindNeighbors(sc, reduction = "harmony", dims = 1:25)
sc <- Seurat::FindClusters(
  sc,
  resolution = 0.35,
  random.seed = PIPELINE_SEEDS$clustering
)
sc <- Seurat::RunUMAP(
  sc,
  reduction = "harmony",
  dims = 1:25,
  seed.use = PIPELINE_SEEDS$clustering
)

sc <- assign_group_metadata(sc)

print(table(sc$samples))
print(table(sc$orig.ident, sc$group))

p_umap_cluster <- Seurat::DimPlot(
  sc,
  label = TRUE,
  repel = TRUE
)
save_plot(p_umap_cluster, file.path(GLOBAL_DIR, "UMAP_clusters"))

p_umap_sample <- Seurat::DimPlot(
  sc,
  group.by = "samples"
)
save_plot(p_umap_sample, file.path(GLOBAL_DIR, "UMAP_by_sample"))

p_feat_contamination <- Seurat::FeaturePlot(
  sc,
  features = "Contamination"
)
save_plot(
  p_feat_contamination,
  file.path(GLOBAL_DIR, "Feature_Contamination"),
  width = 5,
  height = 5
)

save_step_rds(sc, RDS_CLUSTER)

message("全局聚类完成，cluster 数：", length(unique(sc$seurat_clusters)))
