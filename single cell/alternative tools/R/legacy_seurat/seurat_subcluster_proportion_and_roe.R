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

set_pipeline_seed(SEEDS$proportion)

obj <- read_required_rds(RDS_ANNOTATED)
obj <- prepare_subcluster_input(obj)

if (!"sub_cell_type" %in% colnames(obj@meta.data)) {
  stop("当前对象缺少 sub_cell_type，请先运行 03_subcluster_annotation.R。")
}

if (anyNA(obj$sub_cell_type)) {
  stop("sub_cell_type 中仍有 NA，请先补全注释。")
}

if (!"group" %in% colnames(obj@meta.data) || anyNA(obj$group)) {
  stop("group 缺失或存在 NA，请先检查样本分组。")
}

plot_celltype_distribution(
  data = obj,
  sample_col = "samples",
  celltype_col = "sub_cell_type",
  palette_name = "Paired",
  output_dir = SUB_DIR,
  output_file = paste0(SUBCLUSTER_NAME, "_subtype_by_sample.pdf"),
  width = 6,
  height = 6,
  text_size = 7,
  x_angle = 45,
  legend_position = "bottom"
)

plot_celltype_distribution(
  data = obj,
  sample_col = "group",
  celltype_col = "sub_cell_type",
  palette_name = "Paired",
  output_dir = SUB_DIR,
  output_file = paste0(SUBCLUSTER_NAME, "_subtype_by_group.pdf"),
  width = 6,
  height = 6,
  text_size = 7,
  x_angle = 45,
  legend_position = "bottom"
)

count_mat <- table(obj$orig.ident, obj$sub_cell_type)
count_df <- as.data.frame.matrix(count_mat)

write.csv(
  count_df,
  file.path(SUB_DIR, "sample_by_sub_celltype_count.csv"),
  quote = FALSE
)

roe_mat <- Startrac::calTissueDist(
  obj@meta.data,
  byPatient = FALSE,
  colname.patient = "orig.ident",
  colname.cluster = "sub_cell_type",
  colname.tissue = "group",
  method = "chisq",
  min.rowSum = 0
)

saveRDS(roe_mat, RDS_ROE)
write.csv(
  as.data.frame(roe_mat),
  file.path(SUB_DIR, "ROE_matrix.csv"),
  quote = FALSE
)

p_roe <- pheatmap::pheatmap(
  roe_mat,
  scale = "none",
  display_numbers = TRUE,
  number_color = "black",
  number_format = "%.2f",
  color = grDevices::colorRampPalette(c("#2166AC", "#b55c67ff", "#B2182B"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = "grey80",
  cellwidth = 30,
  cellheight = 25,
  silent = TRUE
)

save_pheatmap(
  p_roe,
  filename = file.path(SUB_DIR, "ROE"),
  width = 8,
  height = 6
)

message("亚群比例和 Ro/e 已计算完成。")
