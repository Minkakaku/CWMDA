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

set_pipeline_seed(PIPELINE_SEEDS$proportion)

sc <- read_step_rds(RDS_ANNOTATED)
sc <- assign_group_metadata(sc)

if (!"cell_type" %in% colnames(sc@meta.data)) {
  stop("当前对象没有 cell_type 列，请先运行 04_cell_annotation.R。")
}

if (anyNA(sc$cell_type)) {
  stop("cell_type 中仍有 NA，请先补全注释再计算比例。")
}

if (anyNA(sc$group)) {
  stop("group 中存在 NA，请先检查样本命名是否满足 Con/Exp 规则。")
}

plot_celltype_distribution(
  data = sc,
  sample_col = "samples",
  celltype_col = "cell_type",
  palette_name = "Paired",
  output_dir = PROPORTION_DIR,
  output_file = "sc_proportion_samples.pdf",
  width = 6,
  height = 6,
  text_size = 7,
  x_angle = 45,
  legend_position = "bottom"
)

plot_celltype_distribution(
  data = sc,
  sample_col = "group",
  celltype_col = "cell_type",
  palette_name = "Paired",
  output_dir = PROPORTION_DIR,
  output_file = "sc_proportion_group.pdf",
  width = 6,
  height = 6,
  text_size = 7,
  x_angle = 45,
  legend_position = "bottom"
)

count_mat <- table(sc$orig.ident, sc$cell_type)
count_df <- as.data.frame.matrix(count_mat)

write.csv(
  count_df,
  file = file.path(PROPORTION_DIR, "sample_by_celltype_count.csv"),
  quote = FALSE
)

roe_mat <- Startrac::calTissueDist(
  sc@meta.data,
  byPatient = FALSE,
  colname.patient = "orig.ident",
  colname.cluster = "cell_type",
  colname.tissue = "group",
  method = "chisq",
  min.rowSum = 0
)

saveRDS(roe_mat, RDS_ROE)
write.csv(
  as.data.frame(roe_mat),
  file = file.path(PROPORTION_DIR, "ROE_matrix.csv"),
  quote = FALSE
)

p_roe <- pheatmap::pheatmap(
  roe_mat,
  scale = "none",
  display_numbers = TRUE,
  number_color = "black",
  number_format = "%.2f",
  color = grDevices::colorRampPalette(
    c("#2166AC", "#b55c67ff", "#B2182B")
  )(100),
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
  filename = file.path(PROPORTION_DIR, "ROE"),
  width = 8,
  height = 6
)

save_step_rds(sc, RDS_FINAL)

message("比例统计与 Ro/e 计算完成。")
