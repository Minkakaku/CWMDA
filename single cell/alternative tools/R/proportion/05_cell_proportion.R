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
source(file.path(script_dir, "00_pipeline_common.R"), local = FALSE)

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
  stop("group 中存在 NA，请先检查样本命名是否满足 assign_group_metadata() 中的分组规则。")
}

# =========================================================
# 在这里设置多个比较组别
# 每个元素代表一次单独比较，并分别输出结果
# =========================================================
comparison_list <- list(
  c("Normal", "Microinvasive"),
  c("Normal", "Invasive"),
  c("Microinvasive", "Invasive"),
  # c("Normal", "Microinvasive_Para"),
  # c("Normal", "Invasive_Para"),
  c("Microinvasive_Para", "Microinvasive"),
  c("Invasive_Para", "Invasive")
)

all_groups <- unique(sc$group)

for (selected_groups in comparison_list) {

  if (length(selected_groups) < 2) {
    stop("comparison_list 中每个比较至少需要包含 2 个组别。")
  }

  missing_groups <- setdiff(selected_groups, all_groups)
  if (length(missing_groups) > 0) {
    stop(
      "以下组别在 sc$group 中不存在：",
      paste(missing_groups, collapse = ", "),
      "\n当前可用组别为：",
      paste(sort(all_groups), collapse = ", ")
    )
  }

  sc_sub <- subset(sc, subset = group %in% selected_groups)
  sc_sub$group <- factor(sc_sub$group, levels = selected_groups)

  sample_levels <- unique(sc_sub$samples)
  sc_sub$samples <- factor(sc_sub$samples, levels = sample_levels)

  group_tag <- paste(selected_groups, collapse = "_vs_")
  group_tag <- gsub("[^A-Za-z0-9_]+", "_", group_tag)

  # =========================
  # 1. 按样本画比例图
  # =========================
  plot_celltype_distribution(
    data = sc_sub,
    sample_col = "samples",
    celltype_col = "cell_type",
    palette_name = "Paired",
    output_dir = PROPORTION_DIR,
    output_file = paste0("sc_proportion_samples_", group_tag, ".pdf"),
    width = 6,
    height = 6,
    text_size = 7,
    x_angle = 45,
    legend_position = "bottom"
  )

  # =========================
  # 2. 按组别画比例图
  # =========================
  plot_celltype_distribution(
    data = sc_sub,
    sample_col = "group",
    celltype_col = "cell_type",
    palette_name = "Paired",
    output_dir = PROPORTION_DIR,
    output_file = paste0("sc_proportion_group_", group_tag, ".pdf"),
    width = 6,
    height = 6,
    text_size = 7,
    x_angle = 45,
    legend_position = "bottom"
  )

  # =========================
  # 3. 导出样本 × cell type 计数表
  # =========================
  count_mat <- table(sc_sub$orig.ident, sc_sub$cell_type)
  count_df <- as.data.frame.matrix(count_mat)

  write.csv(
    count_df,
    file = file.path(PROPORTION_DIR, paste0("sample_by_celltype_count_", group_tag, ".csv")),
    quote = FALSE
  )

  message("已完成比较：", paste(selected_groups, collapse = " vs "))
}

save_step_rds(sc, RDS_FINAL)

message("所有比例统计计算完成。")