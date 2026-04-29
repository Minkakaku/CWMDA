rm(list = ls())

get_script_dir_bootstrap <- function() {
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

script_dir <- get_script_dir_bootstrap()
source(file.path(script_dir, "common_utils.R"), local = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

read_metadata <- function(input_meta, input_rds, meta_sep = "\t", meta_has_rownames = TRUE) {
  if (!is.null(input_rds)) {
    obj <- readRDS(input_rds)
    return(obj@meta.data)
  }

  if (is.null(input_meta)) {
    stop("Please provide either --input-meta or --input-rds.")
  }

  if (meta_has_rownames) {
    read.table(
      input_meta,
      sep = meta_sep,
      header = TRUE,
      row.names = 1,
      check.names = FALSE,
      comment.char = "",
      stringsAsFactors = FALSE
    )
  } else {
    read.table(
      input_meta,
      sep = meta_sep,
      header = TRUE,
      check.names = FALSE,
      comment.char = "",
      stringsAsFactors = FALSE
    )
  }
}

read_order_values <- function(path_value) {
  if (is.null(path_value) || path_value == "") {
    return(NULL)
  }

  tbl <- read.table(
    path_value,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  tbl[[1]]
}

calculate_roe <- function(freq_mat) {
  expected_mat <- outer(rowSums(freq_mat), colSums(freq_mat)) / sum(freq_mat)
  roe_mat <- freq_mat / expected_mat
  roe_mat[!is.finite(roe_mat)] <- NA_real_
  roe_mat
}

args <- parse_cli_args()
input_meta <- cli_string(args, "input-meta")
input_rds <- cli_string(args, "input-rds")
out_dir <- ensure_dir(cli_string(args, "out-dir", "roe_result"))
celltype_col <- cli_string(args, "celltype-col", "cell_type")
group_col <- cli_string(args, "group-col", "group")
meta_sep <- cli_string(args, "meta-sep", "\t")
meta_has_rownames <- cli_bool(args, "meta-has-rownames", TRUE)
group_order_file <- cli_string(args, "group-order-file")
celltype_order_file <- cli_string(args, "celltype-order-file")

meta_df <- read_metadata(
  input_meta = input_meta,
  input_rds = input_rds,
  meta_sep = meta_sep,
  meta_has_rownames = meta_has_rownames
)

if (!celltype_col %in% colnames(meta_df)) {
  stop("celltype column not found: ", celltype_col)
}
if (!group_col %in% colnames(meta_df)) {
  stop("group column not found: ", group_col)
}

meta_df <- meta_df[, c(celltype_col, group_col), drop = FALSE]
colnames(meta_df) <- c("Celltype", "Group")
meta_df <- meta_df[stats::complete.cases(meta_df), , drop = FALSE]

group_levels <- read_order_values(group_order_file)
celltype_levels <- read_order_values(celltype_order_file)

if (!is.null(group_levels)) {
  meta_df$Group <- factor(meta_df$Group, levels = group_levels)
}
if (!is.null(celltype_levels)) {
  meta_df$Celltype <- factor(meta_df$Celltype, levels = celltype_levels)
}

freq_mat <- table(meta_df$Celltype, meta_df$Group, useNA = "no")
roe_mat <- calculate_roe(freq_mat)

write.table(
  as.data.frame.matrix(freq_mat),
  file.path(out_dir, "roe_observed_count.tsv"),
  sep = "\t",
  quote = FALSE
)
write.table(
  as.data.frame.matrix(round(roe_mat, 4)),
  file.path(out_dir, "roe_matrix.tsv"),
  sep = "\t",
  quote = FALSE
)

roe_long <- as.data.frame(as.table(roe_mat), stringsAsFactors = FALSE)
colnames(roe_long) <- c("Celltype", "Group", "RoE")
roe_long <- roe_long %>%
  mutate(
    Class = dplyr::case_when(
      RoE >= 1.05 ~ "Enrichment",
      RoE <= 0.95 ~ "Depletion",
      TRUE ~ "No Change"
    )
  )

write.table(
  roe_long,
  file.path(out_dir, "roe_long.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_bubble <- ggplot(roe_long, aes(x = Group, y = Celltype)) +
  geom_point(aes(size = RoE, color = Class), alpha = 0.9) +
  coord_flip() +
  scale_color_manual(
    values = c(
      "Depletion" = "#4a6fe3",
      "No Change" = "grey70",
      "Enrichment" = "#8e063b"
    )
  ) +
  scale_size(range = c(2, 8)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = group_col, y = celltype_col, size = "Ro/e", color = NULL)

save_plot(p_bubble, file.path(out_dir, "roe_bubble"), width = 10, height = 7)

p_heatmap <- pheatmap::pheatmap(
  roe_mat,
  color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(50),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_color = "black",
  border_color = "#ffffff",
  angle_col = 45,
  cellwidth = max(24, 180 / max(1, ncol(roe_mat))),
  cellheight = max(12, 240 / max(1, nrow(roe_mat)))
)

save_pheatmap(
  p_heatmap,
  file.path(out_dir, "roe_heatmap"),
  width = max(7, ncol(roe_mat) * 1.2),
  height = max(4, nrow(roe_mat) * 0.35)
)

message("ROE finished: ", out_dir)
