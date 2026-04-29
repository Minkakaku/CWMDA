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
  library(ggsci)
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

read_order_table <- function(path_value) {
  if (is.null(path_value) || path_value == "") {
    return(NULL)
  }

  read.table(
    path_value,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

build_palette <- function(n_value) {
  if (n_value <= 10) {
    pal_npg()(n_value)
  } else {
    colorRampPalette(pal_npg()(10))(n_value)
  }
}

args <- parse_cli_args()
input_meta <- cli_string(args, "input-meta")
input_rds <- cli_string(args, "input-rds")
out_dir <- ensure_dir(cli_string(args, "out-dir", "proportion_result"))
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

group_order <- read_order_table(group_order_file)
celltype_order <- read_order_table(celltype_order_file)
color_manual <- NULL

if (!is.null(group_order) && ncol(group_order) >= 1) {
  meta_df$Group <- factor(meta_df$Group, levels = group_order[[1]])
}

if (!is.null(celltype_order) && ncol(celltype_order) >= 1) {
  meta_df$Celltype <- factor(meta_df$Celltype, levels = celltype_order[[1]])
  if (ncol(celltype_order) >= 2) {
    color_manual <- stats::setNames(celltype_order[[2]], celltype_order[[1]])
  }
}

count_df <- meta_df %>%
  count(Celltype, Group, name = "Count") %>%
  group_by(Group) %>%
  mutate(GroupProportion = Count / sum(Count)) %>%
  ungroup() %>%
  group_by(Celltype) %>%
  mutate(CelltypeProportion = Count / sum(Count)) %>%
  ungroup() %>%
  mutate(OverallProportion = Count / sum(Count))

write.table(
  count_df,
  file.path(out_dir, "cell_proportion.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

celltype_levels <- unique(as.character(count_df$Celltype))
celltype_levels <- celltype_levels[!is.na(celltype_levels)]

if (is.null(color_manual)) {
  color_manual <- stats::setNames(
    rev(build_palette(length(celltype_levels))),
    rev(celltype_levels)
  )
}

base_theme <- theme_classic() +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.4),
    axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
    legend.title = element_blank()
  )

p_group <- ggplot(count_df, aes(x = Group, y = GroupProportion, fill = Celltype)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = color_manual) +
  labs(title = "Cell proportion by group", x = group_col, y = "Proportion") +
  base_theme

p_celltype <- ggplot(count_df, aes(x = Celltype, y = CelltypeProportion, fill = Group)) +
  geom_col(position = "stack") +
  coord_flip() +
  labs(title = "Group composition within each cell type", x = celltype_col, y = "Proportion") +
  theme_classic() +
  theme(legend.title = element_blank())

p_pie <- ggplot(count_df, aes(x = Group, y = GroupProportion, fill = Celltype)) +
  geom_col(width = 1, color = "white", linewidth = 0.2) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = color_manual) +
  labs(title = "Cell proportion pie view", x = NULL, y = NULL) +
  theme_void() +
  theme(legend.title = element_blank())

save_plot(p_group, file.path(out_dir, "cell_proportion_group"), width = 8, height = 6)
save_plot(p_celltype, file.path(out_dir, "cell_proportion_celltype"), width = 9, height = 7)
save_plot(p_pie, file.path(out_dir, "cell_proportion_pie"), width = 8, height = 6)

message("Cell proportion finished: ", out_dir)
