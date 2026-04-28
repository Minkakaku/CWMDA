rm(list = ls())

script_dir <- {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
  } else if (!is.null(sys.frames()[[1]]$ofile)) {
    dirname(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE))
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }
}

source(file.path(script_dir, "common_utils.R"), local = FALSE)
source(file.path(script_dir, "cibersort_core.R"), local = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

read_gene_matrix <- function(file_path) {
  df <- read_table_by_extension(file_path, header = TRUE)
  if (ncol(df) < 2) {
    stop("Matrix must contain at least one gene column and one sample column: ", file_path)
  }
  rownames(df) <- as.character(df[[1]])
  df[[1]] <- NULL
  df <- df[!duplicated(rownames(df)) & !is.na(rownames(df)) & rownames(df) != "", , drop = FALSE]
  as.data.frame(df, check.names = FALSE)
}

read_marker_list <- function(file_path) {
  df <- read_table_by_extension(file_path, header = TRUE)
  unique(as.character(df[[1]]))
}

args <- parse_cli_args()
sig_matrix_file <- cli_string(args, "sig-matrix")
mixture_file <- cli_string(args, "mixture-file")
mark_matrix_file <- cli_string(args, "mark-matrix")
out_dir <- ensure_dir(cli_string(args, "out-dir", "cibersort_result"))
perm <- cli_numeric(args, "perm", 100)
qn <- cli_bool(args, "qn", FALSE)
min_nonzero_fraction <- cli_numeric(args, "min-nonzero-fraction", 0.5)

if (is.null(sig_matrix_file) || is.null(mixture_file)) {
  stop("Please provide --sig-matrix and --mixture-file.")
}

signature_df <- read_gene_matrix(sig_matrix_file)
if (!is.null(mark_matrix_file)) {
  marker_genes <- read_marker_list(mark_matrix_file)
  keep_genes <- intersect(rownames(signature_df), marker_genes)
  signature_df <- signature_df[keep_genes, , drop = FALSE]
}

mixture_df <- read_gene_matrix(mixture_file)

raw_output_file <- file.path(out_dir, "CIBERSORT-Results.txt")
result_matrix <- run_cibersort_matrix(
  signature_matrix = signature_df,
  mixture_matrix = mixture_df,
  perm = perm,
  qn = qn,
  output_file = raw_output_file
)

save(result_matrix, file = file.path(out_dir, "cibersort_output_obj.Rdata"))

result_df <- as.data.frame(result_matrix, check.names = FALSE)
result_df <- tibble::rownames_to_column(result_df, "Mixture")

write.table(
  result_df,
  file.path(out_dir, "cibersort_result_full.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

fraction_df <- result_df %>%
  select(-`P-value`, -Correlation, -RMSE)

rownames(fraction_df) <- fraction_df$Mixture
fraction_df$Mixture <- NULL

keep_cols <- apply(
  fraction_df,
  2,
  function(x) mean(as.numeric(x) == 0, na.rm = TRUE) < min_nonzero_fraction
)

filtered_df <- fraction_df[, keep_cols, drop = FALSE]
filtered_df <- cbind(Mixture = rownames(filtered_df), filtered_df)

write.table(
  filtered_df,
  file.path(out_dir, "cibersort_result_matrix.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("CIBERSORT finished: ", out_dir)
