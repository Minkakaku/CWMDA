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

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(pheatmap)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

set.seed(123)

read_gmt_fgsea <- function(gmt_file) {
  gmt <- read_table_by_extension(gmt_file, header = FALSE)
  pathways <- split(gmt[, -c(1, 2), drop = FALSE], gmt[[1]])
  lapply(pathways, function(x) unique(stats::na.omit(as.character(unlist(x)))))
}

run_fgsea_one <- function(dep_dt, pathways, gene_col, stat_col, min_size, max_size, nperm) {
  if (!all(c(gene_col, stat_col) %in% colnames(dep_dt))) {
    stop("Dependency table must contain columns: ", gene_col, " and ", stat_col)
  }

  sub_df <- dep_dt[, c(gene_col, stat_col), with = FALSE]
  colnames(sub_df) <- c("gene", "stat")
  sub_df$gene <- as.character(sub_df$gene)
  sub_df$stat <- as.numeric(sub_df$stat)
  sub_df <- sub_df[!is.na(sub_df$gene) & sub_df$gene != "" & !is.na(sub_df$stat), ]
  sub_df <- sub_df[!duplicated(sub_df$gene), ]

  stats_value <- sub_df$stat
  names(stats_value) <- sub_df$gene
  stats_value <- sort(stats_value, decreasing = TRUE)

  fgsea(
    pathways = pathways,
    stats = stats_value,
    minSize = min_size,
    maxSize = max_size,
    nperm = nperm
  )[, .(pathway, NES, ES, pval, padj, size)]
}

calc_drug_score <- function(values, mode = "sum") {
  values <- as.numeric(values)
  values <- values[!is.na(values)]
  if (length(values) == 0) {
    return(NA_real_)
  }
  if (mode == "sum") {
    return(sum(values))
  }
  if (mode == "mean") {
    return(mean(values))
  }
  if (mode == "maxabs") {
    return(values[[which.max(abs(values))]])
  }
  stop("Unknown combine mode: ", mode)
}

scale_0_10 <- function(x) {
  if (length(unique(x)) == 1) {
    return(rep(10, length(x)))
  }
  (x - min(x)) / (max(x) - min(x)) * 10
}

args <- parse_cli_args()
manifest_file <- cli_string(args, "manifest")
gmt_file <- cli_string(args, "gmt-file")
pathways_keep <- cli_csv(args, "pathways")
out_dir <- ensure_dir(cli_string(args, "out-dir", "drug_scoring_outputs"))
gene_col <- cli_string(args, "gene-col", "gene")
stat_col <- cli_string(args, "stat-col", "stat")
min_size <- cli_numeric(args, "min-size", 10)
max_size <- cli_numeric(args, "max-size", 500)
nperm <- cli_numeric(args, "nperm", 10000)
combine_mode <- cli_string(args, "combine-mode", "sum")
top_n <- cli_numeric(args, "top-n", 120)
top_k_label <- cli_numeric(args, "top-k-label", 6)

if (is.null(manifest_file) || is.null(gmt_file)) {
  stop("Please provide --manifest and --gmt-file.")
}

manifest <- read_manifest_table(manifest_file, required_columns = c("drug_name", "dep_file"))
manifest_dir <- dirname(normalizePath(manifest_file, winslash = "/", mustWork = TRUE))
manifest$dep_file <- vapply(
  manifest$dep_file,
  function(path_value) {
    if (grepl("^[A-Za-z]:[/\\\\]|^/|^\\\\\\\\", path_value)) {
      return(path_value)
    }
    normalizePath(file.path(manifest_dir, path_value), winslash = "/", mustWork = FALSE)
  },
  character(1)
)

pathways_all <- read_gmt_fgsea(gmt_file)

if (length(pathways_keep) == 0) {
  pathways_keep <- names(pathways_all)
}

missing_pathways <- setdiff(pathways_keep, names(pathways_all))
if (length(missing_pathways) > 0) {
  stop("These pathways are not found in GMT: ", paste(missing_pathways, collapse = ", "))
}

pathways_use <- pathways_all[pathways_keep]
gsea_list <- vector("list", nrow(manifest))
names(gsea_list) <- manifest$drug_name

for (index_value in seq_len(nrow(manifest))) {
  dep_dt <- data.table::as.data.table(read_table_by_extension(manifest$dep_file[[index_value]], header = TRUE))
  gsea_result <- run_fgsea_one(
    dep_dt = dep_dt,
    pathways = pathways_use,
    gene_col = gene_col,
    stat_col = stat_col,
    min_size = min_size,
    max_size = max_size,
    nperm = nperm
  )
  gsea_result <- gsea_result[match(names(pathways_use), gsea_result$pathway), ]
  gsea_list[[index_value]] <- gsea_result
}

score_long <- do.call(
  rbind,
  lapply(names(gsea_list), function(drug_name) {
    result_df <- gsea_list[[drug_name]]
    data.frame(
      Drug = drug_name,
      Pathway = result_df$pathway,
      NES = result_df$NES,
      ES = result_df$ES,
      pval = result_df$pval,
      padj = result_df$padj,
      size = result_df$size,
      stringsAsFactors = FALSE
    )
  })
)

write.table(
  score_long,
  file.path(out_dir, "drug_fgsea_long.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

score_df <- score_long %>%
  select(Drug, Pathway, NES) %>%
  tidyr::pivot_wider(names_from = Pathway, values_from = NES) %>%
  as.data.frame()

rownames(score_df) <- score_df$Drug
score_df$Drug <- NULL

pdf(file.path(out_dir, "program_scores_heatmap.pdf"), width = 4.5, height = 8)
pheatmap(
  mat = as.matrix(score_df),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  border_color = NA,
  fontsize_row = 7,
  fontsize_col = 9,
  main = "Program scores (NES)"
)
dev.off()

write.table(
  cbind(Drug = rownames(score_df), score_df),
  file.path(out_dir, "program_scores_matrix.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

drug_score <- apply(score_df, 1, calc_drug_score, mode = combine_mode)
drug_names <- rownames(score_df)

pair_df <- expand.grid(
  A = seq_along(drug_names),
  B = seq_along(drug_names)
) |>
  dplyr::filter(A < B) |>
  dplyr::mutate(
    DrugA = drug_names[A],
    DrugB = drug_names[B],
    scoreA = drug_score[DrugA],
    scoreB = drug_score[DrugB],
    synergy_raw = scoreA + scoreB
  ) |>
  dplyr::filter(!is.na(synergy_raw))

pair_top <- pair_df |>
  dplyr::arrange(dplyr::desc(synergy_raw)) |>
  dplyr::slice_head(n = top_n)

pair_top$synergy_0_10 <- scale_0_10(pair_top$synergy_raw)

top_pairs_txt <- pair_top |>
  dplyr::slice_head(n = top_k_label) |>
  dplyr::mutate(
    line = paste0(
      dplyr::row_number(),
      ". ",
      DrugA,
      " + ",
      DrugB,
      " = ",
      sprintf("%.2f", synergy_0_10)
    )
  ) |>
  dplyr::pull(line) |>
  paste(collapse = "\n")

p_pair <- ggplot(pair_top, aes(x = A - 1, y = B - 1)) +
  geom_point(aes(size = synergy_0_10, color = synergy_0_10), alpha = 0.9) +
  scale_size_continuous(range = c(2, 7)) +
  labs(
    title = "Top synergistic pairs as a sparse interaction map",
    x = "Drug index (A)",
    y = "Drug index (B)",
    color = "Integrated synergy score",
    size = "Integrated synergy score"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 12),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  ) +
  annotate(
    geom = "label",
    x = max(pair_top$A - 1) * 0.75,
    y = min(pair_top$B - 1) + 2,
    label = paste0("Top pairs\n(rank: A + B = score)\n\n", top_pairs_txt),
    hjust = 0,
    size = 3
  )

save_plot(p_pair, file.path(out_dir, "top_pairs_sparse_map"), width = 12, height = 5)

write.table(
  pair_top,
  file.path(out_dir, "top_pairs_table.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Drug score FGSEA finished: ", out_dir)
