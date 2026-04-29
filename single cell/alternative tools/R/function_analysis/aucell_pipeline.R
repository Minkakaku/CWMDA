rm(list = ls())
set.seed(1234)

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
  library(AUCell)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

load_gene_sets <- function(gene_set_file = NULL, preset = "cd8_tpex_tex_terminal") {
  if (!is.null(gene_set_file) && gene_set_file != "") {
    if (grepl("\\.gmt$", gene_set_file, ignore.case = TRUE)) {
      lines_value <- readLines(gene_set_file, warn = FALSE)
      gene_sets <- lapply(lines_value, function(line_value) {
        fields <- strsplit(line_value, "\t", fixed = TRUE)[[1]]
        set_names <- fields[[1]]
        genes <- fields[seq.int(3, length(fields))]
        list(name = set_names, genes = genes)
      })
      return(stats::setNames(lapply(gene_sets, `[[`, "genes"), vapply(gene_sets, `[[`, character(1), "name")))
    }

    gene_df <- read.table(
      gene_set_file,
      sep = "\t",
      header = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    if (!all(c("set", "gene") %in% colnames(gene_df))) {
      stop("Gene set file must contain columns named 'set' and 'gene', or be a GMT file.")
    }
    return(split(gene_df$gene, gene_df$set))
  }

  if (preset != "cd8_tpex_tex_terminal") {
    stop("Unsupported preset: ", preset)
  }

  list(
    Tpex = c(
      "Tcf7", "Slamf6", "Pdcd1", "Tox", "Ikzf2", "Myb",
      "Bach2", "Lef1", "Il7r", "Cxcr5", "Nr4a2", "Xcl1", "Sell", "Ltb"
    ),
    TexTerminal = c(
      "Pdcd1", "Havcr2", "Lag3", "Tigit", "Ctla4", "Entpd1",
      "Tox", "Tox2", "Prdm1", "Batf", "Eomes",
      "Cxcl13", "Nr4a1", "Nr4a2", "Nr4a3"
    )
  )
}

match_gene_sets <- function(gene_sets, genes_in_data) {
  gene_map <- genes_in_data
  names(gene_map) <- toupper(genes_in_data)
  lapply(gene_sets, function(gene_list) {
    matched <- gene_map[toupper(gene_list)]
    unique(unname(matched[!is.na(matched)]))
  })
}

args <- parse_cli_args()
input_rds <- cli_string(args, "input-rds")
out_dir <- ensure_dir(cli_string(args, "out-dir", "aucell_result"))
gene_set_file <- cli_string(args, "gene-set-file")
preset <- cli_string(args, "preset", "cd8_tpex_tex_terminal")
assay_name <- cli_string(args, "assay")
reduction_name <- cli_string(args, "reduction", "umap")
ident_col <- cli_string(args, "ident-col")
split_col <- cli_string(args, "split-col")
auc_fraction <- cli_numeric(args, "auc-max-rank-fraction", 0.05)
ncores <- cli_numeric(args, "ncores", 1)

if (is.null(input_rds)) {
  stop("Please provide --input-rds.")
}

seu <- readRDS(input_rds)
if (!is.null(assay_name) && assay_name %in% names(seu@assays)) {
  DefaultAssay(seu) <- assay_name
}

if (!is.null(ident_col) && ident_col %in% colnames(seu@meta.data)) {
  Idents(seu) <- seu[[ident_col, drop = TRUE]]
}

expr_mat <- tryCatch(
  GetAssayData(seu, assay = DefaultAssay(seu), layer = "data"),
  error = function(e) GetAssayData(seu, assay = DefaultAssay(seu), slot = "data")
)

if (!inherits(expr_mat, "dgCMatrix")) {
  expr_mat <- as(as.matrix(expr_mat), "dgCMatrix")
}

gene_sets <- load_gene_sets(gene_set_file = gene_set_file, preset = preset)
gene_sets <- match_gene_sets(gene_sets, rownames(expr_mat))

gene_set_sizes <- vapply(gene_sets, length, numeric(1))
write.table(
  data.frame(
    set = rep(names(gene_sets), gene_set_sizes),
    gene = unlist(gene_sets, use.names = FALSE)
  ),
  file.path(out_dir, "gene_sets_used.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

rankings <- AUCell_buildRankings(expr_mat, nCores = ncores, plotStats = FALSE)
cells_auc <- AUCell_calcAUC(
  gene_sets,
  rankings,
  aucMaxRank = ceiling(auc_fraction * nrow(expr_mat))
)

auc_df <- as.data.frame(t(getAUC(cells_auc)))
colnames(auc_df) <- paste0("AUC_", colnames(auc_df))
auc_df <- auc_df[colnames(seu), , drop = FALSE]
seu <- AddMetaData(seu, auc_df)

write.table(
  data.frame(Cell = rownames(auc_df), auc_df, check.names = FALSE),
  file.path(out_dir, "auc_scores.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

for (auc_col in colnames(auc_df)) {
  feature_plot <- FeaturePlot(
    object = seu,
    features = auc_col,
    reduction = reduction_name,
    split.by = split_col,
    cols = c("lightgrey", "#B83D3D")
  )

  if (!is.null(ident_col) && ident_col %in% colnames(seu@meta.data)) {
    violin_plot <- VlnPlot(
      object = seu,
      features = auc_col,
      group.by = ident_col,
      split.by = split_col,
      pt.size = 0
    )
    combined_plot <- feature_plot / violin_plot
    save_plot(combined_plot, file.path(out_dir, auc_col), width = 10, height = 10)
  } else {
    save_plot(feature_plot, file.path(out_dir, auc_col), width = 8, height = 6)
  }
}

saveRDS(seu, file.path(out_dir, "aucell_result.rds"))
message("AUCell finished: ", out_dir)
