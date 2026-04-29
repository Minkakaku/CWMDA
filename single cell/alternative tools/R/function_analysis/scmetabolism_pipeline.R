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
  library(scMetabolism)
  library(Seurat)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
})

rename_genes_seurat <- function(obj, new_names, gene_use = NULL, assay_name = "RNA") {
  assay_list <- names(obj@assays)
  if (!is.null(gene_use)) {
    obj <- subset(obj, features = gene_use)
  }

  order_name <- function(source_names, target_names, ref_names) {
    target_names <- make.names(target_names, unique = TRUE)
    name_df <- data.frame(source = source_names, target = target_names, stringsAsFactors = FALSE)
    rownames(name_df) <- name_df$source
    name_df[ref_names, , drop = FALSE]
  }

  current_genes <- rownames(obj[[assay_name]])
  name_df <- order_name(source_names = current_genes, target_names = new_names, ref_names = current_genes)

  for (assay_item in assay_list) {
    assay_obj <- obj@assays[[assay_item]]
    feature_names <- rownames(assay_obj)
    assay_df <- order_name(
      source_names = name_df$source,
      target_names = name_df$target,
      ref_names = feature_names
    )

    if (length(assay_obj@counts)) {
      assay_obj@counts@Dimnames[[1]] <- assay_df$target
    }
    if (length(assay_obj@data)) {
      assay_obj@data@Dimnames[[1]] <- assay_df$target
    }
    if (length(assay_obj@scale.data)) {
      rownames(assay_obj@scale.data) <- assay_df$target
    }
    if (length(assay_obj@var.features)) {
      vf_df <- order_name(
        source_names = assay_df$source,
        target_names = assay_df$target,
        ref_names = assay_obj@var.features
      )
      assay_obj@var.features <- vf_df$target
    }

    obj@assays[[assay_item]] <- assay_obj
  }

  if (length(obj@reductions$pca)) {
    pca_loadings <- obj@reductions$pca@feature.loadings
    loading_df <- order_name(
      source_names = name_df$source,
      target_names = name_df$target,
      ref_names = rownames(pca_loadings)
    )
    rownames(obj@reductions$pca@feature.loadings) <- loading_df$target
  }

  try(obj[[assay_name]]@meta.features <- data.frame(row.names = rownames(obj[[assay_name]])), silent = TRUE)
  obj
}

convert_mouse_object_to_human <- function(obj, assay_name = "RNA") {
  if (!requireNamespace("nichenetr", quietly = TRUE)) {
    stop("Mouse conversion requires the nichenetr package.")
  }

  mouse_genes <- rownames(obj[[assay_name]])
  human_genes <- nichenetr::convert_mouse_to_human_symbols(mouse_genes)
  gene_map <- data.frame(
    human = human_genes,
    mouse = mouse_genes,
    stringsAsFactors = FALSE
  )
  gene_map <- stats::na.omit(gene_map)
  gene_map <- gene_map[!duplicated(gene_map$human), , drop = FALSE]
  gene_map <- gene_map[!duplicated(gene_map$mouse), , drop = FALSE]

  converted <- subset(obj, features = gene_map$mouse)
  rename_genes_seurat(
    converted,
    new_names = gene_map$human,
    gene_use = gene_map$mouse,
    assay_name = assay_name
  )
}

args <- parse_cli_args()
input_rds <- cli_string(args, "input-rds")
out_dir <- ensure_dir(cli_string(args, "out-dir", "scmetabolism_result"))
assay_name <- cli_string(args, "assay", "RNA")
species <- cli_string(args, "species", "human")
metabolism_type <- cli_string(args, "metabolism-type", "KEGG")
group_col <- cli_string(args, "group-col")
top_pathways <- cli_numeric(args, "top-pathways", 30)
ncores <- cli_numeric(args, "ncores", 1)

if (is.null(input_rds)) {
  stop("Please provide --input-rds.")
}

species <- match.arg(species, c("human", "mouse"))
seu <- readRDS(input_rds)
DefaultAssay(seu) <- assay_name

analysis_obj <- if (species == "mouse") {
  convert_mouse_object_to_human(seu, assay_name = assay_name)
} else {
  seu
}

metabolism_obj <- sc.metabolism.Seurat(
  obj = analysis_obj,
  method = "AUCell",
  imputation = FALSE,
  ncores = ncores,
  metabolism.type = metabolism_type
)

score_mat <- metabolism_obj@assays[["METABOLISM"]][["score"]]
write.table(
  as.data.frame(score_mat),
  file.path(out_dir, "metabolism_score_matrix.tsv"),
  sep = "\t",
  quote = FALSE
)

phenotype <- if (!is.null(group_col) && group_col %in% colnames(metabolism_obj@meta.data)) {
  as.character(metabolism_obj@meta.data[[group_col]])
} else {
  as.character(Idents(metabolism_obj))
}

metabolism_obj$analysis_group <- phenotype

score_df <- as.data.frame(t(score_mat), check.names = FALSE)
score_df$analysis_group <- phenotype

avg_df <- aggregate(score_df[, setdiff(colnames(score_df), "analysis_group"), drop = FALSE], list(score_df$analysis_group), mean)
rownames(avg_df) <- avg_df$Group.1
avg_df <- t(as.matrix(avg_df[, -1, drop = FALSE]))
avg_df <- avg_df[order(apply(avg_df, 1, stats::var), decreasing = TRUE), , drop = FALSE]

write.table(
  as.data.frame(avg_df),
  file.path(out_dir, "metabolism_group_mean.tsv"),
  sep = "\t",
  quote = FALSE
)

top_n <- min(top_pathways, nrow(avg_df))
heatmap_df <- avg_df[seq_len(top_n), , drop = FALSE]

p_heatmap <- pheatmap(
  heatmap_df,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#1A5592", "white", "#B83D3D"))(100),
  show_colnames = TRUE
)

save_pheatmap(
  p_heatmap,
  file.path(out_dir, "metabolism_heatmap"),
  width = max(7, ncol(heatmap_df) * 1.4),
  height = max(6, nrow(heatmap_df) * 0.25)
)

dot_plot <- DotPlot.metabolism(
  obj = metabolism_obj,
  pathway = rownames(score_mat),
  phenotype = "analysis_group",
  norm = "y"
)

save_plot(dot_plot, file.path(out_dir, "metabolism_dotplot"), width = 11, height = 12)
saveRDS(metabolism_obj, file.path(out_dir, "scmetabolism_result.rds"))

message("scMetabolism finished: ", out_dir)
