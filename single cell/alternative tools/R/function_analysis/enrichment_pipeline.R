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
  library(clusterProfiler)
})

resolve_org_db <- function(organism) {
  pkg_name <- switch(
    organism,
    human = "org.Hs.eg.db",
    mouse = "org.Mm.eg.db",
    stop("Unsupported organism: ", organism)
  )

  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    stop("Package not installed: ", pkg_name)
  }

  get(pkg_name, envir = asNamespace(pkg_name))
}

convert_deg <- function(deg_df, org_db) {
  gene_df <- clusterProfiler::bitr(
    deg_df$gene,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org_db
  )
  dplyr::inner_join(deg_df, gene_df, by = c("gene" = "SYMBOL"))
}

save_enrichment_plot <- function(enr_df, title_text, output_prefix) {
  plot_df <- enr_df[seq_len(min(15, nrow(enr_df))), , drop = FALSE]

  p <- ggplot(plot_df, aes(x = Count, y = reorder(Description, Count), color = -log10(p.adjust))) +
    geom_point(aes(size = Count)) +
    scale_color_gradient(low = "#059CFE", high = "#FE05EB") +
    labs(
      title = title_text,
      x = "Gene count",
      y = NULL,
      color = "-log10(adj.P)",
      size = "Gene count"
    ) +
    theme_classic() +
    theme(legend.position = "bottom")

  save_plot(p, output_prefix, width = 10, height = 8)
}

perform_enrichment_analysis <- function(
  deg_df,
  out_dir,
  label,
  enrichment_type,
  organism,
  org_db
) {
  if (nrow(deg_df) < 10) {
    return(NULL)
  }

  deg_map <- convert_deg(deg_df, org_db)
  if (nrow(deg_map) < 10) {
    return(NULL)
  }

  enrichment_obj <- if (enrichment_type == "GO") {
    enrichGO(
      gene = deg_map$ENTREZID,
      OrgDb = org_db,
      ont = "BP",
      readable = TRUE
    )
  } else {
    enrichKEGG(
      gene = deg_map$ENTREZID,
      organism = if (organism == "human") "hsa" else "mmu"
    )
  }

  enr_df <- as.data.frame(enrichment_obj)
  if (nrow(enr_df) == 0) {
    return(NULL)
  }

  write.table(
    enr_df,
    file.path(out_dir, paste0(tolower(enrichment_type), "_result.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  save_enrichment_plot(
    enr_df = enr_df,
    title_text = label,
    output_prefix = file.path(out_dir, tolower(enrichment_type))
  )
}

parse_comparisons <- function(entries) {
  if (length(entries) == 0) {
    stop("Please provide --comparisons in test:control format, e.g. T1:CON,T2:CON.")
  }

  lapply(entries, function(entry) {
    parts <- trimws(unlist(strsplit(entry, ":", fixed = TRUE)))
    if (length(parts) != 2) {
      stop("Invalid comparison: ", entry, ". Use test:control.")
    }
    list(test = parts[[1]], control = parts[[2]])
  })
}

args <- parse_cli_args()
input_rds <- cli_string(args, "input-rds")
out_dir <- ensure_dir(cli_string(args, "out-dir", "enrichment_result"))
celltype_col <- cli_string(args, "celltype-col", "cell_type")
targets <- cli_csv(args, "targets")
group_col <- cli_string(args, "group-col", "group")
comparisons <- parse_comparisons(cli_csv(args, "comparisons"))
organism <- cli_string(args, "organism", "mouse")
logfc_threshold <- cli_numeric(args, "logfc-threshold", 0.25)
min_pct <- cli_numeric(args, "min-pct", 0.1)
p_adj_threshold <- cli_numeric(args, "p-adj-threshold", 0.05)
only_positive <- cli_bool(args, "only-positive", TRUE)

if (is.null(input_rds)) {
  stop("Please provide --input-rds.")
}

seu <- readRDS(input_rds)
org_db <- resolve_org_db(organism)

if (!celltype_col %in% colnames(seu@meta.data)) {
  stop("celltype column not found: ", celltype_col)
}
if (!group_col %in% colnames(seu@meta.data)) {
  stop("group column not found: ", group_col)
}

Idents(seu) <- seu[[celltype_col, drop = TRUE]]
if (length(targets) == 0) {
  targets <- levels(Idents(seu))
}

for (celltype_name in targets) {
  message("Cell type: ", celltype_name)

  for (comparison in comparisons) {
    test_name <- comparison$test
    control_name <- comparison$control
    comparison_name <- paste0(test_name, "_vs_", control_name)
    cmp_out_dir <- ensure_dir(file.path(out_dir, make.names(celltype_name), comparison_name))

    deg_df <- FindMarkers(
      object = seu,
      subset.ident = celltype_name,
      group.by = group_col,
      ident.1 = test_name,
      ident.2 = control_name,
      only.pos = only_positive,
      logfc.threshold = logfc_threshold,
      min.pct = min_pct
    ) %>%
      tibble::rownames_to_column("gene")

    lfc_col <- if ("avg_log2FC" %in% colnames(deg_df)) {
      "avg_log2FC"
    } else if ("avg_logFC" %in% colnames(deg_df)) {
      "avg_logFC"
    } else {
      stop("Cannot find logFC column in FindMarkers result.")
    }

    deg_df <- deg_df %>%
      filter(p_val_adj <= p_adj_threshold)

    if (!only_positive) {
      deg_df <- deg_df %>%
        filter(abs(.data[[lfc_col]]) >= logfc_threshold)
    }

    write.table(
      deg_df,
      file.path(cmp_out_dir, "deg_findmarkers.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    if (nrow(deg_df) < 10) {
      message("  Skip enrichment for ", comparison_name, ": too few DEGs.")
      next
    }

    perform_enrichment_analysis(
      deg_df = deg_df,
      out_dir = cmp_out_dir,
      label = paste(celltype_name, comparison_name, "GO"),
      enrichment_type = "GO",
      organism = organism,
      org_db = org_db
    )

    perform_enrichment_analysis(
      deg_df = deg_df,
      out_dir = cmp_out_dir,
      label = paste(celltype_name, comparison_name, "KEGG"),
      enrichment_type = "KEGG",
      organism = organism,
      org_db = org_db
    )
  }
}

message("Enrichment finished: ", out_dir)
