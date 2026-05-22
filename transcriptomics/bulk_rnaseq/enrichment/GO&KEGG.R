#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(optparse)
})

get_script_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
        return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
    }
    getwd()
}

source(file.path(get_script_dir(), "..", "utils.R"))

load_orgdb <- function(orgdb_name) {
    if (!requireNamespace(orgdb_name, quietly = TRUE)) {
        stop("OrgDb package is required: ", orgdb_name)
    }
    get(orgdb_name, envir = asNamespace(orgdb_name))
}

select_column <- function(df, col_spec, fallback_names, column_label) {
    idx <- resolve_column(df, col_spec, fallback_names = fallback_names, column_label = column_label)
    colnames(df)[idx]
}

option_list <- list(
    make_option("--diff", type = "character", help = "Differential results table with gene identifiers."),
    make_option("--gene-col", type = "character", help = "Gene identifier column name or index."),
    make_option("--logfc-col", type = "character", help = "logFC column name or index."),
    make_option("--logfc-threshold", type = "double", default = 1),
    make_option("--padj-col", type = "character", help = "Adjusted p-value column name or index."),
    make_option("--padj-threshold", type = "double", default = 0.05),
    make_option("--keytype", type = "character", default = "SYMBOL", help = "Gene identifier type for OrgDb."),
    make_option("--orgdb", type = "character", default = "org.Mm.eg.db", help = "OrgDb package name."),
    make_option("--organism", type = "character", default = "mmu", help = "KEGG organism code."),
    make_option("--outdir", type = "character", default = "results/enrichment")
)

args <- parse_args(OptionParser(option_list = option_list))
if (is.null(args$diff)) {
    stop("--diff is required")
}

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler is required for enrichment analysis")
}

diff_df <- read_table_auto(args$diff)
gene_col <- select_column(diff_df, args$`gene-col`, c("gene_id", "gene", "Geneid", "symbol", "SYMBOL"), "gene column")
logfc_col <- select_column(diff_df, args$`logfc-col`, c("log2FoldChange", "logFC"), "logFC column")

padj_col <- NULL
if (!is.null(args$`padj-col`) || any(c("padj", "FDR", "adj.P.Val") %in% colnames(diff_df))) {
    padj_col <- select_column(diff_df, args$`padj-col`, c("padj", "FDR", "adj.P.Val"), "padj column")
}

filtered <- diff_df[abs(diff_df[[logfc_col]]) >= args$`logfc-threshold`, , drop = FALSE]
if (!is.null(padj_col)) {
    filtered <- filtered[!is.na(filtered[[padj_col]]) & filtered[[padj_col]] <= args$`padj-threshold`, , drop = FALSE]
}

gene_symbols <- unique(filtered[[gene_col]])
gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]
if (length(gene_symbols) == 0) {
    stop("No genes left after filtering. Adjust thresholds or check column names.")
}

orgdb <- load_orgdb(args$orgdb)

ego <- clusterProfiler::enrichGO(
    gene = gene_symbols,
    OrgDb = orgdb,
    keyType = args$keytype,
    ont = "ALL",
    pvalueCutoff = 0.1
)

entrez <- clusterProfiler::bitr(
    gene_symbols,
    fromType = args$keytype,
    toType = "ENTREZID",
    OrgDb = orgdb
)

kk <- clusterProfiler::enrichKEGG(
    gene = entrez$ENTREZID,
    organism = args$organism,
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05
)

dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
write.table(as.data.frame(ego), file.path(args$outdir, "go_enrichment.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(as.data.frame(kk), file.path(args$outdir, "kegg_enrichment.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(gene_symbols, file.path(args$outdir, "gene_list.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

message("GO/KEGG enrichment complete. Outputs written to: ", args$outdir)
