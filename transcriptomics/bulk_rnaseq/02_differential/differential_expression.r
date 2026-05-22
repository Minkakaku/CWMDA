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

write_results <- function(res_df, deg, args, outdir, tool_name) {
    write.table(res_df, file.path(outdir, "all_results.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    write.table(deg, file.path(outdir, "differential_genes.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    summary_df <- data.frame(
        tool = tool_name,
        total_genes = nrow(res_df),
        differential_genes = nrow(deg),
        padj_threshold = args$`padj-threshold`,
        logfc_threshold = args$`logfc-threshold`,
        contrast = paste(args$contrast_parts, collapse = " vs ")
    )
    write.table(summary_df, file.path(outdir, "summary.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
}

filter_deg <- function(res_df, padj_col, logfc_col, args) {
    subset(
        res_df,
        !is.na(res_df[[padj_col]]) & res_df[[padj_col]] <= args$`padj-threshold` &
            abs(res_df[[logfc_col]]) >= args$`logfc-threshold`
    )
}

run_deseq2 <- function(counts, sample_metadata, args, outdir) {
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop("DESeq2 is required for --method deseq2")
    }

    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = counts,
        colData = sample_metadata,
        design = ~ group
    )
    dds <- DESeq2::DESeq(dds)
    res <- DESeq2::results(dds, contrast = c("group", args$contrast_parts[1], args$contrast_parts[2]))

    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)
    res_df <- res_df[order(res_df$pvalue), ]

    deg <- filter_deg(res_df, "padj", "log2FoldChange", args)
    write_results(res_df, deg, args, outdir, "DESeq2")
}

run_edger <- function(counts, sample_metadata, args, outdir) {
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("edgeR is required for --method edger")
    }

    y <- edgeR::DGEList(counts = counts, group = sample_metadata$group)
    keep <- rowSums(edgeR::cpm(y) >= args$`min-counts`) >= args$`min-samples`
    y <- y[keep, , keep.lib.sizes = FALSE]

    y <- edgeR::calcNormFactors(y)

    et <- edgeR::exactTest(y, pair = c(args$contrast_parts[2], args$contrast_parts[1]))
    res <- edgeR::topTags(et, n = Inf)$table

    res$gene_id <- rownames(res)
    res <- res[order(res$PValue), ]

    deg <- filter_deg(res, "FDR", "logFC", args)
    write_results(res, deg, args, outdir, "edgeR")
}

run_limma <- function(counts, sample_metadata, args, outdir) {
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("edgeR is required for --method limma")
    }
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("limma is required for --method limma")
    }

    y <- edgeR::DGEList(counts = counts)
    keep <- rowSums(edgeR::cpm(y) >= args$`min-counts`) >= args$`min-samples`
    y <- y[keep, , keep.lib.sizes = FALSE]

    y <- edgeR::calcNormFactors(y)

    design <- model.matrix(~0 + group, data = sample_metadata)
    colnames(design) <- levels(sample_metadata$group)

    v <- limma::voom(y, design, plot = FALSE)
    fit <- limma::lmFit(v, design)
    contrast_formula <- paste0(args$contrast_parts[1], " - ", args$contrast_parts[2])
    contrast_matrix <- limma::makeContrasts(contrasts = contrast_formula, levels = design)
    fit <- limma::contrasts.fit(fit, contrast_matrix)
    fit <- limma::eBayes(fit)

    res <- limma::topTable(fit, number = Inf, sort.by = "P")
    res$gene_id <- rownames(res)
    res <- res[order(res$P.Value), ]

    deg <- filter_deg(res, "adj.P.Val", "logFC", args)
    write_results(res, deg, args, outdir, "limma")
}

option_list <- list(
    make_option("--counts", type = "character", help = "Count matrix with gene_id in first column."),
    make_option("--input-type", type = "character", default = "matrix", help = "matrix, rsem, or featurecounts."),
    make_option("--input-dir", type = "character", help = "Directory with per-sample count files."),
    make_option("--file-pattern", type = "character", help = "Regex pattern to select count files."),
    make_option("--count-col", type = "character", help = "Count column name or index for per-sample files."),
    make_option("--gene-col", type = "character", help = "Gene ID column name or index."),
    make_option("--sample-sheet", type = "character", help = "Sample sheet with sample_id and group columns."),
    make_option("--sample-id-col", type = "character", default = "sample_id"),
    make_option("--group-col", type = "character", default = "group"),
    make_option("--contrast", type = "character", default = "case,control"),
    make_option("--method", type = "character", default = "deseq2", help = "deseq2, edger, or limma"),
    make_option("--outdir", type = "character", default = "results/differential_expression"),
    make_option("--min-counts", type = "double", default = 1),
    make_option("--min-samples", type = "integer", default = 2),
    make_option("--padj-threshold", type = "double", default = 0.05),
    make_option("--logfc-threshold", type = "double", default = 1)
)

args <- parse_args(OptionParser(option_list = option_list))
if (is.null(args$`sample-sheet`)) {
    stop("--sample-sheet is required")
}

sample_sheet <- read_table_auto(args$`sample-sheet`)
sample_sheet <- validate_sample_sheet(sample_sheet, args$`sample-id-col`, args$`group-col`)
sample_ids <- sample_sheet[[args$`sample-id-col`]]

if (args$`input-type` == "matrix") {
    if (is.null(args$counts)) {
        stop("--counts is required when --input-type is matrix")
    }
    counts <- read_counts_matrix(args$counts, args$`gene-col`)
} else {
    counts <- build_counts_from_dir(
        args$`input-dir`,
        sample_ids,
        args$`input-type`,
        args$`file-pattern`,
        args$`gene-col`,
        args$`count-col`
    )
}

counts <- order_counts_by_samples(counts, sample_ids)
counts <- filter_counts(counts, args$`min-counts`, args$`min-samples`)

args$contrast_parts <- parse_contrast(args$contrast)

sample_metadata <- data.frame(
    group = factor(
        sample_sheet[[args$`group-col`]],
        levels = c(args$contrast_parts[2], args$contrast_parts[1])
    ),
    row.names = sample_ids
)

method <- tolower(args$method)
outdir <- file.path(args$outdir, method)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write_counts_matrix(counts, outdir, "counts_matrix.filtered.tsv")

if (method == "deseq2") {
    run_deseq2(counts, sample_metadata, args, outdir)
} else if (method == "edger") {
    run_edger(counts, sample_metadata, args, outdir)
} else if (method == "limma") {
    run_limma(counts, sample_metadata, args, outdir)
} else {
    stop("Unsupported --method: ", args$method)
}

message("Differential expression complete. Outputs written to: ", outdir)
