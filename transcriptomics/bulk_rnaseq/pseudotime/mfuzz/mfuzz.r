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

option_list <- list(
    make_option("--expression", type = "character", help = "Expression matrix with gene_id in first column."),
    make_option("--sample-sheet", type = "character", help = "Sample sheet with sample_id and time columns."),
    make_option("--sample-id-col", type = "character", default = "sample_id"),
    make_option("--time-col", type = "character", default = "time"),
    make_option("--gene-col", type = "character", help = "Gene ID column name or index."),
    make_option("--clusters", type = "integer", default = 6),
    make_option("--outdir", type = "character", default = "results/mfuzz")
)

args <- parse_args(OptionParser(option_list = option_list))
if (is.null(args$expression) || is.null(args$`sample-sheet`)) {
    stop("--expression and --sample-sheet are required")
}

if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase is required for mfuzz analysis")
}
if (!requireNamespace("Mfuzz", quietly = TRUE)) {
    stop("Mfuzz is required for mfuzz analysis")
}

sample_sheet <- read_table_auto(args$`sample-sheet`)
sample_sheet <- validate_required_columns(sample_sheet, c(args$`sample-id-col`, args$`time-col`))

expr <- read_expression_matrix(args$expression, args$`gene-col`)
expr <- order_counts_by_samples(expr, sample_sheet[[args$`sample-id-col`]])

pheno <- data.frame(
    time = sample_sheet[[args$`time-col`]],
    row.names = sample_sheet[[args$`sample-id-col`]]
)

eset <- Biobase::ExpressionSet(assayData = as.matrix(expr), phenoData = Biobase::AnnotatedDataFrame(pheno))

eset <- Mfuzz::filter.NA(eset)
eset <- Mfuzz::fill.NA(eset, mode = "mean")
eset <- Mfuzz::standardise(eset)

m <- Mfuzz::mestimate(eset)
cl <- Mfuzz::mfuzz(eset, c = args$clusters, m = m)

dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

membership <- data.frame(
    gene_id = rownames(expr),
    cluster = cl$cluster,
    stringsAsFactors = FALSE
)
write.table(membership, file.path(args$outdir, "cluster_membership.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

centers <- as.data.frame(cl$centers)
centers$cluster <- seq_len(nrow(centers))
write.table(centers, file.path(args$outdir, "cluster_centers.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

message("Mfuzz clustering complete. Outputs written to: ", args$outdir)
