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
    make_option("--counts", type = "character", help = "Count matrix with gene_id in first column."),
    make_option("--input-type", type = "character", default = "matrix", help = "matrix, rsem, or featurecounts."),
    make_option("--input-dir", type = "character", help = "Directory with per-sample count files."),
    make_option("--file-pattern", type = "character", help = "Regex pattern to select count files."),
    make_option("--count-col", type = "character", help = "Count column name or index for per-sample files."),
    make_option("--gene-col", type = "character", help = "Gene ID column name or index."),
    make_option("--sample-sheet", type = "character", help = "Sample sheet with sample_id and group columns."),
    make_option("--sample-id-col", type = "character", default = "sample_id"),
    make_option("--group-col", type = "character", default = "group"),
    make_option("--min-counts", type = "double", default = 1),
    make_option("--min-samples", type = "integer", default = 2),
    make_option("--outdir", type = "character", default = "results/expression_qc")
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

outdir <- args$outdir
write_counts_matrix(counts, outdir, "counts_matrix.raw.tsv")

library_sizes <- colSums(counts)
expressed_genes <- colSums(counts > 0)
qc_table <- data.frame(
    sample_id = sample_ids,
    group = sample_sheet[[args$`group-col`]],
    library_size = library_sizes,
    detected_genes = expressed_genes,
    stringsAsFactors = FALSE
)

write.table(qc_table, file.path(outdir, "sample_qc.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

filtered_counts <- filter_counts(counts, args$`min-counts`, args$`min-samples`)
write_counts_matrix(filtered_counts, outdir, "counts_matrix.filtered.tsv")

log2_cpm <- calc_log2_cpm(filtered_counts)
write_matrix(log2_cpm, outdir, "log2cpm.filtered.tsv")

metadata <- sample_sheet[, c(args$`sample-id-col`, args$`group-col`), drop = FALSE]
write.table(metadata, file.path(outdir, "sample_metadata.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

message("QC complete. Outputs written to: ", outdir)
