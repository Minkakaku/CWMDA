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
    make_option("--gene-col", type = "character", help = "Gene ID column name or index."),
    make_option("--power", type = "integer", help = "Soft-thresholding power. If not set, pickSoftThreshold is used."),
    make_option("--network-type", type = "character", default = "signed"),
    make_option("--min-module-size", type = "integer", default = 30),
    make_option("--outdir", type = "character", default = "results/wgcna")
)

args <- parse_args(OptionParser(option_list = option_list))
if (is.null(args$expression)) {
    stop("--expression is required")
}

if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("WGCNA is required for WGCNA analysis")
}

expr <- read_expression_matrix(args$expression, args$`gene-col`)
expr <- t(expr)

WGCNA::allowWGCNAThreads()

gsg <- WGCNA::goodSamplesGenes(expr, verbose = 3)
if (!gsg$allOK) {
    expr <- expr[gsg$goodSamples, gsg$goodGenes]
}

power <- args$power
if (is.null(power)) {
    sft <- WGCNA::pickSoftThreshold(expr, powerVector = 1:20, verbose = 0)
    power <- sft$powerEstimate
    if (is.na(power)) {
        power <- 6
    }
}

net <- WGCNA::blockwiseModules(
    expr,
    power = power,
    networkType = args$`network-type`,
    minModuleSize = args$`min-module-size`,
    reassignThreshold = 0,
    mergeCutHeight = 0.25,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    verbose = 3
)

module_colors <- WGCNA::labels2colors(net$colors)
module_assignment <- data.frame(
    gene_id = colnames(expr),
    module = module_colors,
    stringsAsFactors = FALSE
)

module_eigengenes <- WGCNA::moduleEigengenes(expr, colors = module_colors)$eigengenes
module_eigengenes <- data.frame(sample_id = rownames(module_eigengenes), module_eigengenes, check.names = FALSE)

summary_df <- data.frame(
    power = power,
    network_type = args$`network-type`,
    min_module_size = args$`min-module-size`,
    module_count = length(unique(module_colors))
)

dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
write.table(module_assignment, file.path(args$outdir, "module_assignment.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(module_eigengenes, file.path(args$outdir, "module_eigengenes.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(summary_df, file.path(args$outdir, "summary.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

message("WGCNA complete. Outputs written to: ", args$outdir)
