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

run_rscript <- function(script_path, args) {
    cmd <- c(script_path, args)
    exit_code <- system2("Rscript", cmd)
    if (!identical(exit_code, 0L)) {
        stop("Command failed: Rscript ", paste(cmd, collapse = " "))
    }
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
    make_option("--min-counts", type = "double", default = 1),
    make_option("--min-samples", type = "integer", default = 2),
    make_option("--padj-threshold", type = "double", default = 0.05),
    make_option("--logfc-threshold", type = "double", default = 1),
    make_option("--outdir", type = "character", default = "results/bulk_rnaseq"),
    make_option("--skip-qc", action = "store_true", default = FALSE),
    make_option("--skip-differential", action = "store_true", default = FALSE),
    make_option("--run-enrichment", action = "store_true", default = FALSE),
    make_option("--run-mfuzz", action = "store_true", default = FALSE),
    make_option("--run-wgcna", action = "store_true", default = FALSE),
    make_option("--expression", type = "character", help = "Optional expression matrix for downstream analysis."),
    make_option("--time-col", type = "character", default = "time"),
    make_option("--clusters", type = "integer", default = 6),
    make_option("--wgcna-power", type = "integer", help = "Soft-thresholding power for WGCNA."),
    make_option("--wgcna-network-type", type = "character", default = "signed"),
    make_option("--wgcna-min-module-size", type = "integer", default = 30),
    make_option("--enrichment-keytype", type = "character", default = "SYMBOL"),
    make_option("--enrichment-orgdb", type = "character", default = "org.Mm.eg.db"),
    make_option("--enrichment-organism", type = "character", default = "mmu"),
    make_option("--diff-results", type = "character", help = "Optional differential results table for enrichment.")
)

args <- parse_args(OptionParser(option_list = option_list))

script_dir <- get_script_dir()
qc_outdir <- file.path(args$outdir, "expression_qc")
diff_outdir <- file.path(args$outdir, "differential_expression")
mfuzz_outdir <- file.path(args$outdir, "mfuzz")
wgcna_outdir <- file.path(args$outdir, "wgcna")
enrichment_outdir <- file.path(args$outdir, "enrichment")

if ((!(args$`skip-qc`) || !(args$`skip-differential`) || args$`run-mfuzz`) && is.null(args$`sample-sheet`)) {
    stop("--sample-sheet is required for QC, differential, or Mfuzz steps.")
}

if (!args$`skip-qc`) {
    qc_script <- file.path(script_dir, "01_counts", "main_function.r")
    qc_args <- c(
        "--sample-sheet", args$`sample-sheet`,
        "--sample-id-col", args$`sample-id-col`,
        "--group-col", args$`group-col`,
        "--min-counts", args$`min-counts`,
        "--min-samples", args$`min-samples`,
        "--outdir", qc_outdir,
        "--input-type", args$`input-type`
    )
    if (!is.null(args$counts)) {
        qc_args <- c(qc_args, "--counts", args$counts)
    }
    if (!is.null(args$`input-dir`)) {
        qc_args <- c(qc_args, "--input-dir", args$`input-dir`)
    }
    if (!is.null(args$`file-pattern`)) {
        qc_args <- c(qc_args, "--file-pattern", args$`file-pattern`)
    }
    if (!is.null(args$`count-col`)) {
        qc_args <- c(qc_args, "--count-col", args$`count-col`)
    }
    if (!is.null(args$`gene-col`)) {
        qc_args <- c(qc_args, "--gene-col", args$`gene-col`)
    }
    run_rscript(qc_script, qc_args)
}

counts_input <- if (!args$`skip-qc`) {
    file.path(qc_outdir, "counts_matrix.filtered.tsv")
} else if (!is.null(args$counts)) {
    args$counts
} else {
    stop("Counts input is required when --skip-qc is set.")
}

if (!args$`skip-differential`) {
    diff_script <- file.path(script_dir, "02_differential", "differential_expression.r")
    diff_args <- c(
        "--counts", counts_input,
        "--sample-sheet", args$`sample-sheet`,
        "--sample-id-col", args$`sample-id-col`,
        "--group-col", args$`group-col`,
        "--contrast", args$contrast,
        "--method", args$method,
        "--min-counts", args$`min-counts`,
        "--min-samples", args$`min-samples`,
        "--padj-threshold", args$`padj-threshold`,
        "--logfc-threshold", args$`logfc-threshold`,
        "--outdir", diff_outdir,
        "--input-type", "matrix"
    )
    if (!is.null(args$`gene-col`)) {
        diff_args <- c(diff_args, "--gene-col", args$`gene-col`)
    }
    run_rscript(diff_script, diff_args)
}

expression_input <- if (!is.null(args$expression)) {
    args$expression
} else if (!args$`skip-qc`) {
    file.path(qc_outdir, "log2cpm.filtered.tsv")
} else {
    stop("Expression input is required when --skip-qc is set.")
}

if (args$`run-mfuzz`) {
    mfuzz_script <- file.path(script_dir, "pseudotime", "mfuzz", "mfuzz.r")
    mfuzz_args <- c(
        "--expression", expression_input,
        "--sample-sheet", args$`sample-sheet`,
        "--sample-id-col", args$`sample-id-col`,
        "--time-col", args$`time-col`,
        "--clusters", args$clusters,
        "--outdir", mfuzz_outdir
    )
    if (!is.null(args$`gene-col`)) {
        mfuzz_args <- c(mfuzz_args, "--gene-col", args$`gene-col`)
    }
    run_rscript(mfuzz_script, mfuzz_args)
}

if (args$`run-wgcna`) {
    wgcna_script <- file.path(script_dir, "wgcna", "wgcna.r")
    wgcna_args <- c(
        "--expression", expression_input,
        "--network-type", args$`wgcna-network-type`,
        "--min-module-size", args$`wgcna-min-module-size`,
        "--outdir", wgcna_outdir
    )
    if (!is.null(args$`wgcna-power`)) {
        wgcna_args <- c(wgcna_args, "--power", args$`wgcna-power`)
    }
    if (!is.null(args$`gene-col`)) {
        wgcna_args <- c(wgcna_args, "--gene-col", args$`gene-col`)
    }
    run_rscript(wgcna_script, wgcna_args)
}

if (args$`run-enrichment`) {
    enrichment_script <- file.path(script_dir, "enrichment", "GO&KEGG.R")
    diff_method_dir <- file.path(diff_outdir, tolower(args$method))
    diff_results <- if (!is.null(args$`diff-results`)) {
        args$`diff-results`
    } else {
        file.path(diff_method_dir, "all_results.tsv")
    }
    if (!file.exists(diff_results)) {
        stop("Differential results not found: ", diff_results)
    }
    enrich_args <- c(
        "--diff", diff_results,
        "--logfc-threshold", args$`logfc-threshold`,
        "--padj-threshold", args$`padj-threshold`,
        "--keytype", args$`enrichment-keytype`,
        "--orgdb", args$`enrichment-orgdb`,
        "--organism", args$`enrichment-organism`,
        "--outdir", enrichment_outdir
    )
    if (!is.null(args$`gene-col`)) {
        enrich_args <- c(enrich_args, "--gene-col", args$`gene-col`)
    }
    run_rscript(enrichment_script, enrich_args)
}

message("Bulk RNA-seq analysis complete. Outputs written to: ", args$outdir)
