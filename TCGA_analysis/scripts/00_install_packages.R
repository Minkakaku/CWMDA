options(stringsAsFactors = FALSE)

cran_packages <- c(
  "dplyr", "readr", "tibble", "ggplot2",
  "survival", "survminer", "pheatmap",
  "forcats"
)

bioc_packages <- c(
  "TCGAbiolinks", "SummarizedExperiment", "limma",
  "edgeR", "GSVA", "ConsensusClusterPlus",
  "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_cran) > 0L) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_bioc) > 0L) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

message("Package installation check completed.")
