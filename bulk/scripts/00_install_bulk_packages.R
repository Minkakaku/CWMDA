options(stringsAsFactors = FALSE)

cran_packages <- c(
  "data.table", "ggplot2", "dplyr", "tibble",
  "tidyr", "stringr", "pheatmap", "e1071",
  "readxl"
)

bioc_packages <- c(
  "preprocessCore", "fgsea"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0L) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0L) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

message("Bulk package installation check completed.")
