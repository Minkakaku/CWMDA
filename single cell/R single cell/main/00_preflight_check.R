rm(list = ls())

cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(cmd_file) > 0) {
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", cmd_file[1]), fixed = TRUE), winslash = "/", mustWork = FALSE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
pipeline_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(pipeline_dir, "00_config.R"))) {
  pipeline_dir <- normalizePath(script_dir, winslash = "/", mustWork = FALSE)
}

script_dir <- pipeline_dir
source(file.path(pipeline_dir, "00_config.R"), local = FALSE)

errors <- character(0)
warnings <- character(0)

add_error <- function(...) errors <<- c(errors, paste0(...))
add_warning <- function(...) warnings <<- c(warnings, paste0(...))

missing_packages <- PACKAGES[!vapply(PACKAGES, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  add_error("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

missing_optional <- OPTIONAL_PACKAGES[!vapply(OPTIONAL_PACKAGES, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_optional) > 0) {
  add_warning("Missing optional packages: ", paste(missing_optional, collapse = ", "))
}

if (identical(CLUSTER$batch_method, "harmony") && !requireNamespace("harmony", quietly = TRUE)) {
  add_error("PIPELINE_BATCH_METHOD='harmony' requires package harmony.")
}
if (identical(SUBCLUSTER$batch_method, "harmony") && !requireNamespace("harmony", quietly = TRUE)) {
  add_error("SUBCLUSTER_BATCH_METHOD='harmony' requires package harmony.")
}
if (PROPORTION$make_sankey && !requireNamespace("ggalluvial", quietly = TRUE)) {
  add_warning("PROPORTION_MAKE_SANKEY=TRUE but ggalluvial is not installed; Sankey plots will be skipped.")
}

r_scripts <- list.files(pipeline_dir, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
for (script in r_scripts) {
  tryCatch(
    parse(script),
    error = function(e) add_error("R syntax error in ", script, ": ", conditionMessage(e))
  )
}

if (!dir.exists(DATA_DIR)) {
  add_error("Data directory does not exist: ", DATA_DIR)
} else {
  samples <- list.dirs(DATA_DIR, recursive = FALSE, full.names = FALSE)
  if (length(samples) == 0) {
    add_error("No sample dirs under ", DATA_DIR)
  }

  required_10x <- c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")
  for (sample in samples) {
    matrix_dir <- file.path(DATA_DIR, sample, "outs", "filtered_feature_bc_matrix")
    if (!dir.exists(matrix_dir)) {
      add_warning("Sample missing matrix dir: ", sample)
      next
    }

    missing_files <- required_10x[!file.exists(file.path(matrix_dir, required_10x))]
    if (length(missing_files) > 0) {
      add_error("Sample ", sample, " missing 10x files: ", paste(missing_files, collapse = ", "))
    }
  }
}

if (length(warnings) > 0) {
  message("Preflight warnings:")
  for (x in warnings) message("- ", x)
}

if (length(errors) > 0) {
  message("Preflight errors:")
  for (x in errors) message("- ", x)
  stop("Preflight failed.", call. = FALSE)
}

message("Preflight check passed.")
