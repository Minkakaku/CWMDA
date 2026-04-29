get_current_script_dir <- function() {
  normalize_script_dir <- function(path) {
    dirname(normalizePath(
      path,
      winslash = .Platform$file.sep,
      mustWork = FALSE
    ))
  }

  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)

  if (length(file_arg) > 0) {
    return(normalize_script_dir(gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE)))
  }

  file_arg_index <- match("--file", cmd_args)
  if (!is.na(file_arg_index) && length(cmd_args) >= file_arg_index + 1) {
    return(normalize_script_dir(cmd_args[file_arg_index + 1]))
  }

  normalizePath(getwd(), winslash = .Platform$file.sep, mustWork = FALSE)
}

find_single_cell_project_dir <- function(script_dir) {
  candidates <- c(
    file.path(script_dir, "..", "..", ".."),
    file.path(script_dir, "..", ".."),
    file.path(script_dir, ".."),
    script_dir,
    getwd()
  )

  for (candidate in candidates) {
    candidate <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
    if (dir.exists(file.path(candidate, "04_qs2")) ||
        dir.exists(file.path(candidate, "R single cell")) ||
        file.exists(file.path(candidate, "README.md"))) {
      return(candidate)
    }
  }

  normalizePath(file.path(script_dir, "..", "..", ".."), winslash = "/", mustWork = FALSE)
}

log_timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

format_duration <- function(seconds) {
  seconds <- as.numeric(seconds)

  if (is.na(seconds) || !is.finite(seconds)) {
    return("unknown")
  }

  if (seconds < 60) {
    return(sprintf("%.1fs", seconds))
  }

  minutes <- floor(seconds / 60)
  sprintf("%dm %.1fs", minutes, seconds - minutes * 60)
}

log_info <- function(...) {
  message("[", log_timestamp(), "] [INFO] ", ..., appendLF = TRUE)
}

log_warn <- function(...) {
  message("[", log_timestamp(), "] [WARN] ", ..., appendLF = TRUE)
}

log_step <- function(label, expr) {
  start_time <- Sys.time()
  log_info("START ", label)

  tryCatch(
    {
      result <- force(expr)
      elapsed <- difftime(Sys.time(), start_time, units = "secs")
      log_info("DONE ", label, " (", format_duration(elapsed), ")")
      result
    },
    error = function(e) {
      elapsed <- difftime(Sys.time(), start_time, units = "secs")
      message(
        "[",
        log_timestamp(),
        "] [ERROR] FAILED ",
        label,
        " (",
        format_duration(elapsed),
        "): ",
        conditionMessage(e),
        appendLF = TRUE
      )
      stop(e)
    }
  )
}

require_packages <- function(packages) {
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0) {
    stop(
      "Missing required R packages. Please install them before running: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(lapply(packages, library, character.only = TRUE))
}

first_existing_file <- function(paths) {
  existing <- paths[file.exists(paths)]

  if (length(existing) == 0) {
    stop("No input QS2 file found. Checked: ", paste(paths, collapse = ", "))
  }

  existing[[1]]
}

get_counts_matrix <- function(obj) {
  LayerData(obj, assay = "RNA", layer = "counts")
}

script_dir <- get_current_script_dir()
project_dir <- find_single_cell_project_dir(script_dir)
setwd(project_dir)

require_packages(c("Seurat", "anndata", "Matrix", "qs2"))
log_info("Toh5ad conversion started")

input_override <- Sys.getenv("TOH5AD_INPUT_QS2", unset = "")
if (nzchar(input_override)) {
  input_qs2 <- input_override
  if (!file.exists(input_qs2)) {
    stop("TOH5AD_INPUT_QS2 does not exist: ", input_qs2)
  }
} else {
  input_qs2 <- first_existing_file(c(
    file.path("04_qs2", "05.sc_final.qs2"),
    file.path("04_qs2", "04.sc_annotated.qs2"),
    file.path("04_qs2", "01.sc_annoted.qs2")
  ))
}

output_h5ad <- Sys.getenv("TOH5AD_OUTPUT_H5AD", unset = "")
if (!nzchar(output_h5ad)) {
  output_h5ad <- sub("\\.qs2$", ".h5ad", input_qs2, ignore.case = TRUE)
}

log_info("Input QS2: ", input_qs2)
log_info("Output h5ad: ", output_h5ad)
sc <- log_step("Read Seurat QS2", qs_read(input_qs2))

if (!inherits(sc, "Seurat")) {
  stop("Input QS2 does not contain a Seurat object: ", input_qs2)
}

if (!"RNA" %in% Assays(sc)) {
  stop("Input Seurat object is missing RNA assay.")
}

dimreducs <- intersect(c("pca", "umap"), names(sc@reductions))
seurat_rna <- log_step(
  "DietSeurat RNA assay",
  DietSeurat(
    sc,
    assays = "RNA",
    dimreducs = dimreducs
  )
)

DefaultAssay(seurat_rna) <- "RNA"

seurat_rna <- log_step("JoinLayers RNA assay", JoinLayers(seurat_rna, assay = "RNA"))

counts_mat <- log_step("Extract RNA counts matrix", get_counts_matrix(seurat_rna))
if (nrow(counts_mat) == 0 || ncol(counts_mat) == 0) {
  stop("RNA counts matrix is empty.")
}
log_info("RNA counts matrix: genes=", nrow(counts_mat), "; cells=", ncol(counts_mat))

gene_names <- rownames(counts_mat)
if (is.null(gene_names) || anyNA(gene_names) || any(!nzchar(gene_names))) {
  stop("RNA counts matrix is missing valid gene names.")
}

gene_meta <- data.frame(
  gene_ids = gene_names,
  row.names = gene_names,
  check.names = FALSE
)

ad <- log_step(
  "Build AnnData object",
  AnnData(
    X = t(counts_mat),
    obs = seurat_rna@meta.data,
    var = gene_meta
  )
)

if ("pca" %in% names(seurat_rna@reductions)) {
  ad$obsm$X_pca <- Embeddings(seurat_rna, "pca")
}
if ("umap" %in% names(seurat_rna@reductions)) {
  ad$obsm$X_umap <- Embeddings(seurat_rna, "umap")
}

dir.create(dirname(output_h5ad), recursive = TRUE, showWarnings = FALSE)
log_step("Write h5ad", ad$write_h5ad(output_h5ad))

log_info("Toh5ad conversion finished: ", output_h5ad)
