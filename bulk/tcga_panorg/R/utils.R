`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) {
    return(y)
  }
  if (length(x) == 1L && is.na(x)) {
    return(y)
  }
  x
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

project_path <- function(...) {
  normalizePath(file.path(getwd(), ...), winslash = "/", mustWork = FALSE)
}

split_semicolon_field <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) {
    return(character())
  }
  trimws(unlist(strsplit(x, ";", fixed = TRUE)))
}

clean_tcga_patient_id <- function(x) {
  x <- toupper(gsub("\\.", "-", as.character(x)))
  x <- ifelse(nchar(x) >= 12L, substr(x, 1L, 12L), x)
  x[x %in% c("NA", "")] <- NA_character_
  x
}

first_present_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0L) {
    return(NULL)
  }
  hit[[1]]
}

coalesce_nonempty <- function(df, candidates) {
  if (nrow(df) == 0L) {
    return(character())
  }
  out <- rep(NA_character_, nrow(df))
  for (col in candidates) {
    if (!col %in% colnames(df)) {
      next
    }
    value <- as.character(df[[col]])
    value <- trimws(value)
    value[value == ""] <- NA_character_
    idx <- is.na(out) & !is.na(value)
    out[idx] <- value[idx]
  }
  out
}

days_to_years <- function(x) {
  value <- suppressWarnings(as.numeric(x))
  if (all(is.na(value))) {
    return(rep(NA_real_, length(value)))
  }
  value <- ifelse(abs(value) > 200, abs(value) / 365.25, value)
  value
}

collapse_stage <- function(x) {
  value <- toupper(trimws(as.character(x)))
  value <- gsub("^STAGE ", "", value)
  value <- gsub("[ABC]$", "", value)
  value[value %in% c("", "NA", "[NOT AVAILABLE]", "NOT REPORTED", "STAGE X")] <- NA_character_
  value
}

aggregate_expression_by_gene <- function(mat, genes) {
  keep <- !is.na(genes) & nzchar(genes)
  mat <- mat[keep, , drop = FALSE]
  genes <- genes[keep]
  mode(mat) <- "numeric"
  summed <- rowsum(mat, group = genes, reorder = FALSE)
  counts <- as.numeric(table(genes)[rownames(summed)])
  sweep(summed, 1L, counts, "/")
}

merge_expression_matrices <- function(mats) {
  stopifnot(length(mats) >= 1L)
  all_genes <- unique(unlist(lapply(mats, rownames)))
  merged_list <- lapply(mats, function(mat) {
    idx <- match(all_genes, rownames(mat))
    out <- matrix(NA_real_, nrow = length(all_genes), ncol = ncol(mat))
    rownames(out) <- all_genes
    colnames(out) <- colnames(mat)
    valid <- !is.na(idx)
    out[valid, ] <- mat[idx[valid], , drop = FALSE]
    out
  })
  merged <- do.call(cbind, merged_list)
  rownames(merged) <- all_genes
  merged
}

deduplicate_patient_matrix <- function(mat) {
  patient_ids <- clean_tcga_patient_id(colnames(mat))
  keep <- !duplicated(patient_ids)
  mat <- mat[, keep, drop = FALSE]
  colnames(mat) <- patient_ids[keep]
  mat
}

read_panorgs <- function(file) {
  unique(scan(file, what = "character", quiet = TRUE))
}

write_expression_tsv <- function(mat, file) {
  df <- data.frame(gene_symbol = rownames(mat), mat, check.names = FALSE)
  readr::write_tsv(df, file)
}

read_expression_tsv <- function(file, gene_col = "gene_symbol") {
  df <- readr::read_tsv(file, show_col_types = FALSE)
  if (!gene_col %in% colnames(df)) {
    stop("Expression file is missing gene column: ", gene_col)
  }
  genes <- df[[gene_col]]
  expr <- data.matrix(df[, setdiff(colnames(df), gene_col), drop = FALSE])
  rownames(expr) <- genes
  expr
}

apply_expr_transform <- function(mat, transform = "log2p1") {
  if (is.null(transform) || is.na(transform) || transform %in% c("", "none")) {
    return(mat)
  }
  if (transform == "log2p1") {
    return(log2(mat + 1))
  }
  stop("Unsupported transform: ", transform)
}

scale_rows <- function(mat) {
  scaled <- t(scale(t(mat)))
  scaled[is.na(scaled)] <- 0
  scaled
}

sanitize_feature_name <- function(x) {
  gsub("_+", "_", gsub("[^A-Za-z0-9]+", "_", x))
}

parse_event_vector <- function(x, event_positive = c("1", "TRUE", "DEAD", "DECEASED")) {
  value <- toupper(trimws(as.character(x)))
  out <- ifelse(value %in% toupper(event_positive), 1L, 0L)
  out[is.na(value) | value == ""] <- NA_integer_
  out
}

detect_barcode_column <- function(df) {
  if (nrow(df) == 0L) {
    return(NULL)
  }
  score <- vapply(colnames(df), function(col) {
    value <- as.character(df[[col]])
    value <- value[!is.na(value)]
    if (length(value) == 0L) {
      return(0)
    }
    mean(grepl("^TCGA[-.][A-Z0-9]{2}[-.][A-Z0-9]{4}", toupper(value)))
  }, numeric(1))
  if (max(score) == 0) {
    return(NULL)
  }
  names(score)[which.max(score)]
}

pick_subtype_column <- function(df, barcode_col = NULL) {
  preferred <- c(
    "Subtype_Selected", "subtype_selected", "Subtype", "subtype",
    "mRNA_subtype", "Subtype_mRNA", "iCluster", "Pan.GI.Cluster",
    "paper_BRCA_Subtype_PAM50"
  )
  hit <- preferred[preferred %in% colnames(df)]
  if (length(hit) > 0L) {
    return(hit[[1]])
  }
  candidates <- setdiff(colnames(df), barcode_col %||% character())
  if (length(candidates) == 0L) {
    return(NULL)
  }
  score <- vapply(candidates, function(col) {
    value <- unique(na.omit(trimws(as.character(df[[col]]))))
    n <- length(value)
    if (n < 2L || n > 10L) {
      return(Inf)
    }
    n
  }, numeric(1))
  if (all(!is.finite(score))) {
    return(NULL)
  }
  names(score)[which.min(score)]
}

save_ggplot <- function(plot_obj, file, width = 7, height = 5, dpi = 300) {
  ggplot2::ggsave(filename = file, plot = plot_obj, width = width, height = height, dpi = dpi)
}
