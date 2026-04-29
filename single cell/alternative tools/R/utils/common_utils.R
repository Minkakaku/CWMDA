`%||%` <- function(lhs, rhs) {
  if (is.null(lhs) || length(lhs) == 0) {
    rhs
  } else {
    lhs
  }
}

get_current_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)

  if (length(file_arg) > 0) {
    return(dirname(normalizePath(
      gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE),
      winslash = "/",
      mustWork = FALSE
    )))
  }

  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(
      sys.frames()[[1]]$ofile,
      winslash = "/",
      mustWork = FALSE
    )))
  }

  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

parse_cli_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  parsed <- list()
  if (length(args) == 0) {
    return(parsed)
  }

  index <- 1
  while (index <= length(args)) {
    key <- args[[index]]
    if (!startsWith(key, "--")) {
      stop("Invalid argument format: ", key, ". Use --key value.")
    }

    key <- sub("^--", "", key)
    is_last <- index == length(args)
    next_is_flag <- !is_last && startsWith(args[[index + 1]], "--")

    if (is_last || next_is_flag) {
      parsed[[key]] <- TRUE
      index <- index + 1
    } else {
      parsed[[key]] <- args[[index + 1]]
      index <- index + 2
    }
  }

  parsed
}

cli_string <- function(args, key, default = NULL) {
  value <- args[[key]] %||% default
  if (is.null(value)) {
    return(NULL)
  }
  as.character(value)
}

cli_numeric <- function(args, key, default = NULL) {
  value <- cli_string(args, key, default)
  if (is.null(value)) {
    return(NULL)
  }
  as.numeric(value)
}

cli_bool <- function(args, key, default = FALSE) {
  value <- args[[key]]
  if (is.null(value)) {
    return(default)
  }
  if (is.logical(value)) {
    return(isTRUE(value))
  }
  tolower(as.character(value)) %in% c("1", "true", "t", "yes", "y")
}

cli_csv <- function(args, key, default = NULL) {
  value <- cli_string(args, key, default)
  if (is.null(value) || value == "") {
    return(character(0))
  }
  trimws(unlist(strsplit(value, ",", fixed = TRUE)))
}

ensure_dir <- function(path_value) {
  dir.create(path_value, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path_value, winslash = "/", mustWork = FALSE)
}

save_plot <- function(p, filename, width = 8, height = 6, dpi = 300) {
  ggplot2::ggsave(paste0(filename, ".png"), p, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(paste0(filename, ".pdf"), p, width = width, height = height)
}

save_pheatmap <- function(p, filename, width = 8, height = 6, dpi = 300) {
  grDevices::png(
    filename = paste0(filename, ".png"),
    width = width * dpi,
    height = height * dpi,
    res = dpi
  )
  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  grDevices::dev.off()

  grDevices::pdf(
    file = paste0(filename, ".pdf"),
    width = width,
    height = height
  )
  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  grDevices::dev.off()
}

join_layers_if_needed <- function(obj) {
  if (!"Seurat" %in% loadedNamespaces()) {
    return(obj)
  }
  if (!"JoinLayers" %in% getNamespaceExports("Seurat")) {
    return(obj)
  }

  tryCatch(
    Seurat::JoinLayers(obj),
    error = function(e) obj
  )
}
