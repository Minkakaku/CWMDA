MANIFEST_COLUMNS <- c("sample", "group", "sample_dir", "outs_dir", "matrix_dir", "matrix_type")
OUTS_DIR_NAMES <- c("outs")
MATRIX_DIR_NAMES <- c("filtered_feature_bc_matrix")
MATRIX_FILE_NAMES <- c("filtered_feature_bc_matrix.h5")

as_path <- function(path_value) {
  normalizePath(path.expand(path_value), winslash = "/", mustWork = FALSE)
}

write_tsv <- function(df, output_path) {
  output_path <- as_path(output_path)
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, output_path, sep = "\t", row.names = FALSE, quote = FALSE)
}

read_table_auto <- function(path_value) {
  path <- as_path(path_value)
  ext <- tolower(tools::file_ext(path))
  sep <- if (ext %in% c("tsv", "txt", "xls")) "\t" else ","
  read.table(path, sep = sep, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

bind_rows_base <- function(rows) {
  if (length(rows) == 1L) {
    return(rows[[1L]])
  }
  Reduce(rbind, rows)
}

path_basename_lower <- function(path_value) {
  tolower(basename(path_value))
}

discover_samples <- function(data_root) {
  root <- as_path(data_root)
  if (!dir.exists(root)) {
    stop("Data root not found: ", root, call. = FALSE)
  }

  all_paths <- unique(c(
    list.dirs(root, recursive = TRUE, full.names = TRUE),
    list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  ))
  records <- list()

  for (path in sort(all_paths)) {
    is_matrix_dir <- dir.exists(path) && path_basename_lower(path) %in% MATRIX_DIR_NAMES
    is_matrix_file <- file.exists(path) && !dir.exists(path) && path_basename_lower(path) %in% MATRIX_FILE_NAMES
    if (!is_matrix_dir && !is_matrix_file) {
      next
    }

    outs_dir <- dirname(path)
    if (!tolower(basename(outs_dir)) %in% OUTS_DIR_NAMES) {
      next
    }

    sample_dir <- dirname(outs_dir)
    sample_name <- basename(sample_dir)
    records[[length(records) + 1L]] <- data.frame(
      sample = sample_name,
      group = sample_name,
      sample_dir = sample_dir,
      outs_dir = outs_dir,
      matrix_dir = path,
      matrix_type = if (is_matrix_file) "h5" else "mtx_dir",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  if (length(records) == 0L) {
    candidates <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
    outs_like <- candidates[dir.exists(candidates) & grepl("out", basename(candidates), ignore.case = TRUE)]
    hint <- ""
    if (length(outs_like) > 0L) {
      hint <- paste0(" Nearby outs-like directories: ", paste(head(outs_like, 10L), collapse = "; "))
    }
    stop(
      "No sample/outs/filtered_feature_bc_matrix or filtered_feature_bc_matrix.h5 paths were found under ",
      root,
      ".",
      hint,
      call. = FALSE
    )
  }

  manifest <- bind_rows_base(records)
  manifest <- manifest[!duplicated(manifest$matrix_dir), MANIFEST_COLUMNS, drop = FALSE]
  manifest <- manifest[order(manifest$sample), , drop = FALSE]
  rownames(manifest) <- NULL

  duplicate_samples <- unique(manifest$sample[duplicated(manifest$sample)])
  if (length(duplicate_samples) > 0L) {
    stop(
      "Duplicated sample names were discovered automatically: ",
      paste(duplicate_samples, collapse = ", "),
      ". Rename the sample folders or maintain sample_sheet manually.",
      call. = FALSE
    )
  }

  manifest
}

tokenize_sample_name <- function(sample_name) {
  normalized <- trimws(sample_name)
  normalized <- gsub("([a-z])([A-Z])", "\\1_\\2", normalized, perl = TRUE)
  normalized <- gsub("([A-Za-z])(\\d)", "\\1_\\2", normalized, perl = TRUE)
  normalized <- gsub("(\\d)([A-Za-z])", "\\1_\\2", normalized, perl = TRUE)
  tokens <- unlist(strsplit(normalized, "[_\\-.[:space:]]+", perl = TRUE), use.names = FALSE)
  tokens <- tokens[nzchar(tokens)]
  if (length(tokens) == 0L) sample_name else tokens
}

format_positions <- function(positions) {
  paste(positions, collapse = ",")
}

group_sizes <- function(group_map) {
  sizes <- vapply(group_map, length, integer(1))
  first_sample <- vapply(group_map, function(x) sort(x)[1], character(1))
  paste(sizes[order(-sizes, first_sample)], collapse = ",")
}

empty_suggestions <- function() {
  data.frame(
    suggestion_id = character(0),
    positions = character(0),
    group_count = integer(0),
    redundancy_saved = integer(0),
    balance_score = integer(0),
    groups = character(0),
    group_sizes = character(0),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

enumerate_group_suggestions <- function(sample_names) {
  sample_names <- as.character(sample_names)
  if (length(sample_names) < 2L) {
    return(empty_suggestions())
  }

  tokenized <- lapply(sample_names, tokenize_sample_name)
  names(tokenized) <- sample_names
  max_token_count <- max(vapply(tokenized, length, integer(1)))
  searchable_token_count <- min(max_token_count, 8L)

  suggestions <- list()
  seen_groupings <- character(0)
  suggestion_index <- 1L

  for (position_count in seq_len(searchable_token_count)) {
    combos <- combn(seq.int(0L, max_token_count - 1L), position_count, simplify = FALSE)
    for (positions in combos) {
      if (any(vapply(tokenized, function(tokens) any(positions + 1L > length(tokens)), logical(1)))) {
        next
      }

      group_names <- vapply(
        tokenized,
        function(tokens) paste(tokens[positions + 1L], collapse = "_"),
        character(1)
      )
      group_map <- split(names(group_names), group_names)
      group_count <- length(group_map)
      if (group_count <= 1L || group_count >= length(sample_names)) {
        next
      }

      grouping_signature <- paste(
        vapply(
          sort(names(group_map)),
          function(group_name) paste(group_name, paste(sort(group_map[[group_name]]), collapse = ","), sep = ":"),
          character(1)
        ),
        collapse = ";"
      )
      if (grouping_signature %in% seen_groupings) {
        next
      }
      seen_groupings <- c(seen_groupings, grouping_signature)

      sizes <- vapply(group_map, length, integer(1))
      group_text <- paste(
        vapply(
          sort(names(group_map)),
          function(group_name) paste0(group_name, ":", paste(sort(group_map[[group_name]]), collapse = ",")),
          character(1)
        ),
        collapse = "; "
      )

      suggestions[[length(suggestions) + 1L]] <- data.frame(
        suggestion_id = sprintf("S%02d", suggestion_index),
        positions = format_positions(positions),
        group_count = group_count,
        redundancy_saved = length(sample_names) - group_count,
        balance_score = max(sizes) - min(sizes),
        groups = group_text,
        group_sizes = group_sizes(group_map),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      suggestion_index <- suggestion_index + 1L
    }
  }

  if (length(suggestions) == 0L) {
    return(empty_suggestions())
  }

  suggestion_df <- bind_rows_base(suggestions)
  suggestion_df <- suggestion_df[order(
    suggestion_df$group_count,
    suggestion_df$balance_score,
    -suggestion_df$redundancy_saved,
    suggestion_df$positions,
    suggestion_df$suggestion_id
  ), , drop = FALSE]
  rownames(suggestion_df) <- NULL
  suggestion_df
}

choose_group_assignment <- function(manifest, group_count) {
  suggestions <- enumerate_group_suggestions(manifest$sample)
  sample_count <- nrow(manifest)
  if (group_count < 1L || group_count > sample_count) {
    stop("group_count must be between 1 and ", sample_count, ".", call. = FALSE)
  }

  sample_sheet <- manifest
  if (group_count == 1L) {
    sample_sheet$group <- "group1"
    return(list(sample_sheet = sample_sheet, chosen = NULL, suggestions = suggestions))
  }

  if (group_count == sample_count) {
    sample_sheet$group <- sample_sheet$sample
    return(list(sample_sheet = sample_sheet, chosen = NULL, suggestions = suggestions))
  }

  matched <- suggestions[suggestions$group_count == group_count, , drop = FALSE]
  if (nrow(matched) == 0L) {
    available <- sort(unique(suggestions$group_count))
    available_text <- if (length(available) > 0L) paste(available, collapse = ", ") else "no available suggestions"
    stop(
      "No exact auto-grouping rule can collapse samples into ",
      group_count,
      " groups. Available group counts: ",
      available_text,
      ".",
      call. = FALSE
    )
  }

  matched <- matched[order(matched$balance_score, -matched$redundancy_saved, matched$positions, matched$suggestion_id), , drop = FALSE]
  chosen <- matched[1L, , drop = FALSE]
  positions <- as.integer(strsplit(chosen$positions[[1]], ",", fixed = TRUE)[[1]])
  sample_sheet$group <- vapply(
    sample_sheet$sample,
    function(sample_name) paste(tokenize_sample_name(sample_name)[positions + 1L], collapse = "_"),
    character(1)
  )
  list(sample_sheet = sample_sheet, chosen = chosen, suggestions = suggestions)
}

parse_options <- function(args) {
  options <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    if (i == length(args)) {
      stop("Missing value for argument: ", key, call. = FALSE)
    }
    options[[sub("^--", "", key)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  options
}

require_option <- function(options, key) {
  value <- options[[key]]
  if (is.null(value) || !nzchar(value)) {
    stop("Missing required option --", key, call. = FALSE)
  }
  value
}

print_usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript sample_discovery.R discover --data-root DATA_ROOT --output-manifest manifest.tsv\n",
    "  Rscript sample_discovery.R suggest --manifest manifest.tsv --output-report suggestions.tsv\n",
    "  Rscript sample_discovery.R assign --manifest manifest.tsv --group-count N --output-sample-sheet sample_sheet.tsv [--output-report suggestions.tsv]\n",
    sep = ""
  )
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0L || args[[1]] %in% c("-h", "--help")) {
    print_usage()
    return(invisible(NULL))
  }

  command <- args[[1]]
  options <- parse_options(args[-1])

  if (identical(command, "discover")) {
    manifest <- discover_samples(require_option(options, "data-root"))
    output_manifest <- require_option(options, "output-manifest")
    write_tsv(manifest, output_manifest)
    cat("Discovered ", nrow(manifest), " samples. Manifest written to: ", as_path(output_manifest), "\n", sep = "")
    return(invisible(manifest))
  }

  if (identical(command, "suggest")) {
    manifest <- read_table_auto(require_option(options, "manifest"))
    suggestions <- enumerate_group_suggestions(manifest$sample)
    output_report <- require_option(options, "output-report")
    write_tsv(suggestions, output_report)
    if (nrow(suggestions) == 0L) {
      cat("No auto-grouping rules were found. Empty report written to: ", as_path(output_report), "\n", sep = "")
    } else {
      cat("Generated ", nrow(suggestions), " grouping suggestions. Report written to: ", as_path(output_report), "\n", sep = "")
      cat("Available group counts: ", paste(sort(unique(suggestions$group_count)), collapse = ", "), "\n", sep = "")
    }
    return(invisible(suggestions))
  }

  if (identical(command, "assign")) {
    manifest <- read_table_auto(require_option(options, "manifest"))
    group_count <- as.integer(require_option(options, "group-count"))
    if (is.na(group_count)) {
      stop("--group-count must be an integer.", call. = FALSE)
    }

    result <- choose_group_assignment(manifest, group_count)
    output_sample_sheet <- require_option(options, "output-sample-sheet")
    write_tsv(result$sample_sheet[, MANIFEST_COLUMNS, drop = FALSE], output_sample_sheet)
    cat("sample_sheet written to: ", as_path(output_sample_sheet), "\n", sep = "")

    if (!is.null(options[["output-report"]])) {
      write_tsv(result$suggestions, options[["output-report"]])
      cat("Full suggestion report written to: ", as_path(options[["output-report"]]), "\n", sep = "")
    }

    if (is.null(result$chosen)) {
      cat("The requested grouping uses the default direct assignment rule.\n")
    } else {
      cat(
        "Selected suggestion ",
        result$chosen$suggestion_id[[1]],
        " using token positions [",
        result$chosen$positions[[1]],
        "], collapsed into ",
        result$chosen$group_count[[1]],
        " groups.\n",
        sep = ""
      )
    }
    return(invisible(result))
  }

  print_usage()
  stop("Unknown command: ", command, call. = FALSE)
}

if (identical(environment(), globalenv()) && !interactive()) {
  main()
}
