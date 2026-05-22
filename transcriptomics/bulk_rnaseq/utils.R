read_table_auto <- function(path) {
    ext <- tolower(tools::file_ext(path))
    if (ext == "csv") {
        return(read.csv(path, check.names = FALSE))
    }
    read.delim(path, check.names = FALSE)
}

resolve_column <- function(df, col_spec, fallback_names = NULL, fallback_index = NULL, column_label = "column") {
    if (!is.null(col_spec) && nzchar(col_spec)) {
        suppressWarnings({
            idx <- as.integer(col_spec)
        })
        if (!is.na(idx)) {
            if (idx < 1 || idx > ncol(df)) {
                stop(column_label, " index out of range: ", col_spec)
            }
            return(idx)
        }
        if (col_spec %in% colnames(df)) {
            return(match(col_spec, colnames(df)))
        }
        stop(column_label, " not found: ", col_spec)
    }
    if (!is.null(fallback_names)) {
        for (name in fallback_names) {
            if (name %in% colnames(df)) {
                return(match(name, colnames(df)))
            }
        }
    }
    if (!is.null(fallback_index)) {
        if (fallback_index < 1 || fallback_index > ncol(df)) {
            stop(column_label, " fallback index out of range: ", fallback_index)
        }
        return(fallback_index)
    }
    stop("Unable to resolve ", column_label)
}

parse_contrast <- function(contrast) {
    contrast_parts <- strsplit(contrast, ",")[[1]]
    if (length(contrast_parts) != 2) {
        stop("--contrast must be formatted as 'case,control'")
    }
    contrast_parts
}

validate_sample_sheet <- function(sample_sheet, sample_id_col, group_col) {
    if (!(sample_id_col %in% colnames(sample_sheet))) {
        stop("sample_id column not found in sample sheet")
    }
    if (!(group_col %in% colnames(sample_sheet))) {
        stop("group column not found in sample sheet")
    }
    sample_sheet[!duplicated(sample_sheet[[sample_id_col]]), ]
}

validate_required_columns <- function(sample_sheet, required_cols) {
    missing_cols <- setdiff(required_cols, colnames(sample_sheet))
    if (length(missing_cols) > 0) {
        stop("Missing required columns in sample sheet: ", paste(missing_cols, collapse = ", "))
    }
    sample_sheet
}

get_sample_file <- function(input_dir, sample_id, file_pattern) {
    files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)
    if (length(files) == 0) {
        stop("No files found in input directory: ", input_dir)
    }
    basenames <- basename(files)
    exact_matches <- files[grepl(paste0("^", sample_id), basenames)]
    if (length(exact_matches) == 0) {
        exact_matches <- files[grepl(sample_id, basenames)]
    }
    if (length(exact_matches) == 0) {
        stop("No file found for sample_id ", sample_id, " in ", input_dir)
    }
    if (length(exact_matches) > 1) {
        stop("Multiple files found for sample_id ", sample_id, ": ", paste(basenames[grepl(sample_id, basenames)], collapse = ", "))
    }
    exact_matches[1]
}

read_count_file <- function(file_path, sample_id, gene_col, count_col, count_fallback_names, count_fallback_index) {
    df <- read_table_auto(file_path)
    gene_idx <- resolve_column(df, gene_col, fallback_names = c("gene_id", "gene", "Geneid"), fallback_index = 1, column_label = "gene column")
    count_idx <- resolve_column(
        df,
        count_col,
        fallback_names = count_fallback_names,
        fallback_index = count_fallback_index,
        column_label = "count column"
    )
    data.frame(
        gene_id = df[[gene_idx]],
        !!sample_id := as.numeric(df[[count_idx]]),
        check.names = FALSE
    )
}

merge_counts <- function(counts_list) {
    merged <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), counts_list)
    merged[is.na(merged)] <- 0
    merged
}

build_counts_from_dir <- function(input_dir, sample_ids, input_type, file_pattern, gene_col, count_col) {
    if (is.null(input_dir)) {
        stop("--input-dir is required when using --input-type rsem or featurecounts")
    }
    if (input_type == "rsem") {
        pattern <- ifelse(is.null(file_pattern), "genes\\.results$", file_pattern)
        count_names <- c("expected_count")
        fallback_index <- NULL
    } else if (input_type == "featurecounts") {
        pattern <- ifelse(is.null(file_pattern), "featureCounts", file_pattern)
        count_names <- c("Counts", "count")
        fallback_index <- NULL
    } else {
        stop("Unsupported input type: ", input_type)
    }
    counts_list <- lapply(sample_ids, function(sample_id) {
        file_path <- get_sample_file(input_dir, sample_id, pattern)
        df <- read_table_auto(file_path)
        fallback_index <- if (input_type == "featurecounts") ncol(df) else NULL
        read_count_file(file_path, sample_id, gene_col, count_col, count_names, fallback_index)
    })
    merged <- merge_counts(counts_list)
    merged <- merged[order(merged$gene_id), ]
    rownames(merged) <- merged$gene_id
    merged[, -1, drop = FALSE]
}

read_counts_matrix <- function(counts_path, gene_col) {
    counts_raw <- read_table_auto(counts_path)
    if (ncol(counts_raw) < 2) {
        stop("Counts matrix must have at least 2 columns")
    }
    gene_idx <- resolve_column(
        counts_raw,
        gene_col,
        fallback_names = c("gene_id", "gene", "Geneid"),
        fallback_index = 1,
        column_label = "gene column"
    )
    rownames(counts_raw) <- counts_raw[[gene_idx]]
    counts_raw[, gene_idx] <- NULL
    counts <- as.matrix(counts_raw)
    storage.mode(counts) <- "numeric"
    counts
}

read_expression_matrix <- function(path, gene_col) {
    data_raw <- read_table_auto(path)
    if (ncol(data_raw) < 2) {
        stop("Expression matrix must have at least 2 columns")
    }
    gene_idx <- resolve_column(
        data_raw,
        gene_col,
        fallback_names = c("gene_id", "gene", "Geneid"),
        fallback_index = 1,
        column_label = "gene column"
    )
    rownames(data_raw) <- data_raw[[gene_idx]]
    data_raw[, gene_idx] <- NULL
    data_mat <- as.matrix(data_raw)
    storage.mode(data_mat) <- "numeric"
    data_mat
}

order_counts_by_samples <- function(counts, sample_ids) {
    missing_samples <- setdiff(sample_ids, colnames(counts))
    if (length(missing_samples) > 0) {
        stop("Missing samples in counts matrix: ", paste(missing_samples, collapse = ", "))
    }
    counts[, sample_ids, drop = FALSE]
}

write_counts_matrix <- function(counts, outdir, filename = "counts_matrix.tsv") {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    counts_out <- data.frame(gene_id = rownames(counts), counts, check.names = FALSE)
    write.table(counts_out, file.path(outdir, filename),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
}

write_matrix <- function(matrix, outdir, filename) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    output <- data.frame(gene_id = rownames(matrix), matrix, check.names = FALSE)
    write.table(output, file.path(outdir, filename),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
}

filter_counts <- function(counts, min_counts, min_samples) {
    keep <- rowSums(counts >= min_counts) >= min_samples
    counts[keep, , drop = FALSE]
}

calc_log2_cpm <- function(counts) {
    lib_sizes <- colSums(counts)
    cpm <- sweep(counts, 2, lib_sizes, FUN = "/") * 1e6
    log2(cpm + 1)
}
