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

read_panorgs <- function(file) {
  unique(scan(file, what = "character", quiet = TRUE))
}

norm_chr <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "null", "NULL", "[Not Available]")] <- NA_character_
  x
}

clean_tcga_patient <- function(x) {
  x <- toupper(gsub("\\.", "-", norm_chr(x)))
  x <- ifelse(!is.na(x) & nchar(x) >= 12L, substr(x, 1L, 12L), x)
  x
}

coalesce_text <- function(df, candidates) {
  out <- rep(NA_character_, nrow(df))
  for (col in candidates) {
    if (!col %in% colnames(df)) {
      next
    }
    value <- norm_chr(df[[col]])
    idx <- is.na(out) & !is.na(value)
    out[idx] <- value[idx]
  }
  out
}

coalesce_numeric <- function(df, candidates) {
  out <- rep(NA_real_, nrow(df))
  for (col in candidates) {
    if (!col %in% colnames(df)) {
      next
    }
    value <- suppressWarnings(as.numeric(df[[col]]))
    idx <- is.na(out) & !is.na(value)
    out[idx] <- value[idx]
  }
  out
}

collapse_stage <- function(x) {
  value <- toupper(norm_chr(x))
  value <- gsub("^STAGE ", "", value)
  value <- gsub("[ABC]$", "", value)
  value[value %in% c("STAGE X", "X")] <- NA_character_
  value
}

pick_existing_columns <- function(df, candidates) {
  candidates[candidates %in% colnames(df)]
}

make_sample_patient_columns <- function(meta) {
  meta$sample <- rownames(meta)
  patient_col <- c("patient", "cases", "bcr_patient_barcode", "submitter_id")
  if (any(patient_col %in% colnames(meta))) {
    meta$patient <- clean_tcga_patient(coalesce_text(meta, patient_col))
  } else {
    meta$patient <- clean_tcga_patient(meta$sample)
  }
  meta
}

deduplicate_patients <- function(count, meta) {
  keep <- !duplicated(meta$patient)
  meta <- meta[keep, , drop = FALSE]
  rownames(meta) <- meta$sample
  count <- count[, meta$sample, drop = FALSE]
  stopifnot(identical(colnames(count), rownames(meta)))
  list(count = count, meta = meta)
}

summarize_candidate_column <- function(meta, col) {
  value <- norm_chr(meta[[col]])
  non_na <- sum(!is.na(value))
  uniq <- sort(unique(value[!is.na(value)]))
  data.frame(
    column = col,
    non_missing = non_na,
    coverage = round(non_na / nrow(meta), 4),
    n_levels = length(uniq),
    levels = paste(uniq, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

choose_molecular_entity_column <- function(meta, candidate_columns = character(), preference = "exact5_first") {
  generic_columns <- grep(
    "subtype|cluster|icluster|cms|cris|classifier|pam50",
    colnames(meta),
    ignore.case = TRUE,
    value = TRUE
  )
  generic_columns <- setdiff(
    generic_columns,
    grep("time|days|status|age|sample|patient|barcode|gene|score", generic_columns, ignore.case = TRUE, value = TRUE)
  )

  cols_to_check <- unique(c(candidate_columns, generic_columns))
  cols_to_check <- cols_to_check[cols_to_check %in% colnames(meta)]

  if (length(cols_to_check) == 0L) {
    return(list(column = NULL, summary = data.frame()))
  }

  summary_df <- do.call(rbind, lapply(cols_to_check, function(col) summarize_candidate_column(meta, col)))
  summary_df <- summary_df[summary_df$non_missing > 0, , drop = FALSE]

  if (nrow(summary_df) == 0L) {
    return(list(column = NULL, summary = summary_df))
  }

  chosen <- NULL
  if (identical(preference, "exact5_first")) {
    five_level <- summary_df[summary_df$n_levels == 5 & summary_df$coverage >= 0.2, , drop = FALSE]
    if (nrow(five_level) > 0L) {
      chosen <- five_level$column[which.max(five_level$coverage)]
    }
  }
  if (is.null(chosen)) {
    usable <- summary_df[summary_df$n_levels >= 2 & summary_df$n_levels <= 8 & summary_df$coverage >= 0.2, , drop = FALSE]
    if (nrow(usable) > 0L) {
      chosen <- usable$column[which.max(usable$coverage)]
    }
  }

  list(column = chosen, summary = summary_df[order(-summary_df$coverage, summary_df$n_levels), , drop = FALSE])
}

extract_survival_table <- function(meta) {
  os_time <- ifelse(
    !is.na(coalesce_numeric(meta, c("paper_CLIN.days_to_death", "days_to_death"))),
    coalesce_numeric(meta, c("paper_CLIN.days_to_death", "days_to_death")),
    coalesce_numeric(meta, c("paper_CLIN.days_to_last_followup", "days_to_last_follow_up", "days_to_last_followup"))
  )

  vital_status <- toupper(coalesce_text(meta, c("paper_CLIN.vital_status", "vital_status")))
  os_event <- ifelse(vital_status %in% c("DEAD"), 1L, ifelse(is.na(vital_status), NA_integer_, 0L))

  pfi_time <- coalesce_numeric(meta, c("paper_PFI.time", "PFI.time", "paper_PFS.time", "PFS.time"))
  pfi_event <- coalesce_numeric(meta, c("paper_PFI", "PFI", "paper_PFS", "PFS"))

  if (all(is.na(pfi_time))) {
    pfi_time <- coalesce_numeric(meta, c("paper_CLIN.days_to_last_followup", "days_to_last_follow_up", "days_to_last_followup"))
  }

  if (all(is.na(pfi_event))) {
    tumor_status2 <- coalesce_text(meta, c("paper_CLIN.tumorStatus2", "tumorStatus2"))
    pfi_event <- dplyr::case_when(
      tumor_status2 %in% c("Dead_wTumor", "Alive_wTumor") ~ 1,
      tumor_status2 %in% c("Alive_woTumor") ~ 0,
      TRUE ~ NA_real_
    )
  }

  stage <- collapse_stage(coalesce_text(meta, c("paper_CLIN.clinStage", "paper_clinical_stage", "ajcc_pathologic_stage", "tumor_stage")))
  t_stage <- collapse_stage(coalesce_text(meta, c("ajcc_pathologic_t", "paper_T_stage", "pathologic_T", "pathologic_t")))
  n_stage <- collapse_stage(coalesce_text(meta, c("ajcc_pathologic_n", "paper_N_stage", "pathologic_N", "pathologic_n")))
  m_stage <- collapse_stage(coalesce_text(meta, c("ajcc_pathologic_m", "paper_M_stage", "pathologic_M", "pathologic_m")))

  gleason_score <- coalesce_text(meta, c("paper_gleason_score", "gleason_score", "paper_CLIN.gleason_score"))
  if (all(is.na(gleason_score))) {
    pri <- coalesce_numeric(meta, c("paper_primary_gleason_grade", "primary_gleason_grade"))
    sec <- coalesce_numeric(meta, c("paper_secondary_gleason_grade", "secondary_gleason_grade"))
    gleason_score <- ifelse(!is.na(pri) & !is.na(sec), paste0(pri, "+", sec), NA_character_)
  }

  out <- data.frame(
    sample = meta$sample,
    patient = meta$patient,
    age = coalesce_numeric(meta, c("age_at_index", "paper_age_at_initial_pathologic_diagnosis", "age_at_diagnosis")),
    sex = coalesce_text(meta, c("gender", "sex")),
    tumor_grade = coalesce_text(meta, c("tumor_grade", "paper_tumor_grade", "neoplasm_histologic_grade")),
    stage = stage,
    T_stage = t_stage,
    N_stage = n_stage,
    M_stage = m_stage,
    gleason_score = gleason_score,
    OS_time = suppressWarnings(as.numeric(os_time)),
    OS_event = suppressWarnings(as.integer(os_event)),
    PFI_time = suppressWarnings(as.numeric(pfi_time)),
    PFI_event = suppressWarnings(as.integer(pfi_event)),
    stringsAsFactors = FALSE
  )

  out
}
