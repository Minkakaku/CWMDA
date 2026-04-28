required_download_packages <- function() {
  c(
    "TCGAbiolinks", "SummarizedExperiment",
    "dplyr", "tibble", "readr", "stringr", "purrr"
  )
}

check_required_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0L) {
    stop(
      "Missing packages: ",
      paste(missing, collapse = ", "),
      ". Please run scripts/00_install_packages.R first."
    )
  }
}

extract_gene_symbols_from_se <- function(se_object) {
  row_data <- as.data.frame(SummarizedExperiment::rowData(se_object))
  gene_col <- first_present_col(
    row_data,
    c("gene_name", "external_gene_name", "gene_short_name", "Symbol", "symbol")
  )
  if (is.null(gene_col)) {
    stop("Cannot find a gene symbol column in SummarizedExperiment rowData.")
  }
  as.character(row_data[[gene_col]])
}

download_tcga_expression_project <- function(project,
                                             workflow_type = "HTSeq - FPKM-UQ",
                                             sample_types = "Primary Tumor",
                                             cache_dir = project_path("gdc_cache")) {
  ensure_dir(cache_dir)
  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = workflow_type,
    sample.type = sample_types
  )
  TCGAbiolinks::GDCdownload(query = query, directory = cache_dir)
  se_object <- TCGAbiolinks::GDCprepare(query = query, directory = cache_dir)
  expr <- SummarizedExperiment::assay(se_object)
  genes <- extract_gene_symbols_from_se(se_object)
  expr <- aggregate_expression_by_gene(expr, genes)
  expr <- deduplicate_patient_matrix(expr)
  expr
}

standardize_tcga_clinical <- function(clinical_df, project) {
  df <- tibble::as_tibble(clinical_df)
  patient_id <- clean_tcga_patient_id(
    coalesce_nonempty(
      df,
      c("submitter_id", "bcr_patient_barcode", "case_submitter_id", "patient_id")
    )
  )

  age_value <- coalesce_nonempty(df, c("age_at_diagnosis", "age_at_index", "days_to_birth"))
  age <- days_to_years(age_value)
  if ("days_to_birth" %in% colnames(df) && all(is.na(age))) {
    age <- abs(suppressWarnings(as.numeric(df$days_to_birth))) / 365.25
  }

  vital_status <- coalesce_nonempty(df, c("vital_status"))
  last_follow_up <- coalesce_nonempty(
    df,
    c("days_to_last_follow_up", "days_to_last_known_alive", "days_to_last_known_disease_status")
  )
  days_to_death <- coalesce_nonempty(df, c("days_to_death"))

  last_follow_up <- suppressWarnings(as.numeric(last_follow_up))
  days_to_death <- suppressWarnings(as.numeric(days_to_death))
  os_time <- ifelse(!is.na(days_to_death), days_to_death, last_follow_up)
  os_event <- ifelse(toupper(vital_status) == "DEAD", 1L, 0L)
  os_event[is.na(vital_status)] <- NA_integer_

  primary_gleason <- suppressWarnings(as.numeric(coalesce_nonempty(df, c("primary_gleason_grade"))))
  secondary_gleason <- suppressWarnings(as.numeric(coalesce_nonempty(df, c("secondary_gleason_grade"))))
  gleason_score <- coalesce_nonempty(df, c("gleason_score", "gleason_total", "total_gleason_score"))
  if (all(is.na(gleason_score)) && any(!is.na(primary_gleason)) && any(!is.na(secondary_gleason))) {
    gleason_score <- ifelse(
      !is.na(primary_gleason) & !is.na(secondary_gleason),
      paste0(primary_gleason, "+", secondary_gleason),
      NA_character_
    )
  }

  tibble::tibble(
    patient_id = patient_id,
    disease = dplyr::case_when(
      project %in% c("TCGA-COAD", "TCGA-READ") ~ "CRC",
      project %in% c("TCGA-PRAD") ~ "PRAD",
      TRUE ~ project
    ),
    source_project = project,
    age = age,
    sex = coalesce_nonempty(df, c("gender", "sex")),
    vital_status = vital_status,
    OS_time = os_time,
    OS_event = os_event,
    stage = collapse_stage(coalesce_nonempty(df, c("ajcc_pathologic_stage", "tumor_stage", "pathologic_stage"))),
    T_stage = collapse_stage(coalesce_nonempty(df, c("ajcc_pathologic_t", "pathologic_t", "ajcc_clinical_t", "clinical_t"))),
    N_stage = collapse_stage(coalesce_nonempty(df, c("ajcc_pathologic_n", "pathologic_n", "ajcc_clinical_n", "clinical_n"))),
    M_stage = collapse_stage(coalesce_nonempty(df, c("ajcc_pathologic_m", "pathologic_m", "ajcc_clinical_m", "clinical_m"))),
    gleason_score = gleason_score
  ) |>
    dplyr::filter(!is.na(patient_id) & patient_id != "") |>
    dplyr::distinct(patient_id, .keep_all = TRUE)
}

fetch_tcga_subtypes_table <- function(tumor_codes = character()) {
  subtype_df <- tryCatch(
    TCGAbiolinks::PanCancerAtlas_subtypes(),
    error = function(e) NULL
  )

  if (is.null(subtype_df) || nrow(subtype_df) == 0L) {
    subtype_list <- lapply(unique(tumor_codes), function(code) {
      tryCatch(TCGAbiolinks::TCGAquery_subtype(tumor = code), error = function(e) NULL)
    })
    subtype_df <- dplyr::bind_rows(subtype_list)
  }

  if (is.null(subtype_df) || nrow(subtype_df) == 0L) {
    return(tibble::tibble())
  }

  subtype_df <- tibble::as_tibble(subtype_df)
  barcode_col <- detect_barcode_column(subtype_df)
  subtype_col <- pick_subtype_column(subtype_df, barcode_col = barcode_col)

  if (is.null(barcode_col) || is.null(subtype_col)) {
    return(tibble::tibble())
  }

  tibble::tibble(
    patient_id = clean_tcga_patient_id(subtype_df[[barcode_col]]),
    molecular_entity = trimws(as.character(subtype_df[[subtype_col]])),
    molecular_entity_source = subtype_col
  ) |>
    dplyr::mutate(molecular_entity = dplyr::na_if(molecular_entity, "")) |>
    dplyr::filter(!is.na(patient_id) & patient_id != "") |>
    dplyr::distinct(patient_id, .keep_all = TRUE)
}

read_custom_molecular_entities <- function(custom_file, disease) {
  if (!file.exists(custom_file)) {
    return(tibble::tibble())
  }
  custom_df <- readr::read_csv(custom_file, show_col_types = FALSE)
  required_cols <- c("disease", "patient_id", "molecular_entity")
  if (!all(required_cols %in% colnames(custom_df))) {
    stop(
      "Custom molecular entity file must contain columns: ",
      paste(required_cols, collapse = ", ")
    )
  }
  custom_df |>
    dplyr::filter(.data$disease == disease) |>
    dplyr::transmute(
      patient_id = clean_tcga_patient_id(.data$patient_id),
      molecular_entity_custom = .data$molecular_entity
    ) |>
    dplyr::distinct(patient_id, .keep_all = TRUE)
}

build_tcga_cohort_bundle <- function(cohort_row,
                                     custom_file = project_path("config", "custom_molecular_entities.csv")) {
  projects <- split_semicolon_field(cohort_row$projects)
  tumor_codes <- split_semicolon_field(cohort_row$tumor_codes)
  sample_types <- split_semicolon_field(cohort_row$sample_types)
  workflow_type <- cohort_row$workflow_type

  message("Downloading expression for ", cohort_row$cohort_id, " ...")
  expr_list <- lapply(projects, function(project) {
    download_tcga_expression_project(
      project = project,
      workflow_type = workflow_type,
      sample_types = sample_types
    )
  })

  message("Downloading clinical tables for ", cohort_row$cohort_id, " ...")
  clin_list <- lapply(projects, function(project) {
    clinical_df <- TCGAbiolinks::GDCquery_clinic(project = project, type = "clinical")
    standardize_tcga_clinical(clinical_df, project = project)
  })

  expr <- merge_expression_matrices(expr_list)
  expr <- deduplicate_patient_matrix(expr)
  clin <- dplyr::bind_rows(clin_list) |>
    dplyr::distinct(patient_id, .keep_all = TRUE)

  subtypes <- fetch_tcga_subtypes_table(tumor_codes = tumor_codes)
  if (nrow(subtypes) > 0L) {
    clin <- dplyr::left_join(clin, subtypes, by = "patient_id")
  } else {
    clin$molecular_entity <- NA_character_
    clin$molecular_entity_source <- NA_character_
  }

  custom_entity <- read_custom_molecular_entities(custom_file = custom_file, disease = cohort_row$disease)
  if (nrow(custom_entity) > 0L) {
    clin <- dplyr::left_join(clin, custom_entity, by = "patient_id") |>
      dplyr::mutate(
        molecular_entity_source = ifelse(
          !is.na(.data$molecular_entity_custom),
          "custom_file",
          .data$molecular_entity_source
        ),
        molecular_entity = dplyr::coalesce(.data$molecular_entity_custom, .data$molecular_entity)
      ) |>
      dplyr::select(-.data$molecular_entity_custom)
  }

  common_ids <- intersect(colnames(expr), clin$patient_id)
  expr <- expr[, common_ids, drop = FALSE]
  clin <- clin |>
    dplyr::filter(patient_id %in% common_ids) |>
    dplyr::slice(match(common_ids, patient_id))

  ensure_dir(dirname(project_path(cohort_row$expr_file)))
  ensure_dir(dirname(project_path(cohort_row$clin_file)))
  write_expression_tsv(expr, project_path(cohort_row$expr_file))
  readr::write_tsv(clin, project_path(cohort_row$clin_file))

  invisible(list(expr = expr, clin = clin))
}
