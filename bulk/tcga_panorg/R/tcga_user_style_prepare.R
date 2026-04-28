required_user_style_packages <- function() {
  c(
    "TCGAbiolinks", "SummarizedExperiment", "dplyr", "tibble",
    "readr", "limma", "edgeR", "org.Hs.eg.db", "AnnotationDbi"
  )
}

check_required_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0L) {
    stop("Missing packages: ", paste(missing, collapse = ", "), ". Please run scripts/00_install_packages.R first.")
  }
}

download_tcga_project_user_style <- function(project, cache_dir = project_path("gdc_cache")) {
  ensure_dir(cache_dir)
  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  TCGAbiolinks::GDCdownload(query, files.per.chunk = 50, directory = cache_dir)
  se <- TCGAbiolinks::GDCprepare(query, directory = cache_dir)
  se
}

prepare_count_meta_from_se <- function(se) {
  count <- SummarizedExperiment::assay(se)
  row_df <- as.data.frame(SummarizedExperiment::rowData(se))
  gene_col <- if ("gene_id" %in% colnames(row_df)) "gene_id" else colnames(row_df)[1]
  gene_id <- row_df[[gene_col]]
  gene_id <- sub("\\..*$", "", as.character(gene_id))
  keep_id <- !is.na(gene_id) & gene_id != ""

  count <- count[keep_id, , drop = FALSE]
  gene_id <- gene_id[keep_id]
  rownames(count) <- gene_id
  count <- rowsum(count, group = rownames(count))

  meta <- as.data.frame(SummarizedExperiment::colData(se))
  meta <- make_sample_patient_columns(meta)

  if ("shortLetterCode" %in% colnames(meta)) {
    tumor_flag <- meta$shortLetterCode == "TP"
  } else {
    tumor_flag <- substr(meta$sample, 14, 15) == "01"
  }
  tumor_flag[is.na(tumor_flag)] <- FALSE
  count <- count[, tumor_flag, drop = FALSE]
  meta <- meta[tumor_flag, , drop = FALSE]
  rownames(meta) <- meta$sample

  dedup <- deduplicate_patients(count, meta)
  count <- dedup$count
  meta <- dedup$meta

  stopifnot(identical(colnames(count), rownames(meta)))
  list(count = count, meta = meta)
}

merge_project_bundles <- function(bundle_list) {
  all_genes <- unique(unlist(lapply(bundle_list, function(x) rownames(x$count))))
  count_list <- lapply(bundle_list, function(x) {
    mat <- x$count
    idx <- match(all_genes, rownames(mat))
    out <- matrix(0, nrow = length(all_genes), ncol = ncol(mat))
    rownames(out) <- all_genes
    colnames(out) <- colnames(mat)
    valid <- !is.na(idx)
    out[valid, ] <- mat[idx[valid], , drop = FALSE]
    out
  })

  count <- do.call(cbind, count_list)
  meta <- dplyr::bind_rows(lapply(bundle_list, `[[`, "meta"))
  rownames(meta) <- meta$sample
  meta <- meta[colnames(count), , drop = FALSE]

  dedup <- deduplicate_patients(count, meta)
  count <- dedup$count
  meta <- dedup$meta

  stopifnot(identical(colnames(count), rownames(meta)))
  list(count = count, meta = meta)
}

normalize_count_to_voom <- function(count, meta) {
  dge <- edgeR::DGEList(counts = count)
  dge <- edgeR::calcNormFactors(dge)
  keep <- edgeR::filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  design0 <- stats::model.matrix(~1, data = meta)
  v <- limma::voom(dge, design0, plot = FALSE)
  expr <- v$E
  stopifnot(identical(colnames(expr), rownames(meta)))
  list(dge = dge, voom = v, expr = expr)
}

select_panorg_ensembl <- function(expr, panorgs) {
  map_df <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(panorgs),
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENSEMBL")
  )
  map_df <- unique(map_df[, c("SYMBOL", "ENSEMBL")])
  map_df$ENSEMBL <- sub("\\..*$", "", map_df$ENSEMBL)
  map_df <- map_df[!is.na(map_df$SYMBOL) & !is.na(map_df$ENSEMBL), , drop = FALSE]

  matched <- lapply(unique(panorgs), function(sym) {
    sub_df <- map_df[map_df$SYMBOL == sym, , drop = FALSE]
    hit <- sub_df$ENSEMBL[sub_df$ENSEMBL %in% rownames(expr)]
    if (length(hit) == 0L) {
      return(NULL)
    }
    data.frame(symbol = sym, ensembl = hit[[1]], stringsAsFactors = FALSE)
  })
  matched <- matched[!vapply(matched, is.null, logical(1))]
  if (length(matched) == 0L) {
    return(data.frame(symbol = character(), ensembl = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, matched)
}

prepare_one_tcga_cohort <- function(cohort_row, panorgs, candidate_cfg) {
  cohort_id <- cohort_row$cohort_id
  projects <- split_semicolon_field(cohort_row$projects)

  bundle_list <- lapply(projects, function(project) {
    se <- download_tcga_project_user_style(project)
    prepare_count_meta_from_se(se)
  })
  merged <- merge_project_bundles(bundle_list)
  normed <- normalize_count_to_voom(merged$count, merged$meta)

  meta <- merged$meta
  meta_std <- extract_survival_table(meta)
  candidate_row <- candidate_cfg[candidate_cfg$cohort_id == cohort_id, , drop = FALSE]
  candidate_cols <- if (nrow(candidate_row) == 1L) split_semicolon_field(candidate_row$candidate_columns) else character()
  chosen_entity <- choose_molecular_entity_column(
    meta = meta,
    candidate_columns = candidate_cols,
    preference = cohort_row$molecular_group_preference %||% "exact5_first"
  )

  meta_std$molecular_entity <- if (!is.null(chosen_entity$column)) norm_chr(meta[[chosen_entity$column]]) else NA_character_
  meta_std$molecular_entity_source <- chosen_entity$column %||% NA_character_

  panorg_map <- select_panorg_ensembl(normed$expr, panorgs)

  list(
    cohort_id = cohort_id,
    disease = cohort_row$disease,
    projects = projects,
    count = merged$count,
    expr = normed$expr,
    meta_raw = meta,
    meta_std = meta_std,
    panorg_map = panorg_map,
    molecular_entity_candidates = chosen_entity$summary
  )
}

write_tcga_bundle_outputs <- function(bundle, outdir) {
  ensure_dir(outdir)
  ensure_dir(file.path(outdir, "rds"))
  ensure_dir(file.path(outdir, "tables"))

  saveRDS(bundle, file.path(outdir, "rds", "tcga_bundle.rds"))
  readr::write_csv(bundle$meta_std, file.path(outdir, "tables", "meta_full.csv"))
  readr::write_csv(bundle$panorg_map, file.path(outdir, "tables", "panorg_gene_map.csv"))
  if (nrow(bundle$molecular_entity_candidates) > 0L) {
    readr::write_csv(bundle$molecular_entity_candidates, file.path(outdir, "tables", "molecular_entity_candidates.csv"))
  }
}
