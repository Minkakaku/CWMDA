required_analysis_packages <- function() {
  c(
    "dplyr", "readr", "tibble", "ggplot2", "forcats", "broom",
    "survival", "survminer", "pheatmap", "GSVA", "ConsensusClusterPlus"
  )
}

compute_panoscore <- function(expr_matrix, genes) {
  genes_use <- intersect(genes, rownames(expr_matrix))
  if (length(genes_use) < 5L) {
    stop("Too few PANoRG genes were found in the expression matrix.")
  }

  gene_sets <- list(PANoScore = genes_use)
  if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    score_matrix <- tryCatch(
      {
        param <- GSVA::ssgseaParam(
          exprData = as.matrix(expr_matrix),
          geneSets = gene_sets,
          minSize = 1,
          normalize = TRUE
        )
        GSVA::gsva(param, verbose = FALSE)
      },
      error = function(e) {
        GSVA::gsva(
          as.matrix(expr_matrix),
          gene_sets,
          method = "ssgsea",
          kcdf = "Gaussian",
          abs.ranking = TRUE,
          ssgsea.norm = TRUE,
          verbose = FALSE
        )
      }
    )
  } else {
    score_matrix <- GSVA::gsva(
      as.matrix(expr_matrix),
      gene_sets,
      method = "ssgsea",
      kcdf = "Gaussian",
      abs.ranking = TRUE,
      ssgsea.norm = TRUE,
      verbose = FALSE
    )
  }

  score <- as.numeric(score_matrix[1, ])
  names(score) <- colnames(expr_matrix)
  score
}

run_consensus_cluster <- function(signature_matrix,
                                  outdir,
                                  chosen_k = 2L,
                                  max_k = 6L,
                                  reps = 500L,
                                  seed = 123) {
  ensure_dir(outdir)
  result <- ConsensusClusterPlus::ConsensusClusterPlus(
    d = as.matrix(signature_matrix),
    maxK = max_k,
    reps = reps,
    pItem = 0.8,
    pFeature = 1,
    clusterAlg = "hc",
    distance = "pearson",
    seed = seed,
    title = outdir,
    plot = "png",
    verbose = FALSE
  )
  cluster_raw <- result[[chosen_k]]$consensusClass
  cluster_raw <- cluster_raw[colnames(signature_matrix)]
  list(result = result, cluster = cluster_raw)
}

relabel_clusters_by_score <- function(cluster_raw, panoscore) {
  cluster_raw <- as.character(cluster_raw)
  df <- tibble::tibble(
    sample_id = names(cluster_raw),
    raw_cluster = cluster_raw,
    PANoScore = as.numeric(panoscore[names(cluster_raw)])
  ) |>
    dplyr::group_by(raw_cluster) |>
    dplyr::summarise(mean_score = mean(PANoScore, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(mean_score)

  mapping <- stats::setNames(paste0("C", seq_len(nrow(df))), df$raw_cluster)
  relabeled <- unname(mapping[cluster_raw])
  names(relabeled) <- names(cluster_raw)
  relabeled
}

build_analysis_df <- function(clin_df,
                              panoscore,
                              cluster_label,
                              time_col = "OS_time",
                              event_col = "OS_event") {
  df <- clin_df |>
    dplyr::mutate(
      PANoScore = as.numeric(panoscore[patient_id]),
      cluster = factor(cluster_label[patient_id], levels = c("C1", "C2"))
    )

  score_cut <- stats::median(df$PANoScore, na.rm = TRUE)
  df$PANoScore_group <- factor(
    ifelse(df$PANoScore >= score_cut, "High", "Low"),
    levels = c("Low", "High")
  )

  df[[time_col]] <- suppressWarnings(as.numeric(df[[time_col]]))
  df[[event_col]] <- parse_event_vector(df[[event_col]])

  if ("age" %in% colnames(df)) {
    df$age <- suppressWarnings(as.numeric(df$age))
  }
  if ("sex" %in% colnames(df)) {
    df$sex <- factor(df$sex)
  }
  for (feature in c("stage", "T_stage", "N_stage", "M_stage", "gleason_score", "molecular_entity")) {
    if (feature %in% colnames(df)) {
      df[[feature]] <- factor(df[[feature]])
    }
  }
  df
}

plot_signature_heatmap <- function(signature_matrix, annotation_df, out_file) {
  annotation_col <- annotation_df |>
    dplyr::select(cluster, PANoScore_group) |>
    as.data.frame()
  rownames(annotation_col) <- annotation_df$patient_id
  pheatmap::pheatmap(
    mat = scale_rows(signature_matrix[, annotation_df$patient_id, drop = FALSE]),
    annotation_col = annotation_col,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    border_color = NA,
    fontsize_row = 9,
    filename = out_file,
    width = 10,
    height = 7
  )
}

plot_pca_clusters <- function(signature_matrix, annotation_df, out_file) {
  pca <- stats::prcomp(t(signature_matrix), center = TRUE, scale. = TRUE)
  pca_df <- tibble::tibble(
    patient_id = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  ) |>
    dplyr::left_join(
      annotation_df |>
        dplyr::select(patient_id, cluster, PANoScore),
      by = "patient_id"
    )

  p <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2, color = cluster)) +
    ggplot2::geom_point(size = 2.5, alpha = 0.85) +
    ggplot2::stat_ellipse(linewidth = 0.7, linetype = 2, show.legend = FALSE) +
    ggplot2::labs(title = "PCA of PANoRG expression", x = "PC1", y = "PC2", color = "Cluster") +
    ggplot2::theme_bw(base_size = 12)

  save_ggplot(p, out_file, width = 7, height = 5)
}

plot_score_by_cluster <- function(annotation_df, out_file) {
  p <- ggplot2::ggplot(annotation_df, ggplot2::aes(cluster, PANoScore, fill = cluster)) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
    ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9) +
    ggplot2::geom_jitter(width = 0.12, size = 1.2, alpha = 0.55) +
    ggplot2::labs(title = "PANoScore by cluster", x = NULL, y = "PANoScore") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "none")
  save_ggplot(p, out_file, width = 5, height = 4.5)
}

save_km_plot <- function(df, group_col, time_col, event_col, out_file, title_text) {
  fit_df <- df |>
    dplyr::filter(!is.na(.data[[group_col]]), !is.na(.data[[time_col]]), !is.na(.data[[event_col]]))

  if (nrow(fit_df) < 20L || length(unique(fit_df[[group_col]])) < 2L) {
    return(invisible(NULL))
  }

  fit <- survival::survfit(
    stats::as.formula(paste0("survival::Surv(", time_col, ", ", event_col, ") ~ ", group_col)),
    data = fit_df
  )

  km_plot <- survminer::ggsurvplot(
    fit = fit,
    data = fit_df,
    pval = TRUE,
    risk.table = TRUE,
    ggtheme = ggplot2::theme_bw(base_size = 11),
    title = title_text
  )

  grDevices::pdf(out_file, width = 7.2, height = 6.5)
  print(km_plot)
  grDevices::dev.off()
  invisible(fit)
}

is_valid_cox_variable <- function(x) {
  if (all(is.na(x))) {
    return(FALSE)
  }
  if (is.numeric(x)) {
    return(sum(!is.na(x)) >= 20L && stats::sd(x, na.rm = TRUE) > 0)
  }
  value <- unique(na.omit(as.character(x)))
  length(value) >= 2L
}

extract_cox_table <- function(fit, model_name) {
  coef_df <- as.data.frame(summary(fit)$coefficients)
  ci_df <- as.data.frame(summary(fit)$conf.int)
  if (nrow(coef_df) == 0L) {
    return(tibble::tibble())
  }
  tibble::tibble(
    model = model_name,
    term = rownames(coef_df),
    HR = ci_df$`exp(coef)`,
    conf.low = ci_df$`lower .95`,
    conf.high = ci_df$`upper .95`,
    p.value = coef_df$`Pr(>|z|)`
  )
}

plot_cox_forest <- function(cox_table, out_file, title_text) {
  plot_df <- cox_table |>
    dplyr::filter(is.finite(HR), is.finite(conf.low), is.finite(conf.high), conf.low > 0) |>
    dplyr::mutate(term = factor(term, levels = rev(unique(term))))

  if (nrow(plot_df) == 0L) {
    return(invisible(NULL))
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(term, HR, ymin = conf.low, ymax = conf.high)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
    ggplot2::geom_pointrange() +
    ggplot2::coord_flip() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(title = title_text, x = NULL, y = "Hazard ratio (log scale)") +
    ggplot2::theme_bw(base_size = 11)

  save_ggplot(p, out_file, width = 7, height = max(4.5, 0.38 * nrow(plot_df)))
}

run_cox_models <- function(df, time_col, event_col, clinical_features, outdir) {
  cox_df <- df |>
    dplyr::filter(!is.na(.data[[time_col]]), !is.na(.data[[event_col]]))

  if (nrow(cox_df) < 30L || sum(cox_df[[event_col]], na.rm = TRUE) < 5L) {
    return(invisible(NULL))
  }

  uni_vars <- unique(c("PANoScore", "cluster", clinical_features))
  uni_tables <- lapply(uni_vars, function(var_name) {
    if (!var_name %in% colnames(cox_df) || !is_valid_cox_variable(cox_df[[var_name]])) {
      return(NULL)
    }
    formula_text <- paste0("survival::Surv(", time_col, ", ", event_col, ") ~ ", var_name)
    fit <- tryCatch(
      survival::coxph(stats::as.formula(formula_text), data = cox_df),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(NULL)
    }
    extract_cox_table(fit, model_name = paste0("univariate:", var_name))
  })

  uni_table <- dplyr::bind_rows(uni_tables)
  if (nrow(uni_table) > 0L) {
    readr::write_tsv(uni_table, file.path(outdir, "cox_univariate.tsv"))
    plot_cox_forest(uni_table, file.path(outdir, "cox_univariate_forest.png"), "Univariate Cox")
  }

  multi_vars <- unique(c("PANoScore", clinical_features))
  multi_vars <- multi_vars[multi_vars %in% colnames(cox_df)]
  multi_vars <- multi_vars[vapply(multi_vars, function(x) is_valid_cox_variable(cox_df[[x]]), logical(1))]

  if (length(multi_vars) == 0L) {
    return(invisible(list(univariate = uni_table, multivariate = tibble::tibble())))
  }

  multi_formula <- paste0(
    "survival::Surv(", time_col, ", ", event_col, ") ~ ",
    paste(multi_vars, collapse = " + ")
  )
  multi_fit <- tryCatch(
    survival::coxph(stats::as.formula(multi_formula), data = cox_df),
    error = function(e) NULL
  )
  if (!is.null(multi_fit)) {
    multi_table <- extract_cox_table(multi_fit, model_name = "multivariate")
    readr::write_tsv(multi_table, file.path(outdir, "cox_multivariate.tsv"))
    if (nrow(multi_table) > 0L) {
      plot_cox_forest(multi_table, file.path(outdir, "cox_multivariate_forest.png"), "Multivariate Cox")
    }
  } else {
    multi_table <- tibble::tibble()
  }

  invisible(list(univariate = uni_table, multivariate = multi_table))
}

run_clinical_associations <- function(df, clinical_features, outdir) {
  results <- lapply(clinical_features, function(feature) {
    if (!feature %in% colnames(df)) {
      return(NULL)
    }
    sub_df <- df |>
      dplyr::filter(!is.na(PANoScore), !is.na(.data[[feature]]))

    if (nrow(sub_df) < 10L) {
      return(NULL)
    }

    x <- sub_df[[feature]]
    y <- sub_df$PANoScore

    if (is.numeric(x)) {
      test <- suppressWarnings(stats::cor.test(x, y, method = "spearman"))
      p <- ggplot2::ggplot(sub_df, ggplot2::aes(.data[[feature]], PANoScore)) +
        ggplot2::geom_point(alpha = 0.75, size = 1.8) +
        ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "#B23A48") +
        ggplot2::labs(title = paste("PANoScore vs", feature), x = feature, y = "PANoScore") +
        ggplot2::theme_bw(base_size = 12)
      save_ggplot(
        p,
        file.path(outdir, paste0("PANoScore_by_", sanitize_feature_name(feature), ".png")),
        width = 5.2,
        height = 4.5
      )
      return(tibble::tibble(
        feature = feature,
        variable_type = "numeric",
        n = nrow(sub_df),
        test = "Spearman",
        statistic = unname(test$estimate),
        p.value = test$p.value
      ))
    }

    x <- forcats::fct_drop(factor(x))
    if (nlevels(x) < 2L) {
      return(NULL)
    }
    sub_df$.feature_value <- x
    test <- if (nlevels(x) == 2L) {
      stats::wilcox.test(PANoScore ~ .feature_value, data = sub_df)
    } else {
      stats::kruskal.test(PANoScore ~ .feature_value, data = sub_df)
    }

    p <- ggplot2::ggplot(sub_df, ggplot2::aes(.feature_value, PANoScore, fill = .feature_value)) +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
      ggplot2::geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.9) +
      ggplot2::geom_jitter(width = 0.12, size = 1.1, alpha = 0.45) +
      ggplot2::labs(title = paste("PANoScore by", feature), x = feature, y = "PANoScore") +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(legend.position = "none")
    save_ggplot(
      p,
      file.path(outdir, paste0("PANoScore_by_", sanitize_feature_name(feature), ".png")),
      width = max(5.5, 0.85 * nlevels(x) + 3.5),
      height = 4.5
    )

    tibble::tibble(
      feature = feature,
      variable_type = "categorical",
      n = nrow(sub_df),
      test = if (nlevels(x) == 2L) "Wilcoxon" else "Kruskal-Wallis",
      statistic = unname(test$statistic),
      p.value = test$p.value
    )
  })

  association_table <- dplyr::bind_rows(results)
  if (nrow(association_table) > 0L) {
    readr::write_tsv(association_table, file.path(outdir, "clinical_associations.tsv"))
  }
  invisible(association_table)
}

plot_molecular_entity_distribution <- function(df, out_file) {
  if (!"molecular_entity" %in% colnames(df)) {
    return(invisible(NULL))
  }
  plot_df <- df |>
    dplyr::filter(!is.na(PANoScore), !is.na(molecular_entity))

  if (nrow(plot_df) < 10L || length(unique(plot_df$molecular_entity)) < 2L) {
    return(invisible(NULL))
  }

  plot_df$molecular_entity <- forcats::fct_infreq(plot_df$molecular_entity)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(molecular_entity, PANoScore, fill = molecular_entity)) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
    ggplot2::geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.95) +
    ggplot2::geom_jitter(width = 0.12, size = 1.0, alpha = 0.45) +
    ggplot2::labs(title = "PANoScore across molecular entities", x = NULL, y = "PANoScore") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "none")
  save_ggplot(p, out_file, width = max(6, 0.8 * nlevels(plot_df$molecular_entity) + 3), height = 4.8)
}

run_panorg_analysis <- function(cohort_row,
                                panorgs,
                                chosen_k = 2L,
                                max_k = 6L,
                                reps = 500L,
                                seed = 123) {
  outdir <- project_path(cohort_row$results_dir)
  ensure_dir(outdir)

  expr <- read_expression_tsv(project_path(cohort_row$expr_file), gene_col = "gene_symbol")
  expr <- apply_expr_transform(expr, transform = cohort_row$expr_transform)

  clin <- readr::read_tsv(project_path(cohort_row$clin_file), show_col_types = FALSE) |>
    dplyr::mutate(patient_id = clean_tcga_patient_id(patient_id)) |>
    dplyr::distinct(patient_id, .keep_all = TRUE)

  common_ids <- intersect(colnames(expr), clin$patient_id)
  expr <- expr[, common_ids, drop = FALSE]
  clin <- clin |>
    dplyr::filter(patient_id %in% common_ids) |>
    dplyr::slice(match(common_ids, patient_id))

  signature_genes <- intersect(panorgs, rownames(expr))
  if (length(signature_genes) < 5L) {
    stop("Fewer than 5 PANoRG genes are available in ", cohort_row$cohort_id)
  }

  signature_matrix <- expr[signature_genes, , drop = FALSE]
  panoscore <- compute_panoscore(expr, panorgs)

  cluster_fit <- run_consensus_cluster(
    signature_matrix = scale_rows(signature_matrix),
    outdir = file.path(outdir, "consensus_cluster"),
    chosen_k = chosen_k,
    max_k = max_k,
    reps = reps,
    seed = seed
  )

  cluster_label <- relabel_clusters_by_score(cluster_fit$cluster, panoscore = panoscore)
  analysis_df <- build_analysis_df(
    clin_df = clin,
    panoscore = panoscore,
    cluster_label = cluster_label,
    time_col = cohort_row$time_col,
    event_col = cohort_row$event_col
  ) |>
    dplyr::arrange(cluster, dplyr::desc(PANoScore))

  readr::write_tsv(analysis_df, file.path(outdir, "sample_annotations.tsv.gz"))
  plot_signature_heatmap(signature_matrix, analysis_df, file.path(outdir, "panorg_heatmap.png"))
  plot_pca_clusters(signature_matrix, analysis_df, file.path(outdir, "pca_clusters.png"))
  plot_score_by_cluster(analysis_df, file.path(outdir, "PANoScore_by_cluster.png"))

  save_km_plot(
    df = analysis_df,
    group_col = "cluster",
    time_col = cohort_row$time_col,
    event_col = cohort_row$event_col,
    out_file = file.path(outdir, "km_cluster.pdf"),
    title_text = paste0(cohort_row$cohort_id, ": Cluster survival")
  )

  save_km_plot(
    df = analysis_df,
    group_col = "PANoScore_group",
    time_col = cohort_row$time_col,
    event_col = cohort_row$event_col,
    out_file = file.path(outdir, "km_panoscore.pdf"),
    title_text = paste0(cohort_row$cohort_id, ": PANoScore survival")
  )

  clinical_features <- split_semicolon_field(cohort_row$clinical_features)
  run_cox_models(
    df = analysis_df,
    time_col = cohort_row$time_col,
    event_col = cohort_row$event_col,
    clinical_features = clinical_features,
    outdir = outdir
  )

  run_clinical_associations(
    df = analysis_df,
    clinical_features = clinical_features,
    outdir = outdir
  )

  plot_molecular_entity_distribution(
    df = analysis_df,
    out_file = file.path(outdir, "PANoScore_by_molecular_entity.png")
  )

  invisible(analysis_df)
}
