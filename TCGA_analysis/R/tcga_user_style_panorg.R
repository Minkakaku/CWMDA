required_panorg_packages <- function() {
  c(
    "dplyr", "readr", "tibble", "ggplot2", "survival", "survminer",
    "GSVA", "ConsensusClusterPlus", "pheatmap", "forcats"
  )
}

read_panorgs <- function(file) {
  unique(scan(file, what = "character", quiet = TRUE))
}

compute_panoscore_user_style <- function(expr, panorg_map) {
  if (nrow(panorg_map) < 5L) {
    stop("Too few PANoRG genes could be mapped from SYMBOL to ENSEMBL.")
  }

  gene_ids <- unique(panorg_map$ensembl)
  gene_set <- list(PANoScore = gene_ids)

  score_matrix <- tryCatch(
    {
      param <- GSVA::ssgseaParam(
        exprData = as.matrix(expr),
        geneSets = gene_set,
        minSize = 1,
        normalize = TRUE
      )
      GSVA::gsva(param, verbose = FALSE)
    },
    error = function(e) {
      GSVA::gsva(
        as.matrix(expr),
        gene_set,
        method = "ssgsea",
        kcdf = "Gaussian",
        abs.ranking = TRUE,
        ssgsea.norm = TRUE,
        verbose = FALSE
      )
    }
  )

  score <- as.numeric(score_matrix[1, ])
  names(score) <- colnames(expr)
  score
}

make_signature_expr <- function(expr, panorg_map) {
  sig <- expr[panorg_map$ensembl, , drop = FALSE]
  rownames(sig) <- panorg_map$symbol
  sig
}

scale_rows <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z
}

run_ccp_two_cluster <- function(sig_expr, outdir, seed = 123) {
  ensure_dir(outdir)
  res <- ConsensusClusterPlus::ConsensusClusterPlus(
    d = as.matrix(scale_rows(sig_expr)),
    maxK = 6,
    reps = 500,
    pItem = 0.8,
    pFeature = 1,
    clusterAlg = "hc",
    distance = "pearson",
    seed = seed,
    title = outdir,
    plot = "pdf",
    verbose = FALSE
  )
  cls <- res[[2]]$consensusClass
  cls <- cls[colnames(sig_expr)]
  names(cls) <- colnames(sig_expr)
  cls
}

rename_clusters <- function(cluster_raw, panoscore) {
  score_df <- data.frame(
    sample = names(cluster_raw),
    cluster_raw = as.character(cluster_raw),
    PANoScore = as.numeric(panoscore[names(cluster_raw)]),
    stringsAsFactors = FALSE
  )
  ord <- score_df |>
    dplyr::group_by(cluster_raw) |>
    dplyr::summarise(mean_score = mean(PANoScore, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(mean_score)
  mapping <- setNames(paste0("C", seq_len(nrow(ord))), ord$cluster_raw)
  out <- unname(mapping[as.character(cluster_raw)])
  names(out) <- names(cluster_raw)
  out
}

build_annotation_df <- function(bundle, panoscore, cluster_label) {
  df <- bundle$meta_std
  df$PANoScore <- as.numeric(panoscore[df$sample])
  df$cluster <- factor(cluster_label[df$sample], levels = c("C1", "C2"))
  cutoff <- stats::median(df$PANoScore, na.rm = TRUE)
  df$PANoScore_group <- factor(ifelse(df$PANoScore >= cutoff, "High", "Low"), levels = c("Low", "High"))
  df
}

plot_signature_heatmap_user_style <- function(sig_expr, ann_df, file) {
  annotation_col <- ann_df[, c("cluster", "PANoScore_group"), drop = FALSE]
  rownames(annotation_col) <- ann_df$sample
  pheatmap::pheatmap(
    mat = scale_rows(sig_expr[, ann_df$sample, drop = FALSE]),
    annotation_col = annotation_col,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    border_color = NA,
    filename = file,
    width = 10,
    height = 7
  )
}

plot_pca_user_style <- function(sig_expr, ann_df, file) {
  pca <- stats::prcomp(t(sig_expr), center = TRUE, scale. = TRUE)
  pca_df <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )
  pca_df <- dplyr::left_join(pca_df, ann_df[, c("sample", "cluster")], by = "sample")

  p <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2, color = cluster)) +
    ggplot2::geom_point(size = 2.5, alpha = 0.85) +
    ggplot2::stat_ellipse(linewidth = 0.8, linetype = 2, show.legend = FALSE) +
    ggplot2::labs(title = "PCA of PANoRG expression", x = "PC1", y = "PC2", color = "Cluster") +
    ggplot2::theme_bw(base_size = 12)
  ggplot2::ggsave(file, p, width = 6, height = 5)
}

plot_panoscore_by_group <- function(ann_df, group_col, file, title_text) {
  p <- ggplot2::ggplot(ann_df, ggplot2::aes(.data[[group_col]], PANoScore, fill = .data[[group_col]])) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.15, size = 1.4, alpha = 0.6) +
    ggplot2::labs(x = "", y = "PANoScore", title = title_text) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "none")
  ggplot2::ggsave(file, p, width = 4.5, height = 4)
}

fit_and_save_km <- function(df, time_col, event_col, group_col, file, title_text) {
  sub_df <- df[!is.na(df[[time_col]]) & !is.na(df[[event_col]]) & !is.na(df[[group_col]]), , drop = FALSE]
  if (nrow(sub_df) < 10L || sum(sub_df[[event_col]], na.rm = TRUE) < 5L || length(unique(sub_df[[group_col]])) < 2L) {
    return(NULL)
  }
  sub_df$.time <- as.numeric(sub_df[[time_col]])
  sub_df$.event <- as.integer(sub_df[[event_col]])
  sub_df$.group <- factor(sub_df[[group_col]])
  fit <- survival::survfit(survival::Surv(.time, .event) ~ .group, data = sub_df)
  p <- survminer::ggsurvplot(
    fit,
    data = sub_df,
    pval = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    title = title_text
  )
  ggplot2::ggsave(file, p$plot, width = 6, height = 5)
  fit
}

fit_and_save_cox <- function(df, time_col, event_col, feature_terms, prefix, outdir) {
  sub_df <- df[!is.na(df[[time_col]]) & !is.na(df[[event_col]]), , drop = FALSE]
  if (nrow(sub_df) < 10L || sum(sub_df[[event_col]], na.rm = TRUE) < 5L) {
    return(invisible(NULL))
  }

  sub_df$.time <- as.numeric(sub_df[[time_col]])
  sub_df$.event <- as.integer(sub_df[[event_col]])

  uni_formula <- survival::Surv(.time, .event) ~ PANoScore_group
  uni_fit <- tryCatch(survival::coxph(uni_formula, data = sub_df), error = function(e) NULL)
  if (!is.null(uni_fit)) {
    uni_tbl <- extract_cox_result_table(uni_fit)
    readr::write_csv(uni_tbl, file.path(outdir, "tables", paste0("cox_", prefix, "_univariate.csv")))
    save_cox_forest_plot(
      uni_tbl,
      file.path(outdir, "figures", paste0("Forest_", prefix, "_univariate.pdf")),
      paste(prefix, "Univariate")
    )
    zph_uni <- tryCatch(survival::cox.zph(uni_fit), error = function(e) NULL)
    if (!is.null(zph_uni)) {
      capture.output(zph_uni, file = file.path(outdir, "tables", paste0("cox_", prefix, "_PHtest.txt")))
    }
  }

  multi_terms <- c("PANoScore_group", feature_terms)
  multi_terms <- multi_terms[multi_terms %in% colnames(sub_df)]
  multi_terms <- multi_terms[vapply(multi_terms, function(x) {
    val <- sub_df[[x]]
    if (is.numeric(val)) {
      return(sum(!is.na(val)) >= 10L && stats::sd(val, na.rm = TRUE) > 0)
    }
    length(unique(stats::na.omit(val))) >= 2L
  }, logical(1))]

  if (length(multi_terms) >= 2L) {
    multi_formula <- stats::as.formula(
      paste0("survival::Surv(.time, .event) ~ ", paste(multi_terms, collapse = " + "))
    )
    multi_fit <- tryCatch(survival::coxph(multi_formula, data = sub_df), error = function(e) NULL)
    if (!is.null(multi_fit)) {
      multi_tbl <- extract_cox_result_table(multi_fit)
      readr::write_csv(multi_tbl, file.path(outdir, "tables", paste0("cox_", prefix, "_multivariate.csv")))
      save_cox_forest_plot(
        multi_tbl,
        file.path(outdir, "figures", paste0("Forest_", prefix, "_multivariate.pdf")),
        paste(prefix, "Multivariate")
      )
    }
  }

  invisible(NULL)
}

extract_cox_result_table <- function(fit) {
  coef_df <- as.data.frame(summary(fit)$coefficients)
  ci_df <- as.data.frame(summary(fit)$conf.int)
  if (nrow(coef_df) == 0L || nrow(ci_df) == 0L) {
    return(data.frame())
  }
  data.frame(
    term = rownames(coef_df),
    coef = coef_df[, "coef"],
    exp_coef = ci_df[, "exp(coef)"],
    se_coef = coef_df[, "se(coef)"],
    z = coef_df[, "z"],
    p_value = coef_df[, "Pr(>|z|)"],
    conf_low = ci_df[, "lower .95"],
    conf_high = ci_df[, "upper .95"],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

save_cox_forest_plot <- function(cox_tbl, file, title_text) {
  if (is.null(cox_tbl) || nrow(cox_tbl) == 0L) {
    return(invisible(NULL))
  }

  plot_df <- cox_tbl
  plot_df <- plot_df[
    is.finite(plot_df$exp_coef) &
      is.finite(plot_df$conf_low) &
      is.finite(plot_df$conf_high) &
      plot_df$exp_coef > 0 &
      plot_df$conf_low > 0 &
      plot_df$conf_high > 0,
    ,
    drop = FALSE
  ]

  if (nrow(plot_df) == 0L) {
    return(invisible(NULL))
  }

  plot_df$term <- factor(plot_df$term, levels = rev(plot_df$term))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = term, y = exp_coef, ymin = conf_low, ymax = conf_high)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
    ggplot2::geom_pointrange() +
    ggplot2::coord_flip() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "", y = "Hazard ratio", title = title_text) +
    ggplot2::theme_bw(base_size = 12)

  ggplot2::ggsave(file, p, width = 6.5, height = max(4, 0.45 * nrow(plot_df)))
}

run_clinical_association_user_style <- function(ann_df, features, outdir) {
  out_list <- list()
  for (feature in features) {
    if (!feature %in% colnames(ann_df)) {
      next
    }
    sub_df <- ann_df[!is.na(ann_df$PANoScore) & !is.na(ann_df[[feature]]), , drop = FALSE]
    if (nrow(sub_df) < 10L) {
      next
    }

    if (is.numeric(sub_df[[feature]])) {
      test <- suppressWarnings(stats::cor.test(sub_df[[feature]], sub_df$PANoScore, method = "spearman"))
      p <- ggplot2::ggplot(sub_df, ggplot2::aes(.data[[feature]], PANoScore)) +
        ggplot2::geom_point(size = 1.6, alpha = 0.7) +
        ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "#C0392B") +
        ggplot2::labs(x = feature, y = "PANoScore", title = paste("PANoScore vs", feature)) +
        ggplot2::theme_bw(base_size = 12)
      ggplot2::ggsave(file.path(outdir, "figures", paste0("PANoScore_by_", feature, ".pdf")), p, width = 5, height = 4)
      out_list[[feature]] <- data.frame(
        feature = feature,
        test = "Spearman",
        statistic = unname(test$estimate),
        p.value = test$p.value,
        stringsAsFactors = FALSE
      )
    } else {
      sub_df$.feature_value <- forcats::fct_drop(factor(sub_df[[feature]]))
      if (nlevels(sub_df$.feature_value) < 2L) {
        next
      }
      test <- if (nlevels(sub_df$.feature_value) == 2L) {
        stats::wilcox.test(PANoScore ~ .feature_value, data = sub_df)
      } else {
        stats::kruskal.test(PANoScore ~ .feature_value, data = sub_df)
      }
      p <- ggplot2::ggplot(sub_df, ggplot2::aes(.feature_value, PANoScore, fill = .feature_value)) +
        ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        ggplot2::geom_jitter(width = 0.15, size = 1.2, alpha = 0.55) +
        ggplot2::labs(x = "", y = "PANoScore", title = paste("PANoScore by", feature)) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(legend.position = "none")
      ggplot2::ggsave(file.path(outdir, "figures", paste0("PANoScore_by_", feature, ".pdf")), p, width = 5.5, height = 4)
      out_list[[feature]] <- data.frame(
        feature = feature,
        test = if (nlevels(sub_df$.feature_value) == 2L) "Wilcoxon" else "Kruskal-Wallis",
        statistic = unname(test$statistic),
        p.value = test$p.value,
        stringsAsFactors = FALSE
      )
    }
  }

  out_df <- do.call(rbind, out_list)
  if (!is.null(out_df) && nrow(out_df) > 0L) {
    readr::write_csv(out_df, file.path(outdir, "tables", "clinical_associations.csv"))
  }
}

plot_molecular_entity_user_style <- function(ann_df, outdir) {
  sub_df <- ann_df[!is.na(ann_df$molecular_entity) & !is.na(ann_df$PANoScore), , drop = FALSE]
  if (nrow(sub_df) < 10L || length(unique(sub_df$molecular_entity)) < 2L) {
    return(NULL)
  }
  sub_df$molecular_entity <- forcats::fct_infreq(factor(sub_df$molecular_entity))
  p <- ggplot2::ggplot(sub_df, ggplot2::aes(molecular_entity, PANoScore, fill = molecular_entity)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.15, size = 1.2, alpha = 0.55) +
    ggplot2::labs(x = "", y = "PANoScore", title = "PANoScore by molecular entity") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "none")
  ggplot2::ggsave(file.path(outdir, "figures", "PANoScore_by_molecular_entity.pdf"), p, width = 6.5, height = 4.5)
}

run_panorg_pipeline_user_style <- function(bundle, cohort_row) {
  outdir <- project_path("results", bundle$cohort_id)
  ensure_dir(outdir)
  ensure_dir(file.path(outdir, "figures"))
  ensure_dir(file.path(outdir, "tables"))
  ensure_dir(file.path(outdir, "rds"))

  sig_expr <- make_signature_expr(bundle$expr, bundle$panorg_map)
  panoscore <- compute_panoscore_user_style(bundle$expr, bundle$panorg_map)
  cluster_raw <- run_ccp_two_cluster(sig_expr, file.path(outdir, "consensus_cluster"))
  cluster_label <- rename_clusters(cluster_raw, panoscore)
  ann_df <- build_annotation_df(bundle, panoscore, cluster_label)

  readr::write_csv(ann_df, file.path(outdir, "tables", "sample_annotations.csv"))
  saveRDS(cluster_raw, file.path(outdir, "rds", "consensus_cluster_k2.rds"))

  plot_signature_heatmap_user_style(sig_expr, ann_df, file.path(outdir, "figures", "PANoRG_heatmap.pdf"))
  plot_pca_user_style(sig_expr, ann_df, file.path(outdir, "figures", "PCA_cluster.pdf"))
  plot_panoscore_by_group(ann_df, "cluster", file.path(outdir, "figures", "PANoScore_by_cluster.pdf"), "PANoScore by cluster")
  plot_molecular_entity_user_style(ann_df, outdir)

  fit_and_save_km(ann_df, "OS_time", "OS_event", "cluster", file.path(outdir, "figures", "KM_OS_cluster.pdf"), "OS by cluster")
  fit_and_save_km(ann_df, "OS_time", "OS_event", "PANoScore_group", file.path(outdir, "figures", "KM_OS_PANoScore.pdf"), "OS by PANoScore")
  fit_and_save_km(ann_df, "PFI_time", "PFI_event", "cluster", file.path(outdir, "figures", "KM_PFI_cluster.pdf"), "PFI by cluster")
  fit_and_save_km(ann_df, "PFI_time", "PFI_event", "PANoScore_group", file.path(outdir, "figures", "KM_PFI_PANoScore.pdf"), "PFI by PANoScore")

  clinical_features <- split_semicolon_field(cohort_row$clinical_features)
  fit_and_save_cox(ann_df, "OS_time", "OS_event", clinical_features, "OS", outdir)
  fit_and_save_cox(ann_df, "PFI_time", "PFI_event", clinical_features, "PFI", outdir)
  run_clinical_association_user_style(ann_df, clinical_features, outdir)

  invisible(ann_df)
}
