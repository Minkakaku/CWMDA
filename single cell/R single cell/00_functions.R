get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE), winslash = "/", mustWork = FALSE)
    script_dir <- dirname(script_path)
    candidates <- c(script_dir, dirname(script_dir), getwd(), file.path(getwd(), "pipeline"))
    for (candidate in candidates) {
      candidate <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
      if (file.exists(file.path(candidate, "00_config.R"))) {
        return(candidate)
      }
    }
  }

  if (file.exists(file.path(getwd(), "pipeline", "00_config.R"))) {
    return(normalizePath(file.path(getwd(), "pipeline"), winslash = "/", mustWork = FALSE))
  }

  if (file.exists(file.path(getwd(), "R single cell", "00_config.R"))) {
    return(normalizePath(file.path(getwd(), "R single cell"), winslash = "/", mustWork = FALSE))
  }

  if (file.exists(file.path(getwd(), "00_config.R"))) {
    return(normalizePath(getwd(), winslash = "/", mustWork = FALSE))
  }

  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_dir <- get_script_dir()
source(file.path(script_dir, "00_config.R"), local = FALSE)

missing_packages <- PACKAGES[!vapply(PACKAGES, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages(invisible(lapply(PACKAGES, library, character.only = TRUE)))

set.seed(SEED)

log_msg <- function(...) message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)

run_step <- function(label, expr) {
  start <- Sys.time()
  log_msg("START ", label)
  out <- force(expr)
  log_msg("DONE ", label, " (", round(as.numeric(difftime(Sys.time(), start, units = "mins")), 2), " min)")
  out
}

format_list_value <- function(x) {
  if (is.null(x) || length(x) == 0) return("<none>")
  paste(as.character(x), collapse = ",")
}

save_qs2 <- function(obj, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  run_step(paste("Save", path), qs_save(obj, path))
  cat("结果文件已输出: ", path, "\n", sep = "")
  invisible(path)
}

read_qs2 <- function(path) {
  if (!file.exists(path)) stop("Missing QS2: ", path, call. = FALSE)
  run_step(paste("Read", path), qs_read(path))
}

save_plot <- function(p, path_no_ext, width = 8, height = 6) {
  dir.create(dirname(path_no_ext), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(path_no_ext, ".png"), p, width = width, height = height, dpi = 300, limitsize = FALSE)
  ggsave(paste0(path_no_ext, ".pdf"), p, width = width, height = height, limitsize = FALSE)
}

detect_mito_pattern <- function(features) {
  if (sum(grepl("^MT-", features)) > 0) return("^MT-")
  if (sum(grepl("^mt-", features)) > 0) return("^mt-")
  "^[Mm][Tt]-"
}

detect_hb_pattern <- function() {
  hb_genes <- c(
    "HBA[0-9]*",
    "HBB",
    "HBD",
    "HBE[0-9]*",
    "HBG[0-9]*",
    "HBM",
    "HBQ[0-9]*",
    "HBZ",
    "Hba-[[:alnum:]]+",
    "Hbb-[[:alnum:]]+",
    "Hba[0-9]*",
    "Hbb[0-9]*",
    "Hbd",
    "Hbe[0-9]*",
    "Hbg[0-9]*",
    "Hbm",
    "Hbq[0-9]*",
    "Hbz"
  )
  paste0("^(", paste(hb_genes, collapse = "|"), ")$")
}

add_qc_metrics <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = detect_mito_pattern(rownames(obj)))
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = detect_hb_pattern())
  obj$log10GenesPerUMI <- ifelse(
    obj$nCount_RNA > 1,
    log10(obj$nFeature_RNA) / log10(obj$nCount_RNA),
    NA_real_
  )
  obj
}

read_one_sample <- function(sample) {
  matrix_dir <- file.path(DATA_DIR, sample, "outs", "filtered_feature_bc_matrix")
  if (!dir.exists(matrix_dir)) {
    warning("Skip sample without matrix: ", sample)
    return(NULL)
  }

  counts <- run_step(paste("Read10X", sample), Read10X(matrix_dir))
  if (is.list(counts)) {
    counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1]]
  }

  obj <- CreateSeuratObject(
    counts = counts,
    project = sample,
    min.cells = QC_READ$create_min_cells,
    min.features = QC_READ$create_min_features
  )
  obj <- add_qc_metrics(obj)

  if (ncol(obj) < QC_READ$skip_sample_cells_below) {
    warning("Skip small sample: ", sample)
    return(NULL)
  }

  if (QC_READ$run_doublet) {
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, nfeatures = NORMALIZE$nfeatures, verbose = FALSE)
    sce <- scDblFinder(as.SingleCellExperiment(obj))
    obj$doublet_class <- sce$scDblFinder.class
    obj$doublet_score <- sce$scDblFinder.score
  }

  if (QC_READ$run_decontx) {
    sce <- as.SingleCellExperiment(obj)
    decont <- decontX(counts(sce))
    obj$Contamination <- decont$contamination
  }

  obj
}

merge_samples <- function(obj_list) {
  obj_list <- Filter(Negate(is.null), obj_list)
  if (length(obj_list) == 0) stop("No valid sample objects.", call. = FALSE)
  if (length(obj_list) == 1) return(obj_list[[1]])
  merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = names(obj_list))
}

make_qc_plots <- function(obj) {
  plots <- list(
    violin = VlnPlot(
      obj,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb", "log10GenesPerUMI"),
      ncol = 5,
      pt.size = 0.1
    ),
    scatter_gene_umi = FeatureScatter(obj, "nCount_RNA", "nFeature_RNA"),
    scatter_mt = FeatureScatter(obj, "nCount_RNA", "percent.mt"),
    cells_per_sample = ggplot(obj@meta.data, aes(orig.ident)) +
      geom_bar() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )

  if ("doublet_class" %in% colnames(obj@meta.data)) {
    plots$doublet_ratio <- ggplot(obj@meta.data, aes(orig.ident, fill = doublet_class)) +
      geom_bar(position = "fill") +
      scale_y_continuous(labels = percent) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  if ("Contamination" %in% colnames(obj@meta.data)) {
    plots$contamination <- VlnPlot(obj, features = "Contamination", group.by = "orig.ident", pt.size = 0) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  plots
}

save_qc_plots <- function(obj, stage) {
  plots <- make_qc_plots(obj)
  for (nm in names(plots)) {
    save_plot(
      plots[[nm]],
      file.path(QC_DIR, paste0("QC_", nm, "_", stage)),
      width = ifelse(nm == "violin", 16, 10),
      height = 5
    )
  }
}

apply_qc_filter <- function(sc) {
  keep <- sc$nFeature_RNA > QC_FILTER$nFeature_RNA_min &
    sc$nFeature_RNA < QC_FILTER$nFeature_RNA_max &
    sc$percent.mt < QC_FILTER$percent_mt_max &
    sc$percent.hb < QC_FILTER$percent_hb_max &
    sc$nCount_RNA < QC_FILTER$nCount_RNA_max &
    sc$log10GenesPerUMI > QC_FILTER$log10GenesPerUMI_min

  if (QC_FILTER$keep_singlet_only && "doublet_class" %in% colnames(sc@meta.data)) {
    keep <- keep & sc$doublet_class == "singlet"
  }

  if (!is.null(QC_FILTER$contamination_max) && "Contamination" %in% colnames(sc@meta.data)) {
    keep <- keep & !is.na(sc$Contamination) & sc$Contamination < QC_FILTER$contamination_max
  }

  sc[, keep]
}

assign_group_metadata <- function(sc) {
  sc$samples <- as.character(sc$orig.ident)
  sc$group <- unname(GROUP_BY_SAMPLE[sc$samples])
  if (anyNA(sc$group)) warning("Some samples are not in GROUP_BY_SAMPLE.")
  if (exists("SAMPLE_LEVELS")) {
    sc$samples <- factor(sc$samples, levels = SAMPLE_LEVELS)
  }
  if (exists("GROUP_LEVELS")) {
    sc$group <- factor(sc$group, levels = GROUP_LEVELS)
  }
  sc
}

run_normalize_scale <- function(sc) {
  if (identical(NORMALIZE$method, "SCT")) {
    if (!NORMALIZE$source_assay %in% Assays(sc)) {
      stop("SCT requires source assay: ", NORMALIZE$source_assay, call. = FALSE)
    }

    DefaultAssay(sc) <- NORMALIZE$source_assay
    regress_vars <- intersect(NORMALIZE$regress_vars, colnames(sc@meta.data))
    has_regress <- length(regress_vars) > 0
    has_sct_method <- nzchar(NORMALIZE$sct_method)

    if (has_sct_method) {
      if (identical(NORMALIZE$sct_method, "glmGamPoi") &&
          !requireNamespace("glmGamPoi", quietly = TRUE)) {
        stop("PIPELINE_SCT_METHOD=glmGamPoi requires the glmGamPoi package.", call. = FALSE)
      }
    }

    sct_label <- paste0(
      "SCTransform source_assay=",
      NORMALIZE$source_assay,
      " assay=",
      NORMALIZE$assay,
      " nfeatures=",
      NORMALIZE$nfeatures,
      " regress=",
      format_list_value(regress_vars),
      " method=",
      ifelse(has_sct_method, NORMALIZE$sct_method, "<default>")
    )

    if (has_regress && has_sct_method) {
      sc <- run_step(
        sct_label,
        SCTransform(
          object = sc,
          assay = NORMALIZE$source_assay,
          new.assay.name = NORMALIZE$assay,
          variable.features.n = NORMALIZE$nfeatures,
          vst.flavor = NORMALIZE$sct_vst_flavor,
          return.only.var.genes = NORMALIZE$sct_return_only_var_genes,
          conserve.memory = NORMALIZE$sct_conserve_memory,
          vars.to.regress = regress_vars,
          method = NORMALIZE$sct_method,
          verbose = NORMALIZE$verbose
        )
      )
    } else if (has_regress) {
      sc <- run_step(
        sct_label,
        SCTransform(
          object = sc,
          assay = NORMALIZE$source_assay,
          new.assay.name = NORMALIZE$assay,
          variable.features.n = NORMALIZE$nfeatures,
          vst.flavor = NORMALIZE$sct_vst_flavor,
          return.only.var.genes = NORMALIZE$sct_return_only_var_genes,
          conserve.memory = NORMALIZE$sct_conserve_memory,
          vars.to.regress = regress_vars,
          verbose = NORMALIZE$verbose
        )
      )
    } else if (has_sct_method) {
      sc <- run_step(
        sct_label,
        SCTransform(
          object = sc,
          assay = NORMALIZE$source_assay,
          new.assay.name = NORMALIZE$assay,
          variable.features.n = NORMALIZE$nfeatures,
          vst.flavor = NORMALIZE$sct_vst_flavor,
          return.only.var.genes = NORMALIZE$sct_return_only_var_genes,
          conserve.memory = NORMALIZE$sct_conserve_memory,
          method = NORMALIZE$sct_method,
          verbose = NORMALIZE$verbose
        )
      )
    } else {
      sc <- run_step(
        sct_label,
        SCTransform(
          object = sc,
          assay = NORMALIZE$source_assay,
          new.assay.name = NORMALIZE$assay,
          variable.features.n = NORMALIZE$nfeatures,
          vst.flavor = NORMALIZE$sct_vst_flavor,
          return.only.var.genes = NORMALIZE$sct_return_only_var_genes,
          conserve.memory = NORMALIZE$sct_conserve_memory,
          verbose = NORMALIZE$verbose
        )
      )
    }
    DefaultAssay(sc) <- NORMALIZE$assay
    return(sc)
  }

  if (!identical(NORMALIZE$method, "NormalizeData")) {
    stop("NORMALIZE$method must be NormalizeData or SCT.", call. = FALSE)
  }

  DefaultAssay(sc) <- NORMALIZE$assay
  sc <- run_step("NormalizeData", NormalizeData(sc, verbose = NORMALIZE$verbose))
  sc <- run_step(
    paste0("FindVariableFeatures nfeatures=", NORMALIZE$nfeatures),
    FindVariableFeatures(sc, nfeatures = NORMALIZE$nfeatures, verbose = NORMALIZE$verbose)
  )
  regress_vars <- intersect(NORMALIZE$regress_vars, colnames(sc@meta.data))
  scale_mode <- tolower(NORMALIZE$scale_features)
  if (scale_mode %in% c("variable", "hvg", "variable_features")) {
    scale_features <- VariableFeatures(sc)
  } else if (scale_mode %in% c("all", "all_genes", "all_features")) {
    scale_features <- rownames(sc[[NORMALIZE$assay]])
  } else {
    stop("NORMALIZE$scale_features must be 'variable' or 'all'.", call. = FALSE)
  }

  if (length(scale_features) == 0) {
    stop("No features found before ScaleData.", call. = FALSE)
  }

  if (length(regress_vars) > 0) {
    sc <- run_step(
      paste0(
        "ScaleData method=",
        NORMALIZE$method,
        " assay=",
        NORMALIZE$assay,
        " scale_features=",
        scale_mode,
        " nfeatures=",
        length(scale_features),
        " regress=",
        format_list_value(regress_vars)
      ),
      ScaleData(
        object = sc,
        features = scale_features,
        vars.to.regress = regress_vars,
        verbose = NORMALIZE$verbose
      )
    )
  } else {
    sc <- run_step(
      paste0(
        "ScaleData method=",
      NORMALIZE$method,
      " assay=",
      NORMALIZE$assay,
      " scale_features=",
      scale_mode,
      " nfeatures=",
      length(scale_features),
        " regress=<none>"
      ),
      ScaleData(
        object = sc,
        features = scale_features,
        verbose = NORMALIZE$verbose
      )
    )
  }
  sc
}

valid_cluster_dims <- function(obj, requested_npcs, requested_dims) {
  variable_feature_count <- length(VariableFeatures(obj))
  max_npcs <- max(1L, min(requested_npcs, variable_feature_count, ncol(obj) - 1L))
  dims <- requested_dims[requested_dims <= max_npcs]
  if (length(dims) == 0) {
    stop("No valid clustering dimensions. Check NPCS/NDIMS in 00_config.R.", call. = FALSE)
  }
  list(npcs = max_npcs, dims = dims)
}

ensure_batch_column <- function(obj, batch_col, context_label) {
  if (!batch_col %in% colnames(obj@meta.data)) {
    stop(context_label, " batch column not found: ", batch_col, call. = FALSE)
  }
  batches <- unique(as.character(obj[[batch_col, drop = TRUE]]))
  batches <- batches[!is.na(batches) & nzchar(batches)]
  if (length(batches) < 2) {
    warning(context_label, " has fewer than two batches in ", batch_col, "; using PCA without correction.")
    return(FALSE)
  }
  TRUE
}

run_harmony_correction <- function(obj, settings, assay, context_label) {
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop(context_label, " batch_method='harmony' requires the harmony package.", call. = FALSE)
  }
  if (!ensure_batch_column(obj, settings$batch_group, context_label)) {
    return(list(object = obj, reduction = "pca"))
  }

  obj <- run_step(
    paste0(context_label, " RunHarmony group=", settings$batch_group),
    harmony::RunHarmony(
      obj,
      group.by.vars = settings$batch_group,
      reduction.use = "pca",
      assay.use = assay,
      verbose = FALSE
    )
  )
  list(object = obj, reduction = "harmony")
}

run_cca_correction <- function(obj, settings, normalize_settings, context_label) {
  if (!exists("IntegrateLayers", mode = "function") || !exists("CCAIntegration", mode = "function")) {
    stop(context_label, " batch_method='cca' requires Seurat v5 IntegrateLayers/CCAIntegration.", call. = FALSE)
  }
  if (!ensure_batch_column(obj, settings$batch_group, context_label)) {
    return(list(object = obj, reduction = "pca"))
  }

  assay <- normalize_settings$assay
  DefaultAssay(obj) <- assay
  batch_values <- as.factor(as.character(obj[[settings$batch_group, drop = TRUE]]))
  obj[[assay]] <- split(obj[[assay]], f = batch_values)

  normalization <- if (identical(normalize_settings$method, "SCT")) "SCT" else "LogNormalize"
  obj <- run_step(
    paste0(context_label, " IntegrateLayers CCA group=", settings$batch_group),
    IntegrateLayers(
      object = obj,
      method = CCAIntegration,
      orig.reduction = "pca",
      new.reduction = settings$cca_reduction,
      normalization.method = normalization,
      verbose = FALSE
    )
  )
  obj <- run_step(paste0(context_label, " JoinLayers after CCA"), JoinLayers(obj, assay = assay))
  list(object = obj, reduction = settings$cca_reduction)
}

run_batch_correction <- function(obj, settings, normalize_settings, context_label) {
  method <- settings$batch_method
  if (identical(method, "none")) {
    return(list(object = obj, reduction = "pca"))
  }

  if (identical(method, "harmony")) {
    return(run_harmony_correction(obj, settings, normalize_settings$assay, context_label))
  }

  if (identical(method, "cca")) {
    return(run_cca_correction(obj, settings, normalize_settings, context_label))
  }

  stop("Unknown batch correction method: ", method, call. = FALSE)
}

store_clustering_metadata <- function(obj, context, reduction, dims) {
  if (is.null(obj@misc$pipeline)) {
    obj@misc$pipeline <- list()
  }
  obj@misc$pipeline[[paste0(context, "_reduction")]] <- reduction
  obj@misc$pipeline[[paste0(context, "_dims")]] <- dims
  obj
}

choose_clustering_reduction <- function(obj, context, fallback = c("integrated.cca", "harmony", "pca")) {
  stored <- obj@misc$pipeline[[paste0(context, "_reduction")]]
  if (!is.null(stored) && stored %in% names(obj@reductions)) {
    return(stored)
  }

  available <- fallback[fallback %in% names(obj@reductions)]
  if (length(available) > 0) {
    return(available[1])
  }

  stop("No usable dimensional reduction found for clustering.", call. = FALSE)
}

choose_clustering_dims <- function(obj, context, requested_dims) {
  stored <- obj@misc$pipeline[[paste0(context, "_dims")]]
  dims <- if (!is.null(stored)) stored else requested_dims
  dims[dims > 0]
}

run_graph_cluster_umap <- function(obj, reduction, dims, resolution, context_label, run_graph_cluster = TRUE, run_umap = TRUE) {
  if (run_graph_cluster) {
    obj <- run_step(
      paste0(context_label, " FindNeighbors"),
      FindNeighbors(obj, reduction = reduction, dims = dims, verbose = FALSE)
    )
    obj <- run_step(
      paste0(context_label, " FindClusters"),
      FindClusters(obj, resolution = resolution, random.seed = SEED, verbose = FALSE)
    )
  }

  if (run_umap) {
    obj <- run_step(
      paste0(context_label, " RunUMAP"),
      RunUMAP(obj, reduction = reduction, dims = dims, seed.use = SEED, verbose = FALSE)
    )
  }

  obj
}

run_global_clustering <- function(sc) {
  pca_settings <- valid_cluster_dims(sc, CLUSTER$npcs, CLUSTER$dims)

  sc <- run_step("RunPCA", RunPCA(sc, npcs = pca_settings$npcs, verbose = FALSE))
  save_plot(ElbowPlot(sc, ndims = pca_settings$npcs), file.path(GLOBAL_DIR, "PCA_elbow"))

  corrected <- run_batch_correction(sc, CLUSTER, NORMALIZE, "Global")
  sc <- corrected$object
  reduction <- corrected$reduction
  sc <- run_graph_cluster_umap(
    sc,
    reduction = reduction,
    dims = pca_settings$dims,
    resolution = CLUSTER$resolution,
    context_label = "Global"
  )
  sc <- store_clustering_metadata(sc, "global", reduction, pca_settings$dims)
  sc
}

make_tag <- function(x) {
  tag <- gsub("[^A-Za-z0-9_]+", "_", paste(x, collapse = "_vs_"))
  gsub("^_+|_+$", "", tag)
}

plot_celltype_distribution <- function(
    sc,
    sample_col,
    celltype_col,
    output_file,
    output_dir = PROPORTION_DIR,
    width = 12,
    height = 8) {
  meta <- sc@meta.data[, c(sample_col, celltype_col), drop = FALSE]
  colnames(meta) <- c("Sample", "CellType")
  sample_levels <- if (is.factor(meta$Sample)) levels(droplevels(meta$Sample)) else unique(meta$Sample)
  meta$Sample <- factor(meta$Sample, levels = sample_levels)
  meta$CellType <- factor(meta$CellType)

  palette <- colorRampPalette(brewer.pal(12, "Paired"))(nlevels(meta$CellType))
  names(palette) <- levels(meta$CellType)
  df <- as.data.frame(table(meta$Sample, meta$CellType))
  colnames(df) <- c("Sample", "CellType", "Number")

  p <- ggplot(df, aes(Sample, Number, fill = CellType)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = palette) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(output_dir, output_file), p, width = width, height = height)
}

plot_celltype_sankey <- function(
    sc,
    source_col,
    target_col,
    path_no_ext,
    min_count = 1L,
    width = 10,
    height = 7) {
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    warning("Skip Sankey plot because the optional ggalluvial package is not installed.")
    return(invisible(NULL))
  }

  if (!all(c(source_col, target_col) %in% colnames(sc@meta.data))) {
    warning("Skip Sankey plot because columns are missing: ", source_col, ", ", target_col)
    return(invisible(NULL))
  }

  meta <- sc@meta.data[, c(source_col, target_col), drop = FALSE]
  colnames(meta) <- c("Source", "Target")
  meta <- meta[!is.na(meta$Source) & !is.na(meta$Target), , drop = FALSE]
  df <- as.data.frame(table(meta$Source, meta$Target), stringsAsFactors = FALSE)
  colnames(df) <- c("Source", "Target", "Number")
  df <- df[df$Number >= min_count, , drop = FALSE]
  if (nrow(df) == 0) {
    warning("Skip Sankey plot because no flow reaches min_count=", min_count)
    return(invisible(NULL))
  }

  target_levels <- sort(unique(as.character(df$Target)))
  palette <- colorRampPalette(brewer.pal(12, "Paired"))(length(target_levels))
  names(palette) <- target_levels

  p <- ggplot(
    df,
    aes(axis1 = Source, axis2 = Target, y = Number)
  ) +
    ggalluvial::geom_alluvium(aes(fill = Target), width = 1 / 12, alpha = 0.85) +
    ggalluvial::geom_stratum(width = 1 / 12, fill = "grey95", color = "grey45") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c(source_col, target_col), expand = c(0.08, 0.08)) +
    scale_fill_manual(values = palette) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )

  save_plot(p, path_no_ext, width = width, height = height)
  invisible(path_no_ext)
}

save_pheatmap_plot <- function(p, path_no_ext, width = 8, height = 6, dpi = 300) {
  dir.create(dirname(path_no_ext), recursive = TRUE, showWarnings = FALSE)

  png(paste0(path_no_ext, ".png"), width = width * dpi, height = height * dpi, res = dpi)
  tryCatch({
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
  }, finally = dev.off())

  pdf(paste0(path_no_ext, ".pdf"), width = width, height = height)
  tryCatch({
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
  }, finally = dev.off())
}

prepare_subcluster_input <- function(sc) {
  if (!SUBCLUSTER$sample_col %in% colnames(sc@meta.data)) {
    sc$samples <- as.character(sc$orig.ident)
  } else {
    sc$samples <- as.character(sc[[SUBCLUSTER$sample_col, drop = TRUE]])
  }

  if (SUBCLUSTER$group_col %in% colnames(sc@meta.data)) {
    sc$group <- as.character(sc[[SUBCLUSTER$group_col, drop = TRUE]])
  } else if ("group" %in% colnames(sc@meta.data)) {
    sc$group <- as.character(sc$group)
  } else {
    sc <- assign_group_metadata(sc)
  }

  sc
}

subset_subcluster_target <- function(sc) {
  if (!SUBCLUSTER$parent_col %in% colnames(sc@meta.data)) {
    stop("Missing parent annotation column: ", SUBCLUSTER$parent_col, call. = FALSE)
  }

  values <- as.character(sc[[SUBCLUSTER$parent_col, drop = TRUE]])
  cells <- colnames(sc)[!is.na(values) & values %in% SUBCLUSTER$target_celltype]
  if (length(cells) == 0) {
    stop("No cells found for SUBCLUSTER_TARGET=", paste(SUBCLUSTER$target_celltype, collapse = ","), call. = FALSE)
  }

  subset(sc, cells = cells)
}

run_subcluster_normalize_scale <- function(obj) {
  if (identical(SUBCLUSTER$method, "SCT")) {
    if (!SUBCLUSTER$source_assay %in% Assays(obj)) {
      stop("Subcluster SCT requires source assay: ", SUBCLUSTER$source_assay, call. = FALSE)
    }

    DefaultAssay(obj) <- SUBCLUSTER$source_assay
    regress_vars <- intersect(SUBCLUSTER$regress_vars, colnames(obj@meta.data))
    has_regress <- length(regress_vars) > 0
    has_sct_method <- nzchar(SUBCLUSTER$sct_method)

    if (has_sct_method) {
      if (identical(SUBCLUSTER$sct_method, "glmGamPoi") &&
          !requireNamespace("glmGamPoi", quietly = TRUE)) {
        stop("SUBCLUSTER_SCT_METHOD=glmGamPoi requires the glmGamPoi package.", call. = FALSE)
      }
    }

    sct_label <- paste0(
      "Subcluster SCTransform source_assay=",
      SUBCLUSTER$source_assay,
      " assay=",
      SUBCLUSTER$assay,
      " nfeatures=",
      SUBCLUSTER$nfeatures,
      " regress=",
      format_list_value(regress_vars),
      " method=",
      ifelse(has_sct_method, SUBCLUSTER$sct_method, "<default>")
    )

    if (has_regress && has_sct_method) {
      obj <- run_step(
        sct_label,
        SCTransform(
          object = obj,
          assay = SUBCLUSTER$source_assay,
          new.assay.name = SUBCLUSTER$assay,
          variable.features.n = SUBCLUSTER$nfeatures,
          vst.flavor = SUBCLUSTER$sct_vst_flavor,
          return.only.var.genes = SUBCLUSTER$sct_return_only_var_genes,
          conserve.memory = SUBCLUSTER$sct_conserve_memory,
          vars.to.regress = regress_vars,
          method = SUBCLUSTER$sct_method,
          verbose = TRUE
        )
      )
    } else if (has_regress) {
      obj <- run_step(
        sct_label,
        SCTransform(
          object = obj,
          assay = SUBCLUSTER$source_assay,
          new.assay.name = SUBCLUSTER$assay,
          variable.features.n = SUBCLUSTER$nfeatures,
          vst.flavor = SUBCLUSTER$sct_vst_flavor,
          return.only.var.genes = SUBCLUSTER$sct_return_only_var_genes,
          conserve.memory = SUBCLUSTER$sct_conserve_memory,
          vars.to.regress = regress_vars,
          verbose = TRUE
        )
      )
    } else if (has_sct_method) {
      obj <- run_step(
        sct_label,
        SCTransform(
          object = obj,
          assay = SUBCLUSTER$source_assay,
          new.assay.name = SUBCLUSTER$assay,
          variable.features.n = SUBCLUSTER$nfeatures,
          vst.flavor = SUBCLUSTER$sct_vst_flavor,
          return.only.var.genes = SUBCLUSTER$sct_return_only_var_genes,
          conserve.memory = SUBCLUSTER$sct_conserve_memory,
          method = SUBCLUSTER$sct_method,
          verbose = TRUE
        )
      )
    } else {
      obj <- run_step(
        sct_label,
        SCTransform(
          object = obj,
          assay = SUBCLUSTER$source_assay,
          new.assay.name = SUBCLUSTER$assay,
          variable.features.n = SUBCLUSTER$nfeatures,
          vst.flavor = SUBCLUSTER$sct_vst_flavor,
          return.only.var.genes = SUBCLUSTER$sct_return_only_var_genes,
          conserve.memory = SUBCLUSTER$sct_conserve_memory,
          verbose = TRUE
        )
      )
    }
    DefaultAssay(obj) <- SUBCLUSTER$assay
    return(obj)
  }

  if (!identical(SUBCLUSTER$method, "NormalizeData")) {
    stop("SUBCLUSTER$method must be NormalizeData or SCT.", call. = FALSE)
  }

  DefaultAssay(obj) <- SUBCLUSTER$assay
  obj <- run_step("Subcluster NormalizeData", NormalizeData(obj, verbose = TRUE))
  obj <- run_step(
    paste0("Subcluster FindVariableFeatures nfeatures=", SUBCLUSTER$nfeatures),
    FindVariableFeatures(obj, nfeatures = SUBCLUSTER$nfeatures, verbose = TRUE)
  )

  regress_vars <- intersect(SUBCLUSTER$regress_vars, colnames(obj@meta.data))
  scale_mode <- tolower(SUBCLUSTER$scale_features)
  if (scale_mode %in% c("variable", "hvg", "variable_features")) {
    scale_features <- VariableFeatures(obj)
  } else if (scale_mode %in% c("all", "all_genes", "all_features")) {
    scale_features <- rownames(obj[[SUBCLUSTER$assay]])
  } else {
    stop("SUBCLUSTER$scale_features must be 'variable' or 'all'.", call. = FALSE)
  }

  if (length(scale_features) == 0) {
    stop("No features found before subcluster ScaleData.", call. = FALSE)
  }

  if (length(regress_vars) > 0) {
    obj <- run_step(
      paste0(
        "Subcluster ScaleData method=",
        SUBCLUSTER$method,
        " assay=",
        SUBCLUSTER$assay,
        " scale_features=",
        scale_mode,
        " nfeatures=",
        length(scale_features),
        " regress=",
        format_list_value(regress_vars)
      ),
      ScaleData(
        object = obj,
        features = scale_features,
        vars.to.regress = regress_vars,
        verbose = TRUE
      )
    )
  } else {
    obj <- run_step(
      paste0(
        "Subcluster ScaleData method=",
        SUBCLUSTER$method,
        " assay=",
        SUBCLUSTER$assay,
        " scale_features=",
        scale_mode,
        " nfeatures=",
        length(scale_features),
        " regress=<none>"
      ),
      ScaleData(
        object = obj,
        features = scale_features,
        verbose = TRUE
      )
    )
  }
  obj
}

run_subcluster_clustering <- function(obj) {
  pca_settings <- valid_cluster_dims(obj, SUBCLUSTER$npcs, SUBCLUSTER$dims)

  obj <- run_step("Subcluster RunPCA", RunPCA(obj, npcs = pca_settings$npcs, verbose = FALSE))
  save_plot(ElbowPlot(obj, ndims = pca_settings$npcs), file.path(SUBCLUSTER_DIR, "PCA_elbow"))

  save_plot(
    DimPlot(
      run_step(
        "Subcluster raw RunUMAP",
        RunUMAP(obj, reduction = "pca", dims = pca_settings$dims, seed.use = SEED, verbose = FALSE)
      ),
      group.by = "samples",
      pt.size = 0.2
    ),
    file.path(SUBCLUSTER_DIR, "UMAP_raw_by_sample"),
    9,
    7
  )

  corrected <- run_batch_correction(obj, SUBCLUSTER, SUBCLUSTER, "Subcluster")
  obj <- corrected$object
  reduction <- corrected$reduction
  obj <- run_graph_cluster_umap(
    obj,
    reduction = reduction,
    dims = pca_settings$dims,
    resolution = SUBCLUSTER$resolution,
    context_label = "Subcluster"
  )
  obj <- store_clustering_metadata(obj, "subcluster", reduction, pca_settings$dims)

  obj$subcluster_parent <- paste(SUBCLUSTER$target_celltype, collapse = ",")
  obj$cluster_id <- as.character(Idents(obj))
  obj
}

prepare_marker_assay <- function(obj, assay) {
  DefaultAssay(obj) <- assay
  assay_obj <- obj[[assay]]

  if (inherits(assay_obj, "SCTAssay")) {
    obj <- tryCatch(
      run_step(
        paste0("PrepSCTFindMarkers assay=", assay),
        PrepSCTFindMarkers(obj, assay = assay, verbose = TRUE)
      ),
      error = function(e) {
        warning("PrepSCTFindMarkers skipped: ", conditionMessage(e))
        obj
      }
    )
    assay_obj <- obj[[assay]]
  }

  obj <- JoinLayers(obj, assay = assay)
  obj
}

sort_cluster_ids <- function(x) {
  x <- unique(as.character(x))
  numeric_x <- suppressWarnings(as.numeric(x))
  if (all(!is.na(numeric_x))) {
    return(x[order(numeric_x)])
  }
  sort(x)
}

sort_markers_by_cluster <- function(markers) {
  if (is.null(markers) || nrow(markers) == 0) {
    return(markers)
  }

  logfc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(markers))[1]
  if (!is.na(logfc_col) && "cluster" %in% colnames(markers)) {
    markers <- markers |>
      group_by(cluster) |>
      arrange(desc(.data[[logfc_col]]), .by_group = TRUE)
  }
  markers
}

top_marker_summary <- function(markers, top_n = 10L) {
  if (is.null(markers) || nrow(markers) == 0 || !"cluster" %in% colnames(markers) || !"gene" %in% colnames(markers)) {
    return(character(0))
  }

  markers <- sort_markers_by_cluster(markers)
  genes_by_cluster <- split(as.character(markers$gene), as.character(markers$cluster))
  vapply(
    genes_by_cluster,
    function(genes) paste(head(unique(genes), top_n), collapse = ";"),
    character(1)
  )
}

write_annotation_template <- function(
    cluster_ids,
    label_values,
    path,
    cluster_col = "cluster",
    label_col = "cell_type",
    markers = NULL,
    top_n = 10L) {
  cluster_ids <- sort_cluster_ids(cluster_ids)
  label_values <- unname(label_values[cluster_ids])
  label_values[is.na(label_values)] <- ""

  marker_summary <- top_marker_summary(markers, top_n = top_n)
  top_markers <- unname(marker_summary[cluster_ids])
  top_markers[is.na(top_markers)] <- ""

  template <- data.frame(
    cluster = cluster_ids,
    cell_type = label_values,
    top_markers = top_markers,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  colnames(template)[1:2] <- c(cluster_col, label_col)

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(template, path, row.names = FALSE, quote = FALSE)
  path
}

write_cluster_annotation_template <- function(sc, markers = NULL) {
  cluster_ids <- sort_cluster_ids(sc$seurat_clusters)
  write_annotation_template(
    cluster_ids = cluster_ids,
    label_values = CLUSTER_TO_CELLTYPE,
    path = ANNOTATION$template_path,
    cluster_col = "cluster",
    label_col = "cell_type",
    markers = markers,
    top_n = MARKER_PARAMS$top_n_for_template
  )
}

write_subcluster_annotation_template <- function(obj, markers = NULL) {
  cluster_ids <- sort_cluster_ids(obj$cluster_id)
  write.csv(
    data.frame(cluster_id = cluster_ids, sub_cell_type = "", stringsAsFactors = FALSE),
    SUBCLUSTER$template_path,
    row.names = FALSE,
    quote = FALSE
  )
  if (!is.null(markers) && nrow(markers) > 0) {
    write_annotation_template(
      cluster_ids = cluster_ids,
      label_values = SUBCLUSTER$annotation,
      path = SUBCLUSTER$template_path,
      cluster_col = "cluster_id",
      label_col = "sub_cell_type",
      markers = markers,
      top_n = SUBCLUSTER$marker_params$top_n_for_template
    )
  }
  SUBCLUSTER$template_path
}

read_template_mapping <- function(path, cluster_col, label_col) {
  if (!file.exists(path)) {
    return(character(0))
  }

  template <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c(cluster_col, label_col)
  missing_cols <- setdiff(required_cols, colnames(template))
  if (length(missing_cols) > 0) {
    stop("Annotation template missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  labels <- trimws(as.character(template[[label_col]]))
  ids <- as.character(template[[cluster_col]])
  keep <- !is.na(ids) & nzchar(ids) & !is.na(labels) & nzchar(labels)
  mapping <- stats::setNames(labels[keep], ids[keep])
  mapping[!duplicated(names(mapping), fromLast = TRUE)]
}

resolve_annotation_mapping <- function(config_mapping, template_path, cluster_col, label_col) {
  template_mapping <- read_template_mapping(template_path, cluster_col, label_col)
  if (length(template_mapping) > 0) {
    return(template_mapping)
  }

  config_mapping <- config_mapping[!is.na(config_mapping) & nzchar(config_mapping)]
  config_mapping
}

apply_annotation <- function(
    obj,
    cluster_values,
    config_mapping,
    template_path,
    output_template_fn,
    cluster_col,
    label_col) {
  clusters <- sort_cluster_ids(cluster_values)
  mapping <- resolve_annotation_mapping(config_mapping, template_path, cluster_col, label_col)
  missing_clusters <- setdiff(clusters, names(mapping))
  if (length(missing_clusters) > 0) {
    template <- output_template_fn()
    stop(
      label_col,
      " missing clusters: ",
      paste(missing_clusters, collapse = ", "),
      ". Fill template or 00_config.R, then rerun. Template: ",
      template,
      call. = FALSE
    )
  }

  labels <- unname(mapping[as.character(cluster_values)])
  obj[[label_col]] <- labels
  obj
}

cleanup_removed_clusters <- function(obj, cleanup, context, output_dir) {
  remove_clusters <- cleanup$remove_clusters
  if (is.null(remove_clusters) || length(remove_clusters) == 0) {
    return(obj)
  }

  Idents(obj) <- "seurat_clusters"
  cluster_values <- as.character(Idents(obj))
  present_clusters <- intersect(remove_clusters, unique(cluster_values))
  if (length(present_clusters) == 0) {
    warning("No configured clusters found for removal: ", paste(remove_clusters, collapse = ", "))
    return(obj)
  }

  backup_col <- cleanup$cluster_backup_col
  if (!backup_col %in% colnames(obj@meta.data)) {
    obj[[backup_col]] <- cluster_values
  }

  removed_cells <- colnames(obj)[cluster_values %in% present_clusters]
  removed_summary <- as.data.frame(table(cluster_values[cluster_values %in% present_clusters]))
  colnames(removed_summary) <- c("removed_cluster", "cell_count")
  write.csv(
    removed_summary,
    file.path(output_dir, paste0(context, "_removed_clusters.csv")),
    row.names = FALSE,
    quote = FALSE
  )

  keep_cells <- setdiff(colnames(obj), removed_cells)
  obj <- subset(obj, cells = keep_cells)
  obj$seurat_clusters <- droplevels(factor(as.character(obj$seurat_clusters)))
  Idents(obj) <- "seurat_clusters"
  log_msg(
    "Removed clusters from ",
    context,
    ": ",
    paste(present_clusters, collapse = ", "),
    "; removed cells=",
    length(removed_cells)
  )

  rerun_graph_cluster <- isTRUE(cleanup$rerun_neighbors_clusters)
  rerun_umap <- isTRUE(cleanup$rerun_umap)
  if (!rerun_graph_cluster && !rerun_umap) {
    return(obj)
  }

  reduction <- choose_clustering_reduction(obj, context)
  dims <- choose_clustering_dims(obj, context, cleanup$rerun_dims)
  reduction_dims <- ncol(Embeddings(obj, reduction = reduction))
  dims <- dims[dims <= reduction_dims]
  if (length(dims) == 0) {
    stop("No valid dims left for cleanup rerun using reduction: ", reduction, call. = FALSE)
  }

  obj <- run_graph_cluster_umap(
    obj,
    reduction = reduction,
    dims = dims,
    resolution = cleanup$rerun_resolution,
    context_label = paste0(context, " cleanup"),
    run_graph_cluster = rerun_graph_cluster,
    run_umap = rerun_umap
  )
  store_clustering_metadata(obj, context, reduction, dims)
}

find_all_markers_wrapped <- function(obj, assay, params, label) {
  obj <- prepare_marker_assay(obj, assay)
  markers <- run_step(
    label,
    FindAllMarkers(
      obj,
      assay = assay,
      only.pos = params$only_pos,
      test.use = params$test_use,
      min.pct = params$min_pct,
      logfc.threshold = params$logfc_threshold,
      return.thresh = params$return_thresh
    )
  )
  sort_markers_by_cluster(markers)
}

save_broad_marker_dotplot <- function(sc) {
  broad_markers <- lapply(BROAD_MARKER_GENES, function(x) intersect(x, rownames(sc)))
  broad_markers <- broad_markers[lengths(broad_markers) > 0]

  if (length(broad_markers) > 0) {
    save_plot(
      DotPlot(sc, features = broad_markers, assay = NORMALIZE$assay) + RotatedAxis(),
      file.path(MARKER_DIR, "cell_markers"),
      18,
      6
    )
  }
}

save_annotation_plots <- function(sc) {
  save_plot(
    DimPlot(sc, group.by = "cell_type", label = TRUE, repel = TRUE),
    file.path(ANNOTATION_DIR, "UMAP_by_celltype"),
    14,
    10
  )
  save_plot(
    DimPlot(sc, group.by = "cell_type", split.by = "samples", ncol = 4),
    file.path(ANNOTATION_DIR, "UMAP_split_by_sample"),
    18,
    14
  )
  save_plot(
    DimPlot(sc, group.by = "cell_type", split.by = "group", ncol = 3),
    file.path(ANNOTATION_DIR, "UMAP_split_by_group"),
    16,
    8
  )
}

save_subcluster_annotation_plots <- function(obj) {
  save_plot(
    DimPlot(obj, group.by = "sub_cell_type", label = TRUE, repel = TRUE),
    file.path(SUBCLUSTER_DIR, "UMAP_by_sub_cell_type"),
    9,
    7
  )
  save_plot(
    DimPlot(obj, group.by = "sub_cell_type", split.by = "group", ncol = 3),
    file.path(SUBCLUSTER_DIR, "UMAP_split_by_group"),
    15,
    7
  )
  save_plot(
    DimPlot(obj, group.by = "sub_cell_type", split.by = "samples", ncol = 3),
    file.path(SUBCLUSTER_DIR, "UMAP_split_by_samples"),
    12,
    8
  )
}

save_proportion_tables <- function(obj, sample_col, celltype_col, output_dir, prefix) {
  write.csv(
    as.data.frame.matrix(table(obj[[sample_col, drop = TRUE]], obj[[celltype_col, drop = TRUE]])),
    file.path(output_dir, paste0(prefix, "_count.csv")),
    quote = FALSE
  )
}

save_proportion_outputs <- function(
    obj,
    comparisons,
    sample_col,
    group_col,
    celltype_col,
    output_dir,
    prefix,
    make_sankey = TRUE,
    sankey_min_count = 1L) {
  plot_celltype_distribution(
    obj,
    sample_col = sample_col,
    celltype_col = celltype_col,
    output_file = paste0(prefix, "_by_sample.pdf"),
    output_dir = output_dir
  )
  plot_celltype_distribution(
    obj,
    sample_col = group_col,
    celltype_col = celltype_col,
    output_file = paste0(prefix, "_by_group.pdf"),
    output_dir = output_dir
  )
  save_proportion_tables(obj, sample_col, celltype_col, output_dir, paste0(prefix, "_by_sample"))

  if (make_sankey) {
    plot_celltype_sankey(
      obj,
      source_col = group_col,
      target_col = celltype_col,
      path_no_ext = file.path(output_dir, paste0(prefix, "_sankey_group_to_celltype")),
      min_count = sankey_min_count,
      width = 12,
      height = 8
    )
  }

  for (groups in comparisons) {
    group_values <- as.character(obj[[group_col, drop = TRUE]])
    missing_groups <- setdiff(groups, unique(group_values))
    if (length(missing_groups) > 0) {
      warning("Skip comparison with missing groups: ", paste(groups, collapse = " vs "))
      next
    }

    obj_sub <- subset(obj, cells = colnames(obj)[group_values %in% groups])
    obj_sub[[group_col]] <- factor(obj_sub[[group_col, drop = TRUE]], levels = groups)
    tag <- make_tag(groups)

    plot_celltype_distribution(
      obj_sub,
      sample_col = sample_col,
      celltype_col = celltype_col,
      output_file = paste0(prefix, "_by_sample_", tag, ".pdf"),
      output_dir = output_dir
    )
    plot_celltype_distribution(
      obj_sub,
      sample_col = group_col,
      celltype_col = celltype_col,
      output_file = paste0(prefix, "_by_group_", tag, ".pdf"),
      output_dir = output_dir
    )
    save_proportion_tables(
      obj_sub,
      sample_col,
      celltype_col,
      output_dir,
      paste0(prefix, "_by_sample_", tag)
    )

    if (make_sankey) {
      plot_celltype_sankey(
        obj_sub,
        source_col = group_col,
        target_col = celltype_col,
        path_no_ext = file.path(output_dir, paste0(prefix, "_sankey_group_to_celltype_", tag)),
        min_count = sankey_min_count,
        width = 12,
        height = 8
      )
    }
  }
}

step_qc1 <- function() {
  log_msg("Step 01 QC before filtering started")

  if (!dir.exists(DATA_DIR)) stop("Data dir not found: ", DATA_DIR, call. = FALSE)
  samples <- list.dirs(DATA_DIR, recursive = FALSE, full.names = FALSE)
  if (length(samples) == 0) stop("No sample dirs under ", DATA_DIR, call. = FALSE)

  obj_list <- setNames(lapply(samples, read_one_sample), samples)
  sc <- run_step("Merge samples", merge_samples(obj_list))
  sc <- assign_group_metadata(sc)

  save_qc_plots(sc, "before")
  save_qs2(sc, QS2_PRE_QC)

  log_msg("Step 01 QC before filtering done. Edit QC_FILTER in 00_config.R, then run main/01_qc2.R.")
}

step_qc2 <- function() {
  log_msg("Step 01 QC filtering started")

  sc <- read_qs2(QS2_PRE_QC)
  sc <- run_step("Apply QC_FILTER", apply_qc_filter(sc))
  sc <- assign_group_metadata(sc)

  if (ncol(sc) == 0) {
    stop("No cells left after QC_FILTER. Relax parameters in 00_config.R.", call. = FALSE)
  }

  save_qc_plots(sc, "after")
  save_qs2(sc, QS2_QC)

  log_msg("Step 01 QC filtering done. Kept cells: ", ncol(sc))
}

step_globalcluster <- function() {
  log_msg("Step 02 global clustering started")

  sc <- read_qs2(QS2_QC)
  sc <- assign_group_metadata(sc)
  sc <- run_normalize_scale(sc)
  sc <- run_global_clustering(sc)

  save_plot(DimPlot(sc, label = TRUE, repel = TRUE), file.path(GLOBAL_DIR, "UMAP_clusters"))
  save_plot(DimPlot(sc, group.by = "samples"), file.path(GLOBAL_DIR, "UMAP_by_sample"))
  if ("Contamination" %in% colnames(sc@meta.data)) {
    save_plot(FeaturePlot(sc, features = "Contamination"), file.path(GLOBAL_DIR, "Feature_Contamination"), 5, 5)
  }

  save_qs2(sc, QS2_CLUSTER)
  log_msg("Step 02 done. Clusters: ", length(unique(sc$seurat_clusters)))
}

step_findmarkers <- function() {
  log_msg("Step 03 find markers started")

  sc <- read_qs2(QS2_CLUSTER)
  DefaultAssay(sc) <- MARKER_PARAMS$assay
  Idents(sc) <- "seurat_clusters"

  markers <- find_all_markers_wrapped(sc, MARKER_PARAMS$assay, MARKER_PARAMS, "FindAllMarkers")
  write.csv(markers, file.path(MARKER_DIR, "AllClusters_markers.csv"), row.names = FALSE)
  save_qs2(markers, QS2_MARKERS)

  template <- write_cluster_annotation_template(sc, markers)
  save_broad_marker_dotplot(sc)

  log_msg("Step 03 done. Annotation template: ", template)
}

step_cell_annotation <- function() {
  log_msg("Step 04 cell annotation started")

  sc <- read_qs2(QS2_CLUSTER)
  sc <- assign_group_metadata(sc)
  clusters <- sort_cluster_ids(sc$seurat_clusters)
  mapping_used <- resolve_annotation_mapping(
    CLUSTER_TO_CELLTYPE,
    ANNOTATION$template_path,
    "cluster",
    "cell_type"
  )

  sc <- apply_annotation(
    obj = sc,
    cluster_values = as.character(sc$seurat_clusters),
    config_mapping = CLUSTER_TO_CELLTYPE,
    template_path = ANNOTATION$template_path,
    output_template_fn = function() write_cluster_annotation_template(sc),
    cluster_col = "cluster",
    label_col = "cell_type"
  )
  sc <- cleanup_removed_clusters(sc, ANNOTATION, "global", ANNOTATION_DIR)

  mapping <- data.frame(cluster = names(mapping_used), cell_type = unname(mapping_used))
  write.csv(mapping, ANNOTATION$mapping_path, row.names = FALSE, quote = FALSE)
  write.csv(
    as.data.frame(table(cluster = sc$seurat_clusters, cell_type = sc$cell_type)),
    file.path(ANNOTATION_DIR, "cluster_celltype_count.csv"),
    row.names = FALSE,
    quote = FALSE
  )

  save_annotation_plots(sc)
  save_qs2(sc, QS2_ANNOTATED)

  log_msg("Step 04 done.")
}

step_cell_proportion <- function() {
  log_msg("Step 05 cell proportion started")

  sc <- read_qs2(QS2_ANNOTATED)
  sc <- assign_group_metadata(sc)

  if (!"cell_type" %in% colnames(sc@meta.data)) {
    stop("Missing cell_type. Run main/04_cell_annotation.R first.", call. = FALSE)
  }

  save_proportion_outputs(
    sc,
    comparisons = PROPORTION$comparisons,
    sample_col = "samples",
    group_col = "group",
    celltype_col = "cell_type",
    output_dir = PROPORTION_DIR,
    prefix = "sc_proportion",
    make_sankey = PROPORTION$make_sankey,
    sankey_min_count = PROPORTION$sankey_min_count
  )

  save_qs2(sc, QS2_FINAL)
  log_msg("Step 05 done.")
}

step_subcluster_clustering <- function() {
  log_msg("Step 00 subcluster clustering started: ", SUBCLUSTER_NAME)

  sc <- read_qs2(SUBCLUSTER$input_qs2)
  sc <- prepare_subcluster_input(sc)
  obj <- subset_subcluster_target(sc)

  obj <- run_subcluster_normalize_scale(obj)
  obj <- run_subcluster_clustering(obj)

  save_plot(DimPlot(obj, label = TRUE, repel = TRUE), file.path(SUBCLUSTER_DIR, "UMAP_clusters"), 9, 7)
  save_plot(DimPlot(obj, group.by = "samples"), file.path(SUBCLUSTER_DIR, "UMAP_by_sample"), 9, 7)
  save_plot(DimPlot(obj, group.by = "group"), file.path(SUBCLUSTER_DIR, "UMAP_by_group"), 9, 7)

  save_qs2(obj, QS2_SUBCLUSTER_CLUSTERED)
  log_msg("Step 00 done. Check ", file.path(SUBCLUSTER_DIR, "UMAP_clusters.pdf"))
}

step_subcluster_findmarkers <- function() {
  log_msg("Step 01 subcluster markers started: ", SUBCLUSTER_NAME)

  obj <- read_qs2(QS2_SUBCLUSTER_CLUSTERED)
  Idents(obj) <- "seurat_clusters"
  obj$cluster_id <- as.character(Idents(obj))

  markers <- find_all_markers_wrapped(
    obj,
    SUBCLUSTER$assay,
    SUBCLUSTER$marker_params,
    "Subcluster FindAllMarkers"
  )

  write.csv(markers, file.path(SUBCLUSTER_DIR, "markers_all_clusters.csv"), row.names = FALSE)
  save_qs2(markers, QS2_SUBCLUSTER_MARKERS)
  template <- write_subcluster_annotation_template(obj, markers)

  log_msg("Step 01 done. Subcluster annotation template: ", template)
}

step_subcluster_annotation <- function() {
  log_msg("Step 02 subcluster annotation started: ", SUBCLUSTER_NAME)

  obj <- read_qs2(QS2_SUBCLUSTER_CLUSTERED)
  obj <- prepare_subcluster_input(obj)
  Idents(obj) <- "seurat_clusters"
  obj$cluster_id <- as.character(Idents(obj))
  mapping_used <- resolve_annotation_mapping(
    SUBCLUSTER$annotation,
    SUBCLUSTER$template_path,
    "cluster_id",
    "sub_cell_type"
  )

  obj <- apply_annotation(
    obj = obj,
    cluster_values = obj$cluster_id,
    config_mapping = SUBCLUSTER$annotation,
    template_path = SUBCLUSTER$template_path,
    output_template_fn = function() write_subcluster_annotation_template(obj),
    cluster_col = "cluster_id",
    label_col = "sub_cell_type"
  )
  obj <- cleanup_removed_clusters(obj, SUBCLUSTER, "subcluster", SUBCLUSTER_DIR)

  write.csv(
    data.frame(cluster_id = names(mapping_used), sub_cell_type = unname(mapping_used)),
    file.path(SUBCLUSTER_DIR, "subcluster_annotation_mapping.csv"),
    row.names = FALSE,
    quote = FALSE
  )

  save_subcluster_annotation_plots(obj)
  save_qs2(obj, QS2_SUBCLUSTER_ANNOTATED)

  log_msg("Step 02 done.")
}

step_subcluster_proportion <- function() {
  log_msg("Step 03 subcluster proportion started: ", SUBCLUSTER_NAME)

  obj <- read_qs2(QS2_SUBCLUSTER_ANNOTATED)
  obj <- prepare_subcluster_input(obj)

  if (!"sub_cell_type" %in% colnames(obj@meta.data)) {
    stop("Missing sub_cell_type. Run subcluster/02_subcluster_annotation.R first.", call. = FALSE)
  }

  save_proportion_outputs(
    obj,
    comparisons = SUBCLUSTER$comparisons,
    sample_col = "samples",
    group_col = "group",
    celltype_col = "sub_cell_type",
    output_dir = SUBCLUSTER_DIR,
    prefix = paste0(SUBCLUSTER_NAME, "_subtype"),
    make_sankey = SUBCLUSTER$make_sankey,
    sankey_min_count = PROPORTION$sankey_min_count
  )

  log_msg("Step 03 done.")
}
