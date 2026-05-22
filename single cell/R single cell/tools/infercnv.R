rm(list = ls())

cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(cmd_file) > 0) {
  dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", cmd_file[1]), fixed = TRUE), winslash = "/", mustWork = FALSE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

pipeline_candidates <- unique(c(
  script_dir,
  file.path(script_dir, ".."),
  file.path(script_dir, "..", "main"),
  getwd(),
  file.path(getwd(), ".."),
  file.path(getwd(), "pipeline"),
  file.path(getwd(), "R single cell")
))
pipeline_dir <- ""
for (candidate in pipeline_candidates) {
  candidate <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
  if (file.exists(file.path(candidate, "00_functions.R"))) {
    pipeline_dir <- candidate
    break
  }
}
if (!nzchar(pipeline_dir)) {
  stop("Cannot find 00_functions.R from ", script_dir, call. = FALSE)
}

if (!nzchar(Sys.getenv("PROJECT_DIR", unset = ""))) {
  Sys.setenv(PROJECT_DIR = pipeline_dir)
}
source(file.path(pipeline_dir, "00_functions.R"), local = FALSE)

script_args <- commandArgs(trailingOnly = TRUE)

if (!requireNamespace("infercnv", quietly = TRUE)) {
  stop("Missing R package: infercnv", call. = FALSE)
}
if (!requireNamespace("AnnoProbe", quietly = TRUE) && !nzchar(Sys.getenv("INFERCNV_GENE_ORDER_FILE", unset = ""))) {
  stop("Missing R package: AnnoProbe. Install it or set INFERCNV_GENE_ORDER_FILE.", call. = FALSE)
}
library(infercnv)
if (requireNamespace("AnnoProbe", quietly = TRUE)) library(AnnoProbe)

infercnv_run <- get("run", envir = asNamespace("infercnv"))

is_absolute_path <- function(path) {
  grepl("^[A-Za-z]:[/\\\\]", path) || grepl("^[/\\\\]", path)
}

resolve_existing_file <- function(path, bases = c(getwd(), pipeline_dir, script_dir)) {
  if (!nzchar(path)) return(path)
  candidates <- if (is_absolute_path(path)) {
    path
  } else {
    unique(c(path, file.path(bases, path)))
  }
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0) return(path)
  normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

resolve_output_dir <- function(path) {
  if (is_absolute_path(path)) return(normalizePath(path, winslash = "/", mustWork = FALSE))
  normalizePath(file.path(getwd(), path), winslash = "/", mustWork = FALSE)
}

join_layers_if_available <- function(obj, assay) {
  if (exists("JoinLayers", mode = "function")) {
    return(JoinLayers(obj, assay = assay))
  }
  obj
}

get_counts_matrix <- function(obj, assay) {
  if (exists("LayerData", mode = "function")) {
    out <- tryCatch(LayerData(obj, assay = assay, layer = "counts"), error = function(e) NULL)
    if (!is.null(out)) return(out)
  }
  GetAssayData(obj, assay = assay, slot = "counts")
}

pick_existing_metadata_col <- function(meta, candidates, label) {
  hit <- candidates[candidates %in% colnames(meta)][1]
  if (!is.na(hit) && nzchar(hit)) return(hit)
  stop(
    "Missing ", label, " column. Tried: ", paste(candidates, collapse = ", "),
    ". Run main/04_cell_annotation.R first or set INFERCNV_ANNOTATION_COL.",
    call. = FALSE
  )
}

match_existing_levels <- function(values, requested, label) {
  requested <- unique(as.character(requested))
  existing <- unique(as.character(values))
  matched <- requested[requested %in% existing]
  if (length(matched) == length(requested)) return(matched)

  canonical <- function(x) gsub("[^a-z0-9]+", "", tolower(x))
  existing_by_key <- setNames(existing, canonical(existing))
  filled <- matched
  unresolved <- character(0)
  for (item in setdiff(requested, matched)) {
    key <- canonical(item)
    if (key %in% names(existing_by_key)) {
      filled <- c(filled, existing_by_key[[key]])
    } else {
      unresolved <- c(unresolved, item)
    }
  }
  filled <- unique(filled)
  missing <- unique(unresolved)
  if (length(missing) > 0) {
    stop(
      "Missing ", label, " in annotation column: ", paste(missing, collapse = ", "),
      ". Available labels include: ", paste(head(sort(existing), 30), collapse = ", "),
      call. = FALSE
    )
  }
  filled
}

infercnv_input_arg <- ""
qs2_args <- script_args[grepl("\\.qs2$", script_args, ignore.case = TRUE)]
if (length(qs2_args) > 0) {
  infercnv_input_arg <- qs2_args[1]
} else if (length(script_args) > 0 && file.exists(resolve_existing_file(script_args[1]))) {
  infercnv_input_arg <- script_args[1]
}

first_existing_file <- function(paths) {
  resolved <- vapply(paths, resolve_existing_file, character(1))
  hits <- resolved[file.exists(resolved)]
  if (length(hits) == 0) {
    stop("No input QS2 file found. Checked: ", paste(paths, collapse = ", "), call. = FALSE)
  }
  hits[1]
}

build_gene_order_file <- function(genes, output_file, species) {
  gene_order_override <- env("INFERCNV_GENE_ORDER_FILE", "")
  if (nzchar(gene_order_override)) {
    gene_order_override <- resolve_existing_file(gene_order_override)
    if (!file.exists(gene_order_override)) {
      stop("INFERCNV_GENE_ORDER_FILE does not exist: ", gene_order_override, call. = FALSE)
    }
    return(gene_order_override)
  }

  anno_species <- if (tolower(species) %in% c("human", "hg38", "hg19")) "human" else species
  gene_info <- annoGene(genes, "SYMBOL", anno_species)

  symbol_col <- intersect(c("SYMBOL", "symbol", "gene", "Gene"), colnames(gene_info))[1]
  chr_col <- intersect(c("chr", "chromosome", "seqnames"), colnames(gene_info))[1]
  start_col <- intersect(c("start", "txStart"), colnames(gene_info))[1]
  end_col <- intersect(c("end", "txEnd"), colnames(gene_info))[1]
  required_cols <- c(symbol_col, chr_col, start_col, end_col)
  if (anyNA(required_cols)) {
    stop("AnnoProbe output is missing SYMBOL/chr/start/end columns.", call. = FALSE)
  }

  gene_info <- gene_info[, required_cols]
  colnames(gene_info) <- c("gene", "chr", "start", "end")
  gene_info <- gene_info[!is.na(gene_info$gene) & !is.na(gene_info$chr), ]
  gene_info <- gene_info[!gene_info$chr %in% c("chrM", "MT", "chrX", "X", "chrY", "Y"), ]
  gene_info$chr <- ifelse(grepl("^chr", gene_info$chr), gene_info$chr, paste0("chr", gene_info$chr))
  gene_info$chr_num <- suppressWarnings(as.integer(sub("^chr", "", gene_info$chr)))
  gene_info <- gene_info[!is.na(gene_info$chr_num), ]
  gene_info <- gene_info[order(gene_info$chr_num, gene_info$start), c("gene", "chr", "start", "end")]
  gene_info <- gene_info[!duplicated(gene_info$gene), ]
  gene_info <- gene_info[gene_info$gene %in% genes, ]

  if (nrow(gene_info) == 0) {
    stop("No genes remained after building inferCNV gene order.", call. = FALSE)
  }

  write.table(gene_info, output_file, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  output_file
}

sample_cells_by_group <- function(meta, group_col, max_cells_per_group) {
  if (max_cells_per_group <= 0) return(rownames(meta))

  set.seed(SEED)
  unlist(lapply(split(rownames(meta), meta[[group_col]]), function(x) {
    if (length(x) <= max_cells_per_group) return(x)
    sample(x, max_cells_per_group)
  }), use.names = FALSE)
}

run_epithelial_clustering <- function(obj) {
  DefaultAssay(obj) <- infercnv_assay
  obj <- NormalizeData(obj, verbose = TRUE)
  obj <- FindVariableFeatures(obj, nfeatures = infercnv_nfeatures, verbose = TRUE)
  scale_features <- VariableFeatures(obj)
  regress_vars <- intersect(infercnv_regress_vars, colnames(obj@meta.data))

  if (length(regress_vars) > 0) {
    obj <- ScaleData(obj, features = scale_features, vars.to.regress = regress_vars, verbose = TRUE)
  } else {
    obj <- ScaleData(obj, features = scale_features, verbose = TRUE)
  }

  npcs <- min(infercnv_npcs, length(VariableFeatures(obj)), ncol(obj) - 1)
  dims <- seq_len(min(infercnv_ndims, npcs))
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = npcs, verbose = FALSE)

  reduction <- "pca"
  if (infercnv_use_harmony) {
    if (!requireNamespace("harmony", quietly = TRUE)) {
      stop("INFERCNV_USE_HARMONY=TRUE requires the harmony package.", call. = FALSE)
    }
    obj <- harmony::RunHarmony(obj, group.by.vars = infercnv_harmony_group, reduction.use = "pca")
    reduction <- "harmony"
  }

  obj <- RunUMAP(obj, reduction = reduction, dims = dims, seed.use = SEED)
  obj <- FindNeighbors(obj, reduction = reduction, dims = dims)
  obj <- FindClusters(obj, resolution = infercnv_resolution, random.seed = SEED)
  obj
}

calculate_cnv_score <- function(infercnv_obj, ref_group_names) {
  expr <- infercnv_obj@expr.data
  ref_indices <- unlist(lapply(ref_group_names, function(ref_group_name) {
    indices <- infercnv_obj@reference_grouped_cell_indices[[ref_group_name]]
    if (is.null(indices) || length(indices) == 0) {
      stop("No reference cells found in inferCNV object for: ", ref_group_name, call. = FALSE)
    }
    indices
  }), use.names = FALSE)

  ref_cells <- if (is.numeric(ref_indices)) colnames(expr)[ref_indices] else as.character(ref_indices)
  ref_cells <- intersect(ref_cells, colnames(expr))
  if (length(ref_cells) == 0) {
    stop("Reference cells are not present in inferCNV expression matrix.", call. = FALSE)
  }

  ref_mean <- rowMeans(expr[, ref_cells, drop = FALSE])
  ref_sd <- apply(expr[, ref_cells, drop = FALSE], 1, sd)
  ref_sd[is.na(ref_sd) | ref_sd == 0] <- median(ref_sd[ref_sd > 0], na.rm = TRUE)
  ref_sd[is.na(ref_sd)] <- 0

  low2 <- ref_mean - 4 * ref_sd
  low1 <- ref_mean - 2 * ref_sd
  high1 <- ref_mean + 2 * ref_sd
  high2 <- ref_mean + 4 * ref_sd

  score_mat <- matrix(0, nrow = nrow(expr), ncol = ncol(expr), dimnames = dimnames(expr))
  score_mat[expr <= low1 | expr >= high1] <- 1
  score_mat[expr <= low2 | expr >= high2] <- 2
  colSums(score_mat)
}

input_qs2 <- env("INFERCNV_INPUT_QS2", infercnv_input_arg)
if (!nzchar(input_qs2)) input_qs2 <- first_existing_file(c(QS2_FINAL, QS2_ANNOTATED, QS2_CLUSTER))
input_qs2 <- resolve_existing_file(input_qs2)

output_dir <- resolve_output_dir(env("INFERCNV_OUTPUT_DIR", "11_InferCNV_T_B_Mast_REF"))
infercnv_output_dir <- file.path(output_dir, "03.infercnv_output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

infercnv_assay <- env("INFERCNV_ASSAY", "RNA")
infercnv_annotation_col <- env("INFERCNV_ANNOTATION_COL", "")
infercnv_epithelial_label <- env("INFERCNV_EPITHELIAL_LABEL", "Epithelial")
infercnv_ref_celltypes <- env_csv("INFERCNV_REF_CELLTYPES", c("T_cell", "B_cell", "Mast"))
infercnv_ref_label <- env("INFERCNV_REF_LABEL", "Ref.T_B_Mast")
infercnv_observation_label <- env("INFERCNV_OBSERVATION_LABEL", "Epithelial")
species <- env("INFERCNV_SPECIES", "mouse")
max_ref_cells_per_type <- env_int("INFERCNV_MAX_REF_CELLS_PER_TYPE", 2000L)
max_epithelial_cells <- env_int("INFERCNV_MAX_EPITHELIAL_CELLS", 0L)

infercnv_nfeatures <- env_int("INFERCNV_NFEATURES", 2000L)
infercnv_regress_vars <- env_csv("INFERCNV_REGRESS_VARS", character(0))
infercnv_npcs <- env_int("INFERCNV_NPCS", 30L)
infercnv_ndims <- env_int("INFERCNV_NDIMS", 30L)
infercnv_resolution <- env_num("INFERCNV_RESOLUTION", 0.5)
infercnv_use_harmony <- env_bool("INFERCNV_USE_HARMONY", TRUE)
infercnv_harmony_group <- env("INFERCNV_HARMONY_GROUP", "orig.ident")

cutoff <- env_num("INFERCNV_CUTOFF", 0.1)
cluster_by_groups <- env_bool("INFERCNV_CLUSTER_BY_GROUPS", FALSE)
analysis_mode <- env("INFERCNV_ANALYSIS_MODE", "samples")
hclust_method <- env("INFERCNV_HCLUST_METHOD", "ward.D2")
denoise <- env_bool("INFERCNV_DENOISE", TRUE)
hmm <- env_bool("INFERCNV_HMM", FALSE)
plot_steps <- env_bool("INFERCNV_PLOT_STEPS", FALSE)
leiden_resolution <- env("INFERCNV_LEIDEN_RESOLUTION", "auto")
num_threads <- env_int("INFERCNV_NUM_THREADS", 6L)
run_median_filter <- env_bool("INFERCNV_MEDIAN_FILTER", FALSE)

log_msg("inferCNV epithelial workflow started")
log_msg("Input QS2: ", input_qs2)
log_msg("Reference cell types: ", paste(infercnv_ref_celltypes, collapse = ","))
log_msg("Epithelial label: ", infercnv_epithelial_label)

sc_all <- read_qs2(input_qs2)
if (!inherits(sc_all, "Seurat")) {
  stop("Input QS2 does not contain a Seurat object: ", input_qs2, call. = FALSE)
}
if (!infercnv_assay %in% names(sc_all@assays)) {
  stop("Missing assay for inferCNV: ", infercnv_assay, call. = FALSE)
}
infercnv_annotation_col <- pick_existing_metadata_col(
  sc_all@meta.data,
  unique(c(infercnv_annotation_col, "cell_type", "sub_cell_type")),
  "annotation"
)

DefaultAssay(sc_all) <- infercnv_assay
sc_all <- join_layers_if_available(sc_all, assay = infercnv_assay)

all_celltypes <- as.character(sc_all[[infercnv_annotation_col, drop = TRUE]])
celltype_by_cell <- setNames(all_celltypes, colnames(sc_all))
infercnv_ref_celltypes <- match_existing_levels(all_celltypes, infercnv_ref_celltypes, "reference cell types")
infercnv_epithelial_label <- match_existing_levels(all_celltypes, infercnv_epithelial_label, "epithelial label")
log_msg("Annotation column: ", infercnv_annotation_col)
log_msg("Matched reference cell types: ", paste(infercnv_ref_celltypes, collapse = ","))
log_msg("Matched epithelial label: ", infercnv_epithelial_label)

ref_cells <- colnames(sc_all)[all_celltypes %in% infercnv_ref_celltypes]
epithelial_cells <- colnames(sc_all)[all_celltypes == infercnv_epithelial_label]

if (length(ref_cells) == 0) {
  stop("No reference cells found for cell types: ", paste(infercnv_ref_celltypes, collapse = ","), call. = FALSE)
}
if (length(epithelial_cells) == 0) {
  stop("No epithelial cells found.", call. = FALSE)
}

log_msg("Reference cells before sampling: ", length(ref_cells))
log_msg("All epithelial cells before sampling: ", length(epithelial_cells))

if (max_ref_cells_per_type > 0) {
  set.seed(SEED)
  ref_cells <- unlist(lapply(split(ref_cells, celltype_by_cell[ref_cells]), function(x) {
    if (length(x) <= max_ref_cells_per_type) return(x)
    sample(x, max_ref_cells_per_type)
  }), use.names = FALSE)
}
if (max_epithelial_cells > 0 && length(epithelial_cells) > max_epithelial_cells) {
  set.seed(SEED)
  epithelial_cells <- sample(epithelial_cells, max_epithelial_cells)
}

log_msg("Reference cells used: ", length(ref_cells))
log_msg("Epithelial cells used for inferCNV: ", length(epithelial_cells))

ref_obj <- subset(sc_all, cells = ref_cells)
ref_obj$infercnv_group <- celltype_by_cell[colnames(ref_obj)]
ref_obj$infercnv_source <- "Reference"

epi_input <- subset(sc_all, cells = epithelial_cells)
epi_input$infercnv_group <- infercnv_observation_label
epi_input$infercnv_source <- "Epithelial"

infercnv_input <- merge(ref_obj, y = epi_input)
infercnv_input <- join_layers_if_available(infercnv_input, assay = infercnv_assay)
infercnv_meta <- infercnv_input@meta.data

counts_mat <- get_counts_matrix(infercnv_input, assay = infercnv_assay)
if (nrow(counts_mat) == 0 || ncol(counts_mat) == 0) {
  stop("Counts matrix is empty.", call. = FALSE)
}

gene_file <- build_gene_order_file(rownames(counts_mat), file.path(output_dir, "gene_order.txt"), species)
gene_order <- read.table(gene_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
ordered_genes <- gene_order[[1]]
ordered_genes <- ordered_genes[ordered_genes %in% rownames(counts_mat)]
counts_mat <- counts_mat[ordered_genes, , drop = FALSE]

annotation_file <- file.path(output_dir, "cell_annotations.txt")
annotation_df <- data.frame(
  cell = colnames(counts_mat),
  group = infercnv_meta[colnames(counts_mat), "infercnv_group"]
)
write.table(annotation_df, annotation_file, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.csv(as.data.frame(table(annotation_df$group)), file.path(output_dir, "cell_group_counts.csv"), row.names = FALSE)

log_msg("inferCNV matrix: genes=", nrow(counts_mat), "; cells=", ncol(counts_mat))

infercnv_obj <- run_step(
  "CreateInfercnvObject",
  CreateInfercnvObject(
    raw_counts_matrix = counts_mat,
    annotations_file = annotation_file,
    delim = "\t",
    gene_order_file = gene_file,
    ref_group_names = infercnv_ref_celltypes
  )
)

infercnv_obj <- run_step(
  "infercnv run",
  infercnv_run(
    infercnv_obj,
    cutoff = cutoff,
    out_dir = infercnv_output_dir,
    cluster_by_groups = cluster_by_groups,
    hclust_method = hclust_method,
    analysis_mode = analysis_mode,
    denoise = denoise,
    HMM = hmm,
    plot_steps = plot_steps,
    leiden_resolution = leiden_resolution,
    num_threads = num_threads
  )
)

save_qs2(infercnv_obj, file.path(output_dir, "infercnv_obj.qs2"))

if (run_median_filter) {
  infercnv_obj_medianfiltered <- run_step("infercnv median filtering", apply_median_filtering(infercnv_obj))
  save_qs2(infercnv_obj_medianfiltered, file.path(output_dir, "infercnv_obj_median_filtered.qs2"))
  infercnv_for_score <- infercnv_obj_medianfiltered
} else {
  infercnv_for_score <- infercnv_obj
}

cnv_scores <- calculate_cnv_score(infercnv_for_score, infercnv_ref_celltypes)

all_epi <- subset(sc_all, cells = epithelial_cells)
all_epi <- run_step("Cluster all epithelial cells", run_epithelial_clustering(all_epi))
all_epi$infercnv_group <- paste0("Epi_", as.character(all_epi$seurat_clusters))
all_epi$infercnv_source <- "Epithelial"

save_plot(
  DimPlot(all_epi, group.by = "infercnv_group", label = TRUE, repel = TRUE),
  file.path(output_dir, "all_epithelial_clusters"),
  8,
  6
)

score_df <- data.frame(
  cell = names(cnv_scores),
  totalCNV = as.numeric(cnv_scores),
  infercnv_input_group = annotation_df$group[match(names(cnv_scores), annotation_df$cell)]
)
cluster_df <- data.frame(
  cell = colnames(all_epi),
  epithelial_cluster = as.character(all_epi$infercnv_group)
)
score_df$infercnv_group <- score_df$infercnv_input_group
score_df$infercnv_group[score_df$infercnv_input_group %in% infercnv_ref_celltypes] <- infercnv_ref_label
cluster_match <- match(score_df$cell, cluster_df$cell)
score_df$infercnv_group[!is.na(cluster_match)] <- cluster_df$epithelial_cluster[cluster_match[!is.na(cluster_match)]]
score_df <- score_df[!is.na(score_df$infercnv_group), ]
score_df$infercnv_group <- factor(
  score_df$infercnv_group,
  levels = c(sort(setdiff(unique(score_df$infercnv_group), infercnv_ref_label)), infercnv_ref_label)
)
write.csv(score_df, file.path(output_dir, "infercnv_cell_cnv_scores.csv"), row.names = FALSE)

all_epi$totalCNV <- NA_real_
common_cells <- intersect(colnames(all_epi), score_df$cell)
all_epi$totalCNV[common_cells] <- score_df$totalCNV[match(common_cells, score_df$cell)]
save_qs2(all_epi, file.path(output_dir, "all_epithelial_clustered_with_cnv_score.qs2"))

ref_scores <- score_df$totalCNV[score_df$infercnv_group == infercnv_ref_label]
stat_df <- data.frame()
for (group_name in setdiff(levels(score_df$infercnv_group), infercnv_ref_label)) {
  test_scores <- score_df$totalCNV[score_df$infercnv_group == group_name]
  if (length(test_scores) > 0 && length(ref_scores) > 0) {
    p_value <- wilcox.test(test_scores, ref_scores)$p.value
    stat_df <- rbind(
      stat_df,
      data.frame(
        group = group_name,
        median_cnv = median(test_scores),
        ref_median_cnv = median(ref_scores),
        p_value = p_value,
        p_adj = NA_real_
      )
    )
  }
}
if (nrow(stat_df) > 0) {
  stat_df$p_adj <- p.adjust(stat_df$p_value, method = "BH")
}
write.csv(stat_df, file.path(output_dir, "infercnv_group_wilcox_vs_ref.csv"), row.names = FALSE)

p <- ggplot(score_df, aes(x = infercnv_group, y = totalCNV, color = infercnv_group, fill = infercnv_group)) +
  geom_boxplot(alpha = 0.45, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.15, alpha = 0.35) +
  labs(x = NULL, y = "CNV_score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
save_plot(p, file.path(output_dir, "infercnv_cnv_score_boxplot"), 10, 6)

if (run_median_filter) {
  run_step(
    "infercnv median filtered plot",
    plot_cnv(
      infercnv_for_score,
      out_dir = file.path(output_dir, "infercnv_median_filtered"),
      output_filename = "infercnv.median_filtered",
      x.range = "auto",
      x.center = 1,
      title = "infercnv",
      color_safe_pal = FALSE
    )
  )
}

log_msg("inferCNV epithelial workflow done.")
