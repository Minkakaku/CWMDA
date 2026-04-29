rm(list = ls())

get_script_dir_bootstrap <- function() {
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

script_dir <- get_script_dir_bootstrap()
source(file.path(script_dir, "common_utils.R"), local = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(dplyr)
})

args <- parse_cli_args()
input_rds <- cli_string(args, "input-rds")
out_dir <- ensure_dir(cli_string(args, "out-dir", "cellchat_result"))
group_col <- cli_string(args, "group-col", "group")
celltype_col <- cli_string(args, "celltype-col", "cell_type")
assay_name <- cli_string(args, "assay", "RNA")
species <- cli_string(args, "species", "human")
min_cells <- cli_numeric(args, "min-cells", 10)
group_order <- cli_csv(args, "group-order")
subset_col <- cli_string(args, "subset-col")
subset_values <- cli_csv(args, "subset-values")

if (is.null(input_rds)) {
  stop("Please provide --input-rds.")
}

species <- match.arg(species, c("human", "mouse"))
seu <- readRDS(input_rds)

if (!is.null(subset_col) && length(subset_values) > 0) {
  if (!subset_col %in% colnames(seu@meta.data)) {
    stop("subset column not found: ", subset_col)
  }
  keep_cells <- colnames(seu)[seu[[subset_col, drop = TRUE]] %in% subset_values]
  seu <- subset(seu, cells = keep_cells)
}

if (!group_col %in% colnames(seu@meta.data)) {
  stop("group column not found: ", group_col)
}
if (!celltype_col %in% colnames(seu@meta.data)) {
  stop("celltype column not found: ", celltype_col)
}

if (length(group_order) > 0) {
  seu[[group_col]] <- factor(seu[[group_col, drop = TRUE]], levels = group_order)
}

DefaultAssay(seu) <- assay_name
seu <- NormalizeData(seu, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)

get_data_matrix <- function(obj, assay_name) {
  DefaultAssay(obj) <- assay_name
  tryCatch(
    GetAssayData(obj, layer = "data"),
    error = function(e) GetAssayData(obj, slot = "data")
  )
}

run_cellchat_one <- function(obj, dataset_name) {
  obj <- join_layers_if_needed(obj)
  data_input <- get_data_matrix(obj, assay_name)
  meta <- obj@meta.data[, c(group_col, celltype_col), drop = FALSE]
  colnames(meta) <- c("Group", "CellTypes")

  cellchat <- createCellChat(
    object = data_input,
    meta = meta,
    group.by = "CellTypes"
  )

  if (species == "mouse") {
    cellchat@DB <- CellChatDB.mouse
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- smoothData(cellchat, adj = PPI.mouse)
  } else {
    cellchat@DB <- CellChatDB.human
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- smoothData(cellchat, adj = PPI.human)
  }

  cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  cellchat <- computeNetSimilarity(cellchat, type = "functional")

  dataset_prefix <- file.path(out_dir, dataset_name)
  write.table(
    subsetCommunication(cellchat),
    paste0(dataset_prefix, "_net_lr.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  write.table(
    subsetCommunication(cellchat, slot.name = "netP"),
    paste0(dataset_prefix, "_pathway_lr.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  saveRDS(cellchat, paste0(dataset_prefix, ".rds"))
  cellchat
}

group_levels <- unique(as.character(seu[[group_col, drop = TRUE]]))
group_levels <- group_levels[!is.na(group_levels)]
if (length(group_levels) < 2) {
  stop("CellChat requires at least two groups after filtering.")
}

object_list <- lapply(group_levels, function(level_name) {
  subset(seu, cells = colnames(seu)[seu[[group_col, drop = TRUE]] == level_name])
})
names(object_list) <- group_levels

cellchat_list <- lapply(names(object_list), function(dataset_name) {
  run_cellchat_one(object_list[[dataset_name]], dataset_name)
})
names(cellchat_list) <- group_levels

cellchat_merged <- mergeCellChat(
  cellchat_list,
  add.names = names(cellchat_list),
  cell.prefix = TRUE
)

saveRDS(cellchat_merged, file.path(out_dir, "cellchat_merged.rds"))

count_plot <- compareInteractions(cellchat_merged, show.legend = FALSE)
save_plot(count_plot, file.path(out_dir, "compare_interactions_count"), width = 7, height = 5)

weight_plot <- compareInteractions(cellchat_merged, show.legend = FALSE, measure = "weight")
save_plot(weight_plot, file.path(out_dir, "compare_interactions_weight"), width = 7, height = 5)

weight_max_count <- getMaxWeight(cellchat_list, attribute = c("idents", "count"))
pdf(file.path(out_dir, "circle_count_per_dataset.pdf"), width = max(8, 6 * length(cellchat_list)), height = 6)
par(mfrow = c(1, length(cellchat_list)), xpd = TRUE)
for (idx in seq_along(cellchat_list)) {
  netVisual_circle(
    cellchat_list[[idx]]@net$count,
    weight.scale = TRUE,
    label.edge = FALSE,
    edge.weight.max = weight_max_count[2],
    edge.width.max = 12,
    title.name = paste0("Interaction count - ", names(cellchat_list)[idx])
  )
}
dev.off()

weight_max_weight <- getMaxWeight(cellchat_list, attribute = c("idents", "weight"))
pdf(file.path(out_dir, "circle_weight_per_dataset.pdf"), width = max(8, 6 * length(cellchat_list)), height = 6)
par(mfrow = c(1, length(cellchat_list)), xpd = TRUE)
for (idx in seq_along(cellchat_list)) {
  netVisual_circle(
    cellchat_list[[idx]]@net$weight,
    weight.scale = TRUE,
    label.edge = FALSE,
    edge.weight.max = weight_max_weight[2],
    edge.width.max = 12,
    title.name = paste0("Interaction weight - ", names(cellchat_list)[idx])
  )
}
dev.off()

if (length(cellchat_list) == 2) {
  pdf(file.path(out_dir, "diff_interaction_circle.pdf"), width = 12, height = 6)
  par(mfrow = c(1, 2), xpd = TRUE)
  netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE)
  netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "weight")
  dev.off()

  diff_heatmap_count <- netVisual_heatmap(cellchat_merged)
  diff_heatmap_weight <- netVisual_heatmap(cellchat_merged, measure = "weight")
  pdf(file.path(out_dir, "diff_interaction_heatmap.pdf"), width = 10, height = 5)
  print(diff_heatmap_count + diff_heatmap_weight)
  dev.off()
}

pathway_union <- lapply(cellchat_list, function(x) unique(x@netP$pathways)) |>
  unlist() |>
  unique()

for (pattern_name in c("incoming", "outgoing", "all")) {
  heatmaps <- lapply(seq_along(cellchat_list), function(idx) {
    netAnalysis_signalingRole_heatmap(
      cellchat_list[[idx]],
      pattern = pattern_name,
      signaling = pathway_union,
      title = names(cellchat_list)[idx],
      width = 8,
      height = max(6, length(pathway_union) / 2)
    )
  })
  pdf(
    file.path(out_dir, paste0("signaling_role_", pattern_name, ".pdf")),
    width = 8 * length(heatmaps),
    height = max(6, length(pathway_union) / 2)
  )
  print(Reduce(`+`, heatmaps))
  dev.off()
}

message("CellChat finished: ", out_dir)
