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
  library(nichenetr)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

resolve_resource_file <- function(explicit_path, resources_dir, file_name) {
  if (!is.null(explicit_path) && explicit_path != "") {
    return(explicit_path)
  }
  if (!is.null(resources_dir) && resources_dir != "") {
    return(file.path(resources_dir, file_name))
  }
  NULL
}

load_nichenetr_resources <- function(organism, resources_dir, lr_network_path, ligand_target_path, weighted_networks_path) {
  resource_names <- if (organism == "human") {
    list(
      lr = "lr_network_human_21122021.rds",
      lt = "ligand_target_matrix_nsga2r_final.rds",
      wn = "weighted_networks_nsga2r_final.rds"
    )
  } else {
    list(
      lr = "lr_network_mouse_21122021.rds",
      lt = "ligand_target_matrix_nsga2r_final_mouse.rds",
      wn = "weighted_networks_nsga2r_final_mouse.rds"
    )
  }

  lr_file <- resolve_resource_file(lr_network_path, resources_dir, resource_names$lr)
  lt_file <- resolve_resource_file(ligand_target_path, resources_dir, resource_names$lt)
  wn_file <- resolve_resource_file(weighted_networks_path, resources_dir, resource_names$wn)

  if (is.null(lr_file) || is.null(lt_file) || is.null(wn_file)) {
    stop("Please provide NicheNet resource files with --resources-dir or explicit --lr-network/--ligand-target-matrix/--weighted-networks.")
  }

  list(
    lr_network = readRDS(lr_file) %>% distinct(from, to),
    ligand_target_matrix = readRDS(lt_file),
    weighted_networks = readRDS(wn_file)
  )
}

get_logfc_column <- function(deg_df) {
  if ("avg_log2FC" %in% colnames(deg_df)) {
    return("avg_log2FC")
  }
  if ("avg_logFC" %in% colnames(deg_df)) {
    return("avg_logFC")
  }
  stop("Cannot find logFC column in FindMarkers result.")
}

save_gg <- function(plot_obj, path_prefix, width = 10, height = 8) {
  ggplot2::ggsave(paste0(path_prefix, ".png"), plot_obj, width = width, height = height)
  ggplot2::ggsave(paste0(path_prefix, ".pdf"), plot_obj, width = width, height = height)
}

args <- parse_cli_args()
input_rds <- cli_string(args, "input-rds")
out_dir <- ensure_dir(cli_string(args, "out-dir", "nichenetr_result"))
organism <- cli_string(args, "organism", "human")
resources_dir <- cli_string(args, "resources-dir")
lr_network_path <- cli_string(args, "lr-network")
ligand_target_path <- cli_string(args, "ligand-target-matrix")
weighted_networks_path <- cli_string(args, "weighted-networks")
assay_name <- cli_string(args, "assay", "RNA")
convert_assay <- cli_bool(args, "convert-assay", FALSE)
celltype_col <- cli_string(args, "celltype-col", "cell_type")
condition_col <- cli_string(args, "condition-col", "group")
receiver <- cli_string(args, "receiver")
sender_celltypes <- cli_csv(args, "sender-celltypes")
condition_oi <- cli_string(args, "condition-oi")
condition_reference <- cli_string(args, "condition-reference")
pct <- cli_numeric(args, "pct", 0.05)
lfc_threshold <- cli_numeric(args, "lfc-threshold", 0.25)
p_adj_threshold <- cli_numeric(args, "p-adj-threshold", 0.05)
top_ligands <- cli_numeric(args, "top-ligands", 30)

if (is.null(input_rds)) {
  stop("Please provide --input-rds.")
}
if (is.null(receiver)) {
  stop("Please provide --receiver.")
}
if (length(sender_celltypes) == 0) {
  stop("Please provide --sender-celltypes.")
}
if (is.null(condition_oi) || is.null(condition_reference)) {
  stop("Please provide --condition-oi and --condition-reference.")
}

organism <- match.arg(organism, c("human", "mouse"))
resources <- load_nichenetr_resources(
  organism = organism,
  resources_dir = resources_dir,
  lr_network_path = lr_network_path,
  ligand_target_path = ligand_target_path,
  weighted_networks_path = weighted_networks_path
)

seu <- readRDS(input_rds)
DefaultAssay(seu) <- assay_name
seu <- join_layers_if_needed(seu)

if (convert_assay && exists("Convert_Assay", mode = "function")) {
  seu <- Convert_Assay(seurat_object = seu, convert_to = "V3")
}

if (!celltype_col %in% colnames(seu@meta.data)) {
  stop("celltype column not found: ", celltype_col)
}
if (!condition_col %in% colnames(seu@meta.data)) {
  stop("condition column not found: ", condition_col)
}

Idents(seu) <- seu[[celltype_col, drop = TRUE]]

if (!receiver %in% levels(Idents(seu))) {
  stop("receiver cell type not found in ", celltype_col, ": ", receiver)
}
if (!all(sender_celltypes %in% levels(Idents(seu)))) {
  missing_sender <- sender_celltypes[!sender_celltypes %in% levels(Idents(seu))]
  stop("sender cell types not found in ", celltype_col, ": ", paste(missing_sender, collapse = ", "))
}

expressed_genes_receiver <- get_expressed_genes(receiver, seu, pct = pct)
all_receptors <- unique(resources$lr_network$to)
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- resources$lr_network %>%
  filter(to %in% expressed_receptors) %>%
  pull(from) %>%
  unique()

if (length(expressed_receptors) == 0 || length(potential_ligands) == 0) {
  stop("No expressed receptors or candidate ligands found for receiver: ", receiver)
}

sender_celltypes <- unique(sender_celltypes)
expressed_genes_sender <- sender_celltypes %>%
  lapply(function(celltype_name) get_expressed_genes(celltype_name, seu, pct = pct)) %>%
  unlist() %>%
  unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender)

seu_receiver <- subset(seu, idents = receiver)

deg_receiver <- FindMarkers(
  object = seu_receiver,
  ident.1 = condition_oi,
  ident.2 = condition_reference,
  group.by = condition_col,
  min.pct = pct
) %>%
  tibble::rownames_to_column("gene")

lfc_col <- get_logfc_column(deg_receiver)

write.table(
  deg_receiver,
  file.path(out_dir, "receiver_deg.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

geneset_oi <- deg_receiver %>%
  filter(p_val_adj <= p_adj_threshold & abs(.data[[lfc_col]]) >= lfc_threshold) %>%
  pull(gene)
geneset_oi <- intersect(geneset_oi, rownames(resources$ligand_target_matrix))

if (length(geneset_oi) == 0) {
  stop("No receiver DE genes passed the filtering thresholds.")
}

background_expressed_genes <- intersect(
  expressed_genes_receiver,
  rownames(resources$ligand_target_matrix)
)

ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = resources$ligand_target_matrix,
  potential_ligands = potential_ligands
) %>%
  arrange(desc(aupr_corrected)) %>%
  mutate(rank = rank(desc(aupr_corrected)))

write.table(
  ligand_activities,
  file.path(out_dir, "ligand_activities.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_hist <- ggplot(ligand_activities, aes(x = aupr_corrected)) +
  geom_histogram(color = "black", fill = "darkorange") +
  geom_vline(
    xintercept = min(dplyr::slice_max(ligand_activities, order_by = aupr_corrected, n = top_ligands, with_ties = FALSE)$aupr_corrected),
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  labs(x = paste0("Ligand activity (", condition_oi, ")"), y = "# ligands") +
  theme_classic()
save_gg(p_hist, file.path(out_dir, "ligand_activity_hist"), width = 8, height = 6)

best_upstream_ligands <- ligand_activities %>%
  dplyr::slice_max(order_by = aupr_corrected, n = top_ligands, with_ties = FALSE) %>%
  arrange(desc(aupr_corrected)) %>%
  pull(test_ligand) %>%
  unique()

write.table(
  data.frame(best_upstream_ligands = best_upstream_ligands),
  file.path(out_dir, "top_upstream_ligands.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

plot_ligand_activity <- function(ligand_activity_df, ligand_names) {
  vis_ligand_aupr <- ligand_activity_df %>%
    filter(test_ligand %in% ligand_names) %>%
    tibble::column_to_rownames("test_ligand") %>%
    select(aupr_corrected) %>%
    arrange(aupr_corrected) %>%
    as.matrix()

  make_heatmap_ggplot(
    vis_ligand_aupr,
    "Prioritized ligands",
    "Ligand activity",
    legend_title = "AUPR",
    color = "darkorange"
  ) + theme(axis.text.x.top = element_blank())
}

plot_regulatory_heatmap <- function(link_df, ligand_target_matrix, ligand_names) {
  active_links <- prepare_ligand_target_visualization(
    ligand_target_df = link_df,
    ligand_target_matrix = ligand_target_matrix,
    cutoff = 0.33
  )
  order_ligands <- intersect(ligand_names, colnames(active_links)) %>% rev()
  order_targets <- intersect(unique(link_df$target), rownames(active_links))
  vis_ligand_target <- t(active_links[order_targets, order_ligands, drop = FALSE])

  make_heatmap_ggplot(
    vis_ligand_target,
    "Prioritized ligands",
    "Predicted target genes",
    color = "purple",
    legend_title = "Regulatory potential"
  ) +
    scale_fill_gradient2(low = "whitesmoke", high = "purple")
}

plot_receptor_heatmap <- function(link_df, ligand_names) {
  vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
    link_df,
    ligand_names,
    order_hclust = "both"
  )
  make_heatmap_ggplot(
    t(vis_ligand_receptor_network),
    y_name = "Ligands",
    x_name = "Receptors",
    color = "mediumvioletred",
    legend_title = "Prior interaction potential"
  )
}

p_ligand_activity <- plot_ligand_activity(ligand_activities, best_upstream_ligands)
save_gg(p_ligand_activity, file.path(out_dir, "ligand_activity"), width = 10, height = 8)

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(
    get_weighted_ligand_target_links,
    geneset = geneset_oi,
    ligand_target_matrix = resources$ligand_target_matrix,
    n = 100
  ) %>%
  bind_rows() %>%
  tidyr::drop_na()

write.table(
  active_ligand_target_links_df,
  file.path(out_dir, "active_ligand_target_links.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_target <- plot_regulatory_heatmap(
  link_df = active_ligand_target_links_df,
  ligand_target_matrix = resources$ligand_target_matrix,
  ligand_names = best_upstream_ligands
)
save_gg(p_target, file.path(out_dir, "predicted_target_genes"), width = 16, height = 10)

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands,
  expressed_receptors,
  resources$lr_network,
  resources$weighted_networks$lr_sig
)

write.table(
  ligand_receptor_links_df,
  file.path(out_dir, "ligand_receptor_links.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_receptor <- plot_receptor_heatmap(ligand_receptor_links_df, best_upstream_ligands)
save_gg(p_receptor, file.path(out_dir, "prior_interaction_potential"), width = 12, height = 8)

if (length(potential_ligands_focused) > 0) {
  ligand_activities_sender <- ligand_activities %>%
    filter(test_ligand %in% potential_ligands_focused)

  best_sender_ligands <- ligand_activities_sender %>%
    dplyr::slice_max(order_by = aupr_corrected, n = top_ligands, with_ties = FALSE) %>%
    arrange(desc(aupr_corrected)) %>%
    pull(test_ligand) %>%
    unique()

  p_sender_activity <- plot_ligand_activity(ligand_activities_sender, best_sender_ligands)
  save_gg(p_sender_activity, file.path(out_dir, "sender_focused_ligand_activity"), width = 10, height = 8)

  active_sender_target_links_df <- best_sender_ligands %>%
    lapply(
      get_weighted_ligand_target_links,
      geneset = geneset_oi,
      ligand_target_matrix = resources$ligand_target_matrix,
      n = 100
    ) %>%
    bind_rows() %>%
    tidyr::drop_na()

  p_sender_target <- plot_regulatory_heatmap(
    link_df = active_sender_target_links_df,
    ligand_target_matrix = resources$ligand_target_matrix,
    ligand_names = best_sender_ligands
  )
  save_gg(p_sender_target, file.path(out_dir, "sender_focused_target_genes"), width = 16, height = 10)

  sender_receptor_links_df <- get_weighted_ligand_receptor_links(
    best_sender_ligands,
    expressed_receptors,
    resources$lr_network,
    resources$weighted_networks$lr_sig
  )
  p_sender_receptor <- plot_receptor_heatmap(sender_receptor_links_df, best_sender_ligands)
  save_gg(p_sender_receptor, file.path(out_dir, "sender_focused_receptors"), width = 12, height = 8)

  sender_cells <- colnames(seu)[seu[[celltype_col, drop = TRUE]] %in% sender_celltypes]
  p_dotplot <- DotPlot(
    subset(seu, cells = sender_cells),
    features = rev(best_sender_ligands)
  ) +
    coord_flip()
  save_gg(p_dotplot, file.path(out_dir, "sender_focused_dotplot"), width = 12, height = 8)

  celltype_order <- levels(Idents(seu))
  de_table_top_ligands <- lapply(
    celltype_order[celltype_order %in% sender_celltypes],
    get_lfc_celltype,
    seurat_obj = seu,
    condition_colname = condition_col,
    condition_oi = condition_oi,
    condition_reference = condition_reference,
    celltype_col = celltype_col,
    min.pct = 0,
    logfc.threshold = 0,
    features = best_sender_ligands
  )

  if (length(de_table_top_ligands) > 0) {
    de_table_top_ligands <- Reduce(f = full_join, x = de_table_top_ligands) %>%
      tibble::column_to_rownames("gene")

    vis_ligand_lfc <- as.matrix(de_table_top_ligands[rev(best_sender_ligands), , drop = FALSE])
    p_lfc <- make_threecolor_heatmap_ggplot(
      vis_ligand_lfc,
      "Prioritized ligands",
      "LFC in sender",
      low_color = "midnightblue",
      mid_color = "white",
      mid = stats::median(vis_ligand_lfc),
      high_color = "red",
      legend_title = "LFC"
    )
    save_gg(p_lfc, file.path(out_dir, "sender_focused_lfc"), width = 12, height = 8)
  }
}

message("NicheNet finished: ", out_dir)
