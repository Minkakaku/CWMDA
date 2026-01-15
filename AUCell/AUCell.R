############################################################
## AUCell (Tpex / Tex terminal)
## Visualization: scplotter::FeatureStatPlot
## FULL pipeline from cd8_clean_final.rds
############################################################

rm(list = ls())
set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(AUCell)
  library(Matrix)
  library(scplotter)
  library(ggplot2)
  library(patchwork)
})

## ===============================
## 0) Paths
## ===============================
out_dir <- "output"

auc_dir <- file.path(out_dir, "03_CD8_Subcluster/07_CD8_AUCell")
auc_rds <- file.path(auc_dir, "rds")
auc_fig <- file.path(auc_dir, "figures")

dir.create(auc_rds, recursive = TRUE, showWarnings = FALSE)
dir.create(auc_fig, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, filename, w = 10, h = 5, dpi = 300) {
  ggsave(paste0(filename, ".png"), p, width = w, height = h, dpi = dpi)
  ggsave(paste0(filename, ".pdf"), p, width = w, height = h)
}

## ===============================
## 1) Load clean CD8 object
## ===============================
seu_path <- file.path(out_dir, "rds/cd8_clean_final.rds")
stopifnot(file.exists(seu_path))

seu <- readRDS(seu_path)
FeaturePlot(seu,c("Cd8a","Cd8"))

stopifnot(inherits(seu, "Seurat"))
stopifnot("umap" %in% names(seu@reductions))

## Required meta columns
stopifnot("cd8_cell_type" %in% colnames(seu@meta.data))

## Optional split column
split_col <- if ("Phase" %in% colnames(seu@meta.data)) {
  "Phase"
} else if ("group" %in% colnames(seu@meta.data)) {
  "group"
} else {
  NULL
}

## ===============================
## 2) Expression matrix for AUCell
## ===============================
assay_use <- DefaultAssay(seu)
slot_use  <- "data"

expr <- tryCatch(
  GetAssayData(seu, assay = assay_use, layer = slot_use),
  error = function(e) GetAssayData(seu, assay = assay_use, slot = slot_use)
)

if (!inherits(expr, "dgCMatrix")) {
  expr <- as(as.matrix(expr), "dgCMatrix")
}

genes_in_data <- rownames(expr)

## ===============================
## 3) Gene sets (mouse, CD8-focused)
## ===============================
Tpex_geneset <- c(
  "Tcf7","Slamf6","Pdcd1","Tox","Ikzf2","Myb",
  "Bach2","Lef1","Il7r","Cxcr5","Nr4a2","Xcl1","Sell","Ltb"
)

Tex_terminal_geneset <- c(
  "Pdcd1","Havcr2","Lag3","Tigit","Ctla4","Entpd1",
  "Tox","Tox2","Prdm1","Batf","Eomes",
  "Cxcl13","Nr4a1","Nr4a2","Nr4a3"
)

match_genes <- function(gs, genes) {
  map <- genes
  names(map) <- toupper(genes)
  hit <- map[toupper(gs)]
  unique(unname(hit[!is.na(hit)]))
}

geneSets <- list(
  Tpex        = match_genes(Tpex_geneset, genes_in_data),
  TexTerminal = match_genes(Tex_terminal_geneset, genes_in_data)
)

sizes <- sapply(geneSets, length)
print(sizes)
if (any(sizes < 5)) {
  warning("Some gene sets have <5 genes after filtering.")
}

## ===============================
## 4) AUCell
## ===============================
rankings <- AUCell_buildRankings(
  expr,
  nCores = max(1, parallel::detectCores() - 1),
  plotStats = FALSE
)

cellsAUC <- AUCell_calcAUC(
  geneSets,
  rankings,
  aucMaxRank = ceiling(0.05 * nrow(expr))
)

auc_df <- as.data.frame(t(getAUC(cellsAUC)))
colnames(auc_df) <- paste0("AUC_", colnames(auc_df))
auc_df <- auc_df[colnames(seu), , drop = FALSE]

seu <- AddMetaData(seu, auc_df)

## ===============================
## 5) FeatureStatPlot — Tpex
## ===============================
p_tpex_dim <- FeatureStatPlot(
  object    = seu,
  plot_type = "dim",
  features  = "AUC_Tpex",
  reduction = "umap",
  split_by  = split_col,
  highlight = TRUE,
  order     = "high-top",
  palette   = "Blues",
  theme     = "theme_blank",
  title     = "Tpex AUCell score (UMAP)"
)

p_tpex_vln <- FeatureStatPlot(
  object    = seu,
  plot_type = "violin",
  features  = "AUC_Tpex",
  ident   = "cd8_cell_type",
  split_by  = split_col,
  palette   = "Blues",
  title     = "Tpex AUCell score by CD8 cell type"
)

p_tpex <- p_tpex_dim / p_tpex_vln
p_tpex
## ===============================
## 6) FeatureStatPlot — Tex terminal
## ===============================
p_tex_dim <- FeatureStatPlot(
  object    = seu,
  plot_type = "dim",
  features  = "AUC_TexTerminal",
  reduction = "umap",
  split_by  = split_col,
  highlight = TRUE,
  order     = "high-top",
  palette   = "Reds",
  theme     = "theme_blank",
  title     = "Tex terminal AUCell score (UMAP)"
)

p_tex_vln <- FeatureStatPlot(
  object    = seu,
  plot_type = "violin",
  features  = "AUC_TexTerminal",
  ident   = "cd8_cell_type",
  split_by  = split_col,
  palette   = "Reds",
  title     = "Tex terminal AUCell score by CD8 cell type"
)

p_tex <- p_tex_dim / p_tex_vln

## ===============================
## 7) Combine & save
## ===============================

save_plot(
  p_tpex,
  file.path(auc_fig, "FeatureStatPlot_Tpex_AUCell"),
  w = 10, h = 10
)

save_plot(
  p_tex,
  file.path(auc_fig, "FeatureStatPlot_TexTerminal_AUCell"),
  w = 10, h = 10
)


saveRDS(
  seu,
  file = file.path(auc_rds, "cd8_clean_final_with_Tpex_Tex_AUCell.rds")
)

message("DONE: AUCell + FeatureStatPlot")
message("Figures saved to: ", auc_fig)
message("RDS saved to: ", auc_rds)
