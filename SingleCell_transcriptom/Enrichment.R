rm(list=ls())

## =========================
## libraries
## =========================
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(RColorBrewer)

## =========================
## plotting utils
## =========================
save_plot <- function(p, filename_no_ext, width = 10, height = 6, dpi = 300) {
  ggsave(paste0(filename_no_ext, ".pdf"),
         plot = p, width = width, height = height, dpi = dpi)
  ggsave(paste0(filename_no_ext, ".png"),
         plot = p, width = width, height = height, dpi = dpi)
}

## =========================
## gene id conversion
## =========================
convert_deg <- function(deg, org_db = org.Mm.eg.db) {
  gene_df <- bitr(
    deg$gene,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org_db
  )
  inner_join(deg, gene_df, by = c("gene" = "SYMBOL"))
}

## =========================
## enrichment + plot
## =========================
perform_enrichment_analysis <- function(
  deg,
  out_dir,
  label,
  enrichment_type = "GO",
  org_db = org.Mm.eg.db
) {

  if (nrow(deg) < 10) return(NULL)

  deg2 <- convert_deg(deg, org_db)
  if (nrow(deg2) < 10) return(NULL)

  if (enrichment_type == "GO") {
    enr <- enrichGO(
      gene     = deg2$ENTREZID,
      OrgDb   = org_db,
      ont     = "BP",
      readable = TRUE
    )
    file_name <- "GO"
  } else {
    enr <- enrichKEGG(
      gene     = deg2$ENTREZID,
      organism = "mmu"
    )
    file_name <- "KEGG"
  }

  enr_df <- as.data.frame(enr)
  if (nrow(enr_df) == 0) return(NULL)

  write.csv(
    enr_df,
    file.path(out_dir, paste0(file_name, "_result.csv")),
    row.names = FALSE
  )

  plot_df <- enr_df[1:min(15, nrow(enr_df)), ]

  p <- ggplot(
    plot_df,
    aes(1, reorder(Description, -pvalue))
  ) +
    geom_point(
      aes(size = Count, color = -log2(pvalue))
    ) +
    scale_colour_gradient(low = "#059CFE", high = "#FE05EB") +
    scale_size(range = c(3, 20)) +
    labs(
      title = label,
      x = "", y = "",
      color = "-log2(p)",
      size  = "Gene number"
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks  = element_blank(),
      panel.grid  = element_blank(),
      panel.background = element_blank(),
      axis.text = element_text(size = 12, face = "bold"),
      legend.position = "bottom"
    ) +
    scale_y_discrete(
      position = "right",
      labels = function(x) stringr::str_wrap(x, 40)
    )

  save_plot(
    p,
    file.path(out_dir, file_name),
    height = 10
  )
}

## =========================
## parameterized FindMarkers + enrichment
## =========================
perform_analysis <- function(
  seu,
  cell_type_col,
  targets = NULL,          # NULL = all subtypes
  group_col = "group",
  comparisons,             # list(c("T1","T2"), ...)
  outdir = "06_Enrich"
) {

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  ## set cell identity
  Idents(seu) <- seu[[cell_type_col]][, 1]

  if (is.null(targets)) {
    targets <- levels(Idents(seu))
  }

  ## loop cell subtypes
  for (ct in targets) {

    message("Cell subtype: ", ct)

    ## loop comparisons
    for (cmp in comparisons) {

      ctrl <- cmp[1]
      test <- cmp[2]
      cmp_name <- paste0(test, "_vs_", ctrl)

      message("  └─ ", cmp_name)

      out_dir <- file.path(outdir, make.names(ct), cmp_name)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      ## differential expression
      deg <- FindMarkers(
        seu,
        subset.ident = ct,
        group.by = group_col,
        ident.1 = test,
        ident.2 = ctrl,
        only.pos = TRUE,
        logfc.threshold = 0.25,
        min.pct = 0.1
      ) %>%
        rownames_to_column("gene") %>%
        filter(p_val_adj < 0.05)

      if (nrow(deg) < 10) {
        message("    skip (too few DEGs)")
        next
      }

      write.csv(
        deg,
        file.path(out_dir, "DEG_FindMarkers.csv"),
        row.names = FALSE
      )

      ## enrichment
      perform_enrichment_analysis(
        deg,
        out_dir,
        label = paste(ct, cmp_name, "GO"),
        enrichment_type = "GO"
      )

      perform_enrichment_analysis(
        deg,
        out_dir,
        label = paste(ct, cmp_name, "KEGG"),
        enrichment_type = "KEGG"
      )
    }
  }
}

## =========================
## run
## =========================
seu <- readRDS("03_Subcluster\\Epithelial\\subcluster_annotated.rds")
seu$sub_cell_type
cell_type_col <- "sub_cell_type" # 主群是cell_type
targets <- names(table(seu$sub_cell_type))
comparisons <- list(
  c("T1", "T2"),
  c("CON", "T1"),
  c("CON", "T2")
)
Sys.setenv(http_proxy = "http://127.0.0.1:7897")
Sys.setenv(https_proxy = "http://127.0.0.1:7897")
perform_analysis(
  seu = seu,
  cell_type_col = cell_type_col,
  targets = targets,
  group_col = "group",
  comparisons = comparisons,
  outdir = "06_Enrich"
)
