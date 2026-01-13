############################################################
## AUCell (Tpex / Tex terminal)
## Visualization: scplotter::FeatureStatPlot
## FULL pipeline from cd8_clean_final.rds
############################################################

rm(list = ls())
set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(scplotter)
  library(ggplot2)
  library(patchwork)
})

## ===============================
## 0) Paths
## ===============================
out_dir <- "output"

auc_dir <- file.path(out_dir, "03_CD8_Subcluster/08_target_gene_show")
auc_rds <- file.path(auc_dir, "rds")
auc_fig <- file.path(auc_dir, "figures")

dir.create(auc_rds, recursive = TRUE, showWarnings = FALSE)
dir.create(auc_fig, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, filename, w = 10, h = 5, dpi = 300) {
  ggsave(paste0(filename, ".png"), p, width = w, height = h, dpi = dpi)
  ggsave(paste0(filename, ".pdf"), p, width = w, height = h)
}

sc = readRDS("output/rds/sc_global_annotations.rds")

sc$group <- NA_character_
sc$sample = sc$orig.ident
sc$group[sc$sample %in% c("Control1", "Control2")] <- "Control"
sc$group[sc$sample %in% c("MC_1", "MC_2", "MC_3")] <- "MC"
sc$group <- factor(sc$group, levels = c("Control", "MC"))



markers <- "Ar"
sc = JoinLayers(sc)
p= FeatureStatPlot(sc, markers, reduction = "umap", plot_type = "dim",
                facet_by = "group", split_by = TRUE, ncol = 1,lower_quantile = 0,
                upper_quantile = 0.8)
save_plot(p, file.path(out_dir, "03_CD8_Subcluster/08_target_gene_show/figures/Ar_Featureplot"), 8, 6)


p = CellDimPlot(sc, group_by = "cell_type", reduction = "umap",label = T)
save_plot(p, file.path(out_dir, "03_CD8_Subcluster/08_target_gene_show/figures/UMAP_by_annotations"), 10, 10)

p = CellDimPlot(sc, group_by = "cell_type", reduction = "umap",label = T,split_by ="group")
save_plot(p, file.path(out_dir, "03_CD8_Subcluster/08_target_gene_show/figures/UMAP_by_annotations_splitgroup"), 20, 10)


cd8 = readRDS("output/rds/cd8_clean_final.rds")
allmarkers <- Seurat::FindAllMarkers(cd8,group.by ="cd8_cell_type")
write.csv(allmarkers,"output/03_CD8_Subcluster/08_target_gene_show/figures/cd8markers.csv")
table(allmarkers$cluster)

isg = allmarkers[which(allmarkers$cluster=="ISG high CD8"& allmarkers$gene %in% c('Ifit3',	'Isg15',	'Tnf'	,'Irf7')),]
teff = allmarkers[which(allmarkers$cluster=="Teff CD8"& allmarkers$gene %in% c('Ly6c2',	'Gzmk',	'Klrc1',	'Nkg7')),]
Tex = allmarkers[which(allmarkers$cluster=="Tex CD8"& allmarkers$gene %in% c('Tox',	'Pdcd1',	'Tigit',	'Serpinb6a')),]
Tpro= allmarkers[which(allmarkers$cluster=="Tpro CD8"& allmarkers$gene %in% c('Mki67',	'Cdk1',	'Ccnb2'	,'Birc5')),]
Treg = allmarkers[which(allmarkers$cluster=="Treg CD8"& allmarkers$gene %in% c('Foxp3'	,'Ctla4',	'Il2ra',	'Tnfrsf4')),]
Xcl1= allmarkers[which(allmarkers$cluster=="Xcl1 high Teff CD8"& allmarkers$gene %in% c('Xcl1',	'Cxcl10',	'Klf2',	'Fos')),]

modimarkers = rbind(isg,teff,Tex,Tpro,Treg,Xcl1)

# 显示图形
save_plot(pp, file.path(out_dir, "03_CD8_Subcluster/08_target_gene_show/figures/cd8markers"))


p_cd8_dot_top <- DotPlot(cd8, features = split (modimarkers$gene,modimarkers$cluster)) + RotatedAxis()
save_plot(
  p_cd8_dot_top,
  file.path(out_dir, "03_CD8_Subcluster/08_target_gene_show/figures/cd8markers")
)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)   # 小鼠
  library(dplyr)
  library(ggplot2)
})

# 只保留显著 & 上调（防止噪音）
allmarkers_filt <- allmarkers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.25)

gene_list_cd8 <- list(
  ISG_high_CD8 = allmarkers_filt %>%
    filter(cd8_cell_type == "ISG high CD8") %>% pull(gene) %>% unique(),
  
  Teff_CD8 = allmarkers_filt %>%
    filter(cd8_cell_type == "Teff CD8") %>% pull(gene) %>% unique(),
  
  Tex_CD8 = allmarkers_filt %>%
    filter(cd8_cell_type == "Tex CD8") %>% pull(gene) %>% unique(),
  
  Tpro_CD8 = allmarkers_filt %>%
    filter(cd8_cell_type == "Tpro CD8") %>% pull(gene) %>% unique(),
  
  Treg_CD8 = allmarkers_filt %>%
    filter(cd8_cell_type == "Treg CD8") %>% pull(gene) %>% unique(),
  
  Xcl1_high_Teff_CD8 = allmarkers_filt %>%
    filter(cd8_cell_type == "Xcl1 high Teff CD8") %>% pull(gene) %>% unique()
)

# 看一下每群基因数（很重要）
sapply(gene_list_cd8, length)
## ===============================
## GO enrichment (BP)
## ===============================
go_kegg_fig_dir <- file.path(out_dir, "03_CD8_Subcluster/09_GO_KEGG/figures")
dir.create(go_kegg_fig_dir, recursive = TRUE, showWarnings = FALSE)
go_res_list <- lapply(
  gene_list_cd8,
  function(genes) {
    enrichGO(
      gene          = genes,
      OrgDb         = org.Mm.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = TRUE
    )
  }
)

# 保存结果
for (nm in names(go_res_list)) {
  go_df <- as.data.frame(go_res_list[[nm]])
  if (nrow(go_df) > 0) {
    write.csv(
      go_df,
      file = file.path(
        out_dir,
        "03_CD8_Subcluster/09_GO_KEGG",
        paste0("GO_BP_", nm, ".csv")
      ),
      row.names = FALSE
    )
  }
}
## ===============================
## KEGG enrichment
## ===============================

kegg_res_list <- lapply(
  gene_list_cd8,
  function(genes) {
    
    # SYMBOL → ENTREZ
    eg <- bitr(
      genes,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = org.Mm.eg.db
    )
    
    if (is.null(eg) || nrow(eg) < 10) return(NULL)
    
    enrichKEGG(
      gene          = eg$ENTREZID,
      organism      = "mmu",
      pvalueCutoff  = 0.05,
      use_internal_data =T
    )
  }
)

# 保存结果
for (nm in names(kegg_res_list)) {
  if (is.null(kegg_res_list[[nm]])) next
  kegg_df <- as.data.frame(kegg_res_list[[nm]])
  if (nrow(kegg_df) > 0) {
    write.csv(
      kegg_df,
      file = file.path(
        out_dir,
        "03_CD8_Subcluster/09_GO_KEGG",
        paste0("KEGG_", nm, ".csv")
      ),
      row.names = FALSE
    )
  }
}

## ===============================
## Loop plotting: GO BP (save_plot)
## ===============================

for (nm in names(go_res_list)) {
  
  ego <- go_res_list[[nm]]
  ego <- simplify(
    ego,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message("Skip GO BP: ", nm)
    next
  }
  
  ## DotPlot
  p_dot <- dotplot(
    ego,
    showCategory = 15,
    title = paste0("GO BP – ", nm)
  ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  save_plot(
    p_dot,
    file.path(go_kegg_fig_dir, paste0("GO_BP_", nm, "_dotplot")),
    w = 8,
    h = 6
  )
  
  ## BarPlot
  p_bar <- barplot(
    ego,
    showCategory = 15,
    title = paste0("GO BP – ", nm)
  ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  save_plot(
    p_bar,
    file.path(go_kegg_fig_dir, paste0("GO_BP_", nm, "_barplot")),
    w = 8,
    h = 6
  )
}
## ===============================
## Loop plotting: KEGG (save_plot)
## ===============================

for (nm in names(kegg_res_list)) {
  
  ekk <- kegg_res_list[[nm]]
  if (is.null(ekk) || nrow(as.data.frame(ekk)) == 0) {
    message("Skip KEGG: ", nm)
    next
  }
  
  ## DotPlot
  p_dot <- dotplot(
    ekk,
    showCategory = 15,
    title = paste0("KEGG – ", nm)
  ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  save_plot(
    p_dot,
    file.path(go_kegg_fig_dir, paste0("KEGG_", nm, "_dotplot")),
    w = 8,
    h = 6
  )
  
  ## BarPlot
  p_bar <- barplot(
    ekk,
    showCategory = 15,
    title = paste0("KEGG – ", nm)
  ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  save_plot(
    p_bar,
    file.path(go_kegg_fig_dir, paste0("KEGG_", nm, "_barplot")),
    w = 8,
    h = 6
  )
}

p = FeatureStatPlot(cd8, features = c("Ar"),
                ident = "cd8_cell_type", plot_type = "box", facet_scales = "free_y")
save_plot(
  p,
  file.path(out_dir, "03_CD8_Subcluster/08_target_gene_show/figures/Ar_violin")
)
dev.off()
