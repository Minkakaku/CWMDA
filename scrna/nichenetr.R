rm(list=ls())
library(nichenetr)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(Seurat)
library(rliger)
library(scCustomize)
library(qs)

seuratObjV5 <- readRDS("Myeloid_fibroblast.rds")
seuratObj <- Convert_Assay(seurat_object = seuratObjV5, convert_to = "V3")

organism <- "human"

if(organism == "human"){  
  lr_network <- readRDS("lr_network_human_21122021.rds")  
  ligand_target_matrix <- readRDS("ligand_target_matrix_nsga2r_final.rds") 
  weighted_networks <- readRDS("weighted_networks_nsga2r_final.rds")
  } else if(organism == "mouse"){  
  lr_network <- readRDS("./mouse/lr_network_mouse_21122021.rds")  
  ligand_target_matrix <- readRDS("./mouse/ligand_target_matrix_nsga2r_final_mouse.rds")  
  weighted_networks <- readRDS("./mouse/weighted_networks_nsga2r_final_mouse.rds")
  }

lr_network <- lr_network %>% distinct(from, to)

table(seuratObj$cell_type)

receiver = "Fibro_L3"
Idents(seuratObj) = seuratObj$cell_type
expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.05)

all_receptors <- unique(lr_network$to)
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

sender_celltypes <- c("Myeloid Cells")
# 使用 lapply 获取每个 sender 细胞类型给定表达比例基因
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

condition_oi <-  "OA"
condition_reference <- "Nor"
seurat_obj_receiver <- subset(seuratObj, idents = receiver)
Idents(seurat_obj_receiver) = seurat_obj_receiver$Group
DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "Group", 
                                  min.pct = 0.05) %>% rownames_to_column("gene")
geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))

p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) +
  geom_histogram(color="black", fill="darkorange")  +   
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) +
  labs(x="ligand activity (OA)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity
ggsave("p_hist_lig_activity.png")
ggsave("p_hist_lig_activity.pdf")


best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
best_upstream_ligands
write.csv(best_upstream_ligands,"Top30 upstream_ligands.csv")

LigandActivity <- function(ligand_activities = ligand_activities, best_upstream_ligands = best_upstream_ligands){
  vis_ligand_aupr <- ligand_activities %>%     
    filter(test_ligand %in% best_upstream_ligands) %>%   
    column_to_rownames("test_ligand") %>%    
    select(aupr_corrected) %>%    
    arrange(aupr_corrected) %>%   
    as.matrix(ncol = 1)
  heatmapplot <- (make_heatmap_ggplot(vis_ligand_aupr, "Prioritized ligands", "Ligand activity", legend_title = "AUPR", color = "darkorange") + theme(axis.text.x.top = element_blank())) 
  vis_ligand_aupr_ncol <- ncol(vis_ligand_aupr)
  return(list(heatmapplot = heatmapplot, vis_ligand_aupr_ncol = vis_ligand_aupr_ncol))
}
LigandActivity(ligand_activities = ligand_activities, best_upstream_ligands = best_upstream_ligands)
ggsave("LigandActivity.png",width = 10)
ggsave("LigandActivity.pdf",width = 10)


active_ligand_target_links_df <- best_upstream_ligands %>%  lapply(get_weighted_ligand_target_links,
                                                                   geneset = geneset_oi,
                                                                   ligand_target_matrix = ligand_target_matrix,
                                                                   n = 100) %>%  bind_rows() %>% drop_na()
write.csv(active_ligand_target_links_df,"active_ligand_target_links_df.csv")
active_ligand_target_links <- prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df,
                                                                  ligand_target_matrix = ligand_target_matrix,
                                                                  cutoff = 0.33) 

RegulatoryPlot <- function(best_upstream_ligands = best_upstream_ligands, active_ligand_target_links = active_ligand_target_links){
  order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() 
  order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))
  vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])
  regheatmap <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes", color = "purple", legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple")
  vis_ligand_targetncol <- ncol(vis_ligand_target)
  return(list(regheatmap = regheatmap, vis_ligand_targetncol = vis_ligand_targetncol))
}
RegulatoryPlot(best_upstream_ligands = best_upstream_ligands, active_ligand_target_links = active_ligand_target_links)
ggsave("Predicted target genes.png",width = 30)
ggsave("Predicted target genes.pdf",width = 30)

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(  best_upstream_ligands, expressed_receptors,  lr_network, weighted_networks$lr_sig) 
write.csv(ligand_receptor_links_df,"ligand_receptor_links_df.csv")

receptorplot <- function(ligand_receptor_links_df = ligand_receptor_links_df, best_upstream_ligands = best_upstream_ligands){
  vis_ligand_receptor_network <- prepare_ligand_receptor_visualization( ligand_receptor_links_df, best_upstream_ligands, order_hclust = "both") 
  receptorp <- (make_heatmap_ggplot(t(vis_ligand_receptor_network), y_name = "Ligands", x_name = "Receptors",  color = "mediumvioletred", legend_title = "Prior interaction potential"))  
  return(receptorp)
}
receptorplot(ligand_receptor_links_df = ligand_receptor_links_df, best_upstream_ligands = best_upstream_ligands)
ggsave("Prior interaction potential.png",width = 15)
ggsave("Prior interaction potential.pdf",width = 15)

ligand_activities_all <- ligand_activities 
best_upstream_ligands_all <- best_upstream_ligands
ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>%  pull(test_ligand) %>% unique()
p_ligand_auprs <- LigandActivity(ligand_activities = ligand_activities, best_upstream_ligands = best_upstream_ligands)
p_ligand_aupr <- p_ligand_auprs$heatmapplot
p_ligand_aupr
ggsave("Sender-focused-LigandActivity.png",width = 10)
ggsave("Sender-focused-LigandActivity.pdf",width = 10)

active_ligand_target_links_df <- best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,     
                                                                  geneset = geneset_oi,    
                                                                  ligand_target_matrix = ligand_target_matrix,    
                                                                  n = 100) %>% bind_rows() %>% drop_na()
active_ligand_target_links <- prepare_ligand_target_visualization( ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33) 
p_ligand_targets <- RegulatoryPlot(best_upstream_ligands = best_upstream_ligands, active_ligand_target_links = active_ligand_target_links)
p_ligand_target <- p_ligand_targets$regheatmap
p_ligand_target
ggsave("Sender-focused-Predicted target genes.png",width = 30)
ggsave("Sender-focused-Predicted target genes.pdf",width = 30)

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(  best_upstream_ligands, expressed_receptors,  lr_network, weighted_networks$lr_sig) 
p_ligand_receptor <- receptorplot(ligand_receptor_links_df = ligand_receptor_links_df, best_upstream_ligands = best_upstream_ligands)
p_ligand_receptor
ggsave("Sender-focused-Prior interaction potential.png",width = 15)
ggsave("Sender-focused-Prior interaction potential.pdf",width = 15)


# Dotplot of sender-focused 模式
p_dotplot <- DotPlot(subset(seuratObj, cell_type %in% sender_celltypes), features = rev(best_upstream_ligands), cols = "RdYlBu") + coord_flip() + scale_y_discrete(position = "right")
p_dotplot
ggsave("Dotplot of sender-focusedl.png")
ggsave("Dotplot of sender-focused.pdf")

celltype_order <- levels(Idents(seuratObj)) 
DE_table_top_ligands <- lapply(  celltype_order[celltype_order %in% sender_celltypes],  get_lfc_celltype,   seurat_obj = seuratObj,  condition_colname = "Group",  condition_oi = condition_oi,  condition_reference = condition_reference,  celltype_col = "cell_type",  min.pct = 0, logfc.threshold = 0,  features = best_upstream_ligands ) 
DE_table_top_ligands <- DE_table_top_ligands %>%  reduce(., full_join) %>%   column_to_rownames("gene") 
vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), , drop = FALSE])
p_lfc <- make_threecolor_heatmap_ggplot(vis_ligand_lfc,                                "Prioritized ligands", "LFC in Sender",                                low_color = "midnightblue", mid_color = "white",                                mid = median(vis_ligand_lfc), high_color = "red",                                legend_title = "LFC")
p_lfc
ggsave("Case - Control sender-focusedl.png")
ggsave("Case - Control sender-focused.pdf")

