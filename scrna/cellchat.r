rm(list=ls())
# contamination sample process ----
## load lib ----
# function 
source("E:\\SOP\\src\\packages_loaded.R")
packages_loaded(c("Seurat",
                  "SeuratWrappers",
                  "ggplot2",
                  "dplyr",
                  "scDblFinder",
                  "patchwork",
                  "dplyr",
                  "clusterProfiler",
                  "org.Hs.eg.db",
                  "decontX",
                  "RColorBrewer",
                  "reshape2",
                  "patchwork",
                  "scales",
                  "BiocParallel"))

major  = readRDS("seu_final_annotated.rds")
fb = readRDS("scRNA_fib_annotated.rds")

m_seurat <- subset(major, cell_type == "Myeloid Cells")
m_seurat$cell_type

fb$cell_type <- paste0("Fibro_", fb$cell_type_sub)
combined_seurat <- merge(m_seurat, y = fb)
combined_seurat <- NormalizeData(combined_seurat)
combined_seurat <- ScaleData(combined_seurat)


# factor
table(combined_seurat$Group)
levels(combined_seurat$Group)
combined_seurat$Group <- factor(
  combined_seurat$Group,
  levels = c(
    "Nor" ,
    "OA" 
  )
)
table(combined_seurat$Group)

rm(fb)
rm(combined_seurat)
rm(m_seurat)
rm(macrophage_seurat)
rm(major)
# cellchat
library(CellChat)

Idents(combined_seurat) <- combined_seurat$Group
Nor <- subset(combined_seurat, idents = "Nor")  # 对照组
Nor = JoinLayers(Nor)
OA <- subset(combined_seurat, idents = "OA")    # 脑转移组
OA = JoinLayers(OA)
table(Nor$cell_type)
table(OA$cell_type)
## Nor
data.input <- GetAssayData(Nor, layer = 'data')
meta <- Nor@meta.data[, c("Group", "cell_type")]
colnames(meta) <- c("Group", "CellTypes")
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "CellTypes")
CellChatDB <- CellChatDB.human  # 使用人类细胞通讯数据库
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)  # 数据预处理
cellchat <- identifyOverExpressedGenes(cellchat)  # 筛选高表达基因
cellchat <- identifyOverExpressedInteractions(cellchat)  # 筛选高表达互作
cellchat <- smoothData(cellchat, adj =PPI.human)  # 投影到人类蛋白质互作网络
cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE)  # 计算通讯概率
cellchat <- filterCommunication(cellchat, min.cells = 10)  # 过滤低频互作
cellchat <- computeCommunProbPathway(cellchat)  # 计算通路层面的通讯概率
cellchat <- aggregateNet(cellchat)  # 聚合网络
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")  # 计算网络中心性
cellchat <- computeNetSimilarity(cellchat, type = "functional")  # 计算功能相似性
cco.Nor <- cellchat
saveRDS(cco.Nor,"mfbcco_Nor.rds")

## OA
data.input <- GetAssayData(OA, layer = 'data')
meta <- OA@meta.data[, c("Group", "cell_type")]
colnames(meta) <- c("Group", "CellTypes")
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "CellTypes")
CellChatDB <- CellChatDB.human  # 使用人类细胞通讯数据库
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)  # 数据预处理
cellchat <- identifyOverExpressedGenes(cellchat)  # 筛选高表达基因
cellchat <- identifyOverExpressedInteractions(cellchat)  # 筛选高表达互作
cellchat <- smoothData(cellchat, adj =PPI.human)  # 投影到人类蛋白质互作网络
cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE)  # 计算通讯概率
cellchat <- filterCommunication(cellchat, min.cells = 10)  # 过滤低频互作
cellchat <- computeCommunProbPathway(cellchat)  # 计算通路层面的通讯概率
cellchat <- aggregateNet(cellchat)  # 聚合网络
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")  # 计算网络中心性
cellchat <- computeNetSimilarity(cellchat, type = "functional")  # 计算功能相似性
cco.OA <- cellchat
saveRDS(cco.OA,"mfbcco_OA.rds")

cco.OA = readRDS("mfbcco_OA.rds")
cco.Nor = readRDS("mfbcco_Nor.rds")
cco.list <- list(Nor = cco.Nor, OA = cco.OA)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)

gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "weight")
p <- gg1 + gg2
ggsave("比较细胞通讯相互作用的总数和相互作用的强度.png", p, width = 6, height = 4)
ggsave("比较细胞通讯相互作用的总数和相互作用的强度.pdf", p, width = 6, height = 4)


pdf("不同细胞群间相互作用次数或相互作用强度的差异.pdf", width = 7, height = 7)
par(mfrow = c(1,2), xpd=TRUE)
p = netVisual_diffInteraction(cellchat, weight.scale = T,comparison = c(1,2))
p = netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",comparison = c(1,2))
p
dev.off()


pdf("不同细胞群间相互作用次数或相互作用强度的差异(热图).pdf", width =8, height =4)
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()

pdf("singlesample_interaction_strength.pdf", width =15, height =15)
weight.max <- getMaxWeight(cco.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_circle(cco.list[[i]]@net$count, weight.scale = T, arrow.size = 1,label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(cco.list)[i]))
}
dev.off()

num.link <- sapply(cco.list, function(x) {rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link))
gg <- list()
for (i in 1:length(cco.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(cco.list[[i]], title = names(cco.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)
ggsave("Number_of_interaction_strength.pdf",width = 15,height = 7)

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional", do.parallel = FALSE)
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
ggsave("Umap_of_interaction_strength.pdf",width = 15,height = 15)


rankSimilarity(cellchat, type = "functional",comparison1 =c(1,2),comparison2 =c(1,2))


gg1 <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE)
gg1 + gg2
ggsave("Histon_of_interaction_strength.pdf",width = 10,height = 10)


library(ComplexHeatmap)
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)

ht1 <- netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "all", signaling = pathway.union, 
                                         title = names(cco.list)[1], width = 12, height = 28)
ht2 <- netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "all", signaling = pathway.union, 
                                         title = names(cco.list)[2], width = 12, height = 28)

pdf("Cellchat_pathway_compared.pdf", width = 25, height = 15, pointsize = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()


library(ComplexHeatmap)
i = 1
pathway.union <- union(cco.list [[i]]@netP$pathways, cco.list [[i+1]]@netP$pathways)
##outgoing
ht1 = netAnalysis_signalingRole_heatmap(cco.list [[i]], pattern = "outgoing", signaling = pathway.union, title = names(cco.list )[i], width = 10, height =20)
ht2 = netAnalysis_signalingRole_heatmap(cco.list [[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(cco.list )[i+1], width = 10, height = 20)
pdf("比较与每个细胞群相关的传出outgoing信号.pdf", width =15, height =15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()



##incoming
ht1 = netAnalysis_signalingRole_heatmap( cco.list[[i]], pattern = "incoming", 
                                        signaling = pathway.union, title = names( cco.list)[i],
                                        width = 10, height = 20, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap( cco.list[[i+1]], pattern = "incoming",
                                        signaling = pathway.union, title = names( cco.list)[i+1],
                                        width = 10, height = 20, color.heatmap = "GnBu")
#ht3 = netAnalysis_signalingRole_heatmap( cco.list[[i+2]], pattern = "incoming", 
#                                        signaling = pathway.union, title = names( cco.list)[i+2],
#                                        width = 5, height = 6, color.heatmap = "GnBu")
#draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm"))
pdf("比较与每个细胞群相关的传入incoming信号.pdf", width =15, height =15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()



##all
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[i]], pattern = "all", 
                                        signaling = pathway.union, title = names(cco.list)[i], 
                                        width = 10, height = 20, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[i+1]], pattern = "all",
                                        signaling = pathway.union, title = names(cco.list)[i+1], 
                                        width = 10, height = 20, color.heatmap = "OrRd")
#ht3 = netAnalysis_signalingRole_heatmap(cco.list[[i+2]], pattern = "all", 
#                                        signaling = pathway.union, title = names(cco.list)[i+2], 
#                                        width = 5, height = 6, color.heatmap = "OrRd")
#draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm"))
pdf("比较与每个细胞群相关的所有信号.pdf", width =15, height =15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

# 通过比较通讯概率来识别故障信号
pdf("通过比较通讯概率来识别故障信号(所有细胞).pdf", width =15, height =60)
netVisual_bubble(cellchat,  comparison = c(1, 2), angle.x = 45)
dev.off()

pdf("通过比较通讯概率来识别故障信号(L-myleoid).pdf", width =15, height =15)
netVisual_bubble(cellchat, sources.use = c(1:3), targets.use =7,  comparison = c(1, 2),angle.x = 45)
dev.off()
pdf("通过比较通讯概率来识别故障信号(myleoid-L).pdf", width =15, height =15)
netVisual_bubble(cellchat, sources.use = c(1:3), targets.use =7,  comparison = c(1, 2),angle.x = 45)
dev.off()


print(paste(1:length(levels( cellchat@idents$joint )),levels( cellchat@idents$joint ),sep=":")) #合并的

pdf("通过比较通讯概率来识别故障信号(L3-myleoid).pdf", width =15, height =15)
gg1 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = 7,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in OA", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = 7,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in OA", angle.x = 45, remove.isolate = T)
gg1 + gg2
dev.off()

pdf("通过比较通讯概率来识别故障信号(myleoid-L3).pdf", width =15, height =15)
gg1 <- netVisual_bubble(cellchat, sources.use = 7, targets.use = 3,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in OA", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 7, targets.use = 3,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in OA", angle.x = 45, remove.isolate = T)
gg1 + gg2
dev.off()

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "OA"

# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset

# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "datasets", 
                                       pos.dataset = pos.dataset,
                                       features.name = features.name,
                                       only.pos = FALSE,
                                       thresh.pc = 0.1,
                                       thresh.fc = 0.1,
                                       thresh.p = 1)

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)

# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net,
                              datasets = "OA",ligand.logFC = 0.2,
                              receptor.logFC = NULL)

# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net,
                                datasets = "Nor",
                                ligand.logFC = -0.1,
                                receptor.logFC = -0.1)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[,"interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat,
                        pairLR.use = pairLR.use.up,
                        # sources.use = 4, targets.use = c(5:11),
                        comparison = c(1, 2),
                        angle.x = 90,
                        remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(cco.list)[2]))

pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat,
                        pairLR.use = pairLR.use.down,
                        # sources.use = 4, targets.use = c(5:11),
                        comparison = c(1, 2),
                        angle.x = 90, remove.isolate = T,
                title.name = paste0("Down-regulated signaling in ", names(cco.list)[2]))
pdf("通过差异分析来识别故障信号(所有).pdf", width =30, height =20)
gg1 + gg2
dev.off()

gg1 <- netVisual_bubble(cellchat,
                        pairLR.use = pairLR.use.up,
                        sources.use = 3, targets.use = 7,
                        comparison = c(1, 2),
                        angle.x = 90,
                        remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(cco.list)[2]))

pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat,
                        pairLR.use = pairLR.use.down,
                        sources.use = 3, targets.use = 7,
                        comparison = c(1, 2),
                        angle.x = 90, remove.isolate = T,
                        title.name = paste0("Down-regulated signaling in ", names(cco.list)[2]))
pdf("通过差异分析来识别故障信号(L3-myleoid).pdf", width =7, height =10)
gg1 + gg2
dev.off()

gg1 <- netVisual_bubble(cellchat,
                        pairLR.use = pairLR.use.up,
                        sources.use = 7, targets.use = 3,
                        comparison = c(1, 2),
                        angle.x = 90,
                        remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(cco.list)[2]))

pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat,
                        pairLR.use = pairLR.use.down,
                        sources.use = 7, targets.use = 3,
                        comparison = c(1, 2),
                        angle.x = 90, remove.isolate = T,
                        title.name = paste0("Down-regulated signaling in ", names(cco.list)[2]))
pdf("通过差异分析来识别故障信号(myleoid-L3).pdf", width =7, height =10)
gg1 + gg2
dev.off()

# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
pdf("使用和弦图可视化上调和下调的信号配体受体对(myleoid-L3).pdf", width =10, height =10)
p = netVisual_chord_gene(cco.list[[2]],
                     sources.use = 7,targets.use = 3,
                     slot.name = 'net',net = net.up, 
                     lab.cex = 0.8,small.gap = 3.5,
                     title.name = paste0("Up-regulated signaling in ",names(cco.list)[2]))
p = netVisual_chord_gene(cco.list[[1]],
                     sources.use = 7, targets.use =3,
                     slot.name = 'net', net = net.down,
                     lab.cex = 0.8, small.gap = 3.5,
                     title.name = paste0("Down-regulated signaling in ", names(cco.list)[2]))
p
dev.off()

par(mfrow = c(1,2), xpd=TRUE)
pdf("使用和弦图可视化上调和下调的信号配体受体对(L3-myleoid).pdf", width =10, height =10)
p = netVisual_chord_gene(cco.list[[2]],
                         sources.use = 3,targets.use = 7,
                         slot.name = 'net',net = net.up, 
                         lab.cex = 0.8,small.gap = 3.5,
                         title.name = paste0("Up-regulated signaling in ",names(cco.list)[2]))
p = netVisual_chord_gene(cco.list[[1]],
                         sources.use = 3, targets.use =7,
                         slot.name = 'net', net = net.down,
                         lab.cex = 0.8, small.gap = 3.5,
                         title.name = paste0("Down-regulated signaling in ", names(cco.list)[2]))
p
dev.off()



pdf("上调和下调的信号配体受体对词云(上调).pdf", width =6, height =6)
# visualize the enriched ligands in the second condition
computeEnrichmentScore(net.up, species = 'human')
dev.off()
pdf("上调和下调的信号配体受体对词云(下调).pdf", width =6, height =6)
# visualize the enriched ligands in the first condition
computeEnrichmentScore(net.down, species = 'human')
dev.off()