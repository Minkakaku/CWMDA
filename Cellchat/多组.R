rm(list = ls())
library(cowplot)
library(Seurat)
library(RColorBrewer)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(patchwork)
library(getopt)
library(CellChat)
options(stringsAsFactors = FALSE)
outdir = "/lvdata/wzb/scRNA/FWSC20240315/2024.09.19/Cellchat"
inrds = "/lvdata/wzb/scRNA/FWSC20240315/2024.09.19/Cellchat/h5_trans.rds"
species = 'mouse'
# pos.dataset = "iMCD"
# ref.dataset = "NC"



setwd(outdir)


compare.list = factor(c(pos.dataset,ref.dataset),levels = c(ref.dataset,pos.dataset))
Seuset = readRDS(inrds)
cellchat.list = SplitObject(Seuset,split.by = "group")
meta.list = lapply(cellchat.list,FUN = function(x){
  return(x@meta.data)
})
group = as.character(unique(Seuset$group))
cellchat.objet.list = lapply(group,FUN = function(x){
  createCellChat(cellchat.list[[x]],meta = meta.list[[x]],group.by = "cell_type")
})
names(cellchat.objet.list) = group

runCellchat = function(cellchat,species = "human",use.ident = "cell_type",filter.cells = 10,name_ct){
  print(name_ct)
  cellchat = setIdent(cellchat, ident.use = use.ident) # 将 "labels" 设为默认细胞标记类型，这个可以根据自己的数据自定义
  if (species=="human") {
    CellChatDB = CellChatDB.human  
  }else{
    CellChatDB = CellChatDB.mouse  
  }
  write.table(CellChatDB$interaction,'cellchatDB.xls',sep = "\t")
  showDatabaseCategory(CellChatDB)
  # 我们这里使用“Secreted Signaling”部分做后续的细胞通讯分析
  cellchat@DB = subsetDB(CellChatDB,search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact")) 
  cellchat = subsetData(cellchat)
  # future::plan("multiprocess", workers = 4)
  cellchat = identifyOverExpressedGenes(cellchat)###首先识别过表达基因（配体——受体）
  cellchat = identifyOverExpressedInteractions(cellchat)###然后识别过表达配受体之间过表达的相互作用   （绕绕绕绕绕....）
  # project gene expression data onto PPI network (optional)
  if(species=="human"){
    cellchat = projectData(cellchat, PPI.human)
  }else{
    cellchat = projectData(cellchat, PPI.mouse)
  }
  # group.new = levels(cellchat@idents)
  # cell.chat = liftCellChat(cellchat,group.new)
  cellchat@idents = droplevels(cellchat@idents, exclude = setdiff(levels(cellchat@idents),unique(cellchat@idents)))
  cellchat = computeCommunProb(cellchat,type="triMean",trim = 0)###这应该是核心的一行代码
  cellchat = filterCommunication(cellchat, min.cells = filter.cells)##过滤掉小于10个细胞的胞间通讯网络
  cellchat = computeCommunProbPathway(cellchat)
  df.net = subsetCommunication(cellchat)
  relation = paste0(df.net[,1],"->",df.net[,2])
  df.net = cbind(relation,df.net)
  colnames(df.net)[1]="relation"
  write.table(df.net,paste0(name_ct,".","net_lr.txt"),sep="\t",quote=F,row.names=F,col.names=T)
  #保存信号通路细胞通讯网络
  cellchat = computeCommunProbPathway(cellchat)
  df.netp = subsetCommunication(cellchat,slot.name="netP")
  write.table(df.netp,paste0(name_ct,".","netPathways_lr.txt"),sep="\t",quote=F,row.names=F,col.names=T)
  cellchat = aggregateNet(cellchat)
  cellchat = netAnalysis_computeCentrality(cellchat)
  return(cellchat)
}

cellchat.objet.list.new = lapply(group,FUN = function(x){
  runCellchat(cellchat.objet.list[[x]],species = species,name_ct=x)
})

names(cellchat.objet.list.new) = group
idents.list = lapply(cellchat.objet.list.new,FUN = function(x){
  length(levels(x@idents))
})


# 这里判断是否包含的细胞类型数目一致不一致的话选择里面数量最多的，其他的看齐
if (length(unique(unlist(idents.list)))>1) {
  maxgroup = names(idents.list)[idents.list == max(unique(unlist(idents.list)))]
  group.max = levels(cellchat.objet.list.new[[maxgroup]]@idents)
  cellchat.objet.list.new = lapply(cellchat.objet.list.new, function(x){
    liftCellChat(x, group.max)
  })

}

cellchat = mergeCellChat(cellchat.objet.list.new,add.names = names(cellchat.objet.list.new))
celltype = unique(cellchat@meta$cell_type) %>% as.character()
save(cellchat,file = "cellchat_object.RData")
save(cellchat.objet.list.new,file = "cellchat_object_list.RData")


plot_save = function(filename,plt,dpi=300,w=7,h=7){
  ggsave(paste0(filename,'.pdf'),plt,height = h,width = w)
  
}
fig_save = function(filename,plt,dpi=300,w=7,h=7){
  pdf(paste0(filename,'.pdf'),height = h,width = w)
  draw(plt)
  dev.off()
}
recor_save = function(filename,plt,dpi=300,w=7,h=7){
  pdf(paste0(filename,'.pdf'),height = h,width = w)
  replayPlot(plt)
  dev.off()
}
compare.list = c(2,1)
gg1 = compareInteractions(cellchat, show.legend = F,measure = "count")
gg1
plot_save(filename = "Number of inferred interactions.barplot",plt = gg1)
gg2 = compareInteractions(cellchat, show.legend = F, measure = "weight")
plot_save("Interaction strength.barplot",gg2)
library(ComplexHeatmap)
for (i in 1:length(group)) {
  for (j in (1:length(group))[(1:length(group)) !=i]) {
    compare.list = c(i,j)
    
    # cellchat@netP$
    gg3 = netVisual_diffInteraction(cellchat,  comparison =compare.list,
                                    weight.scale = T,measure = "count")
    dev.off()
    recor_save(paste0(group[i],"VS",group[j],"_Differential number of interactions.net"),gg3)
    gg4 = netVisual_diffInteraction(cellchat, comparison = compare.list,
                                    weight.scale = T, measure = "weight")
    recor_save(paste0(group[i],"VS",group[j],"_Differential interaction strength.net"),gg4)
    gg5 = netVisual_heatmap(cellchat, comparison = compare.list)
    fig_save(paste0(group[i],"VS",group[j],"_Differential number of interactions.heatmap"),gg5)
    gg6 = netVisual_heatmap(cellchat, measure = "weight", comparison = compare.list)
    fig_save(paste0(group[i],"VS",group[j],"_Differential interaction strength.heatmap"),gg6)
    
    gg7 = rankNet(cellchat,
                  mode = "comparison",
                  measure = "weight",
                  sources.use = NULL,
                  targets.use = NULL, comparison =compare.list ,
                  stacked = T,
                  do.stat = TRUE)
    plot_save(paste0(group[i],"VS",group[j],"_Relative information flow.weight"),gg7)
    gg8 = rankNet(cellchat,
                  mode = "comparison",
                  measure = "weight",
                  sources.use = NULL,
                  targets.use = NULL, comparison = compare.list,
                  stacked = F,
                  do.stat = TRUE)
    plot_save(paste0(group[i],"VS",group[j],"_Information flow.weight"),gg8)
    
    gg9 = rankNet(cellchat,
                  mode = "comparison",
                  measure = "count",
                  sources.use = NULL,
                  targets.use = NULL,
                  comparison =compare.list,
                  stacked = T,
                  do.stat = TRUE)
    plot_save(paste0(group[i],"VS",group[j],"_Relative information flow.count"),gg9)
    gg10 = rankNet(cellchat,
                   mode = "comparison",
                   measure = "count",
                   sources.use = NULL,
                   targets.use = NULL, comparison = compare.list,
                   stacked = F,
                   do.stat = TRUE)
    plot_save(paste0(group[i],"VS",group[j],"_Information flow.count"),gg10)
    
    
    
    dif_pathway = as.character(gg8$data[c(1:3,(nrow(gg8$data)-3):nrow(gg8$data)),1])
    
    # 比较群体间的传入传出区别
    library(ComplexHeatmap)
    pathway.list = lapply(cellchat@netP, function(x){
      unique(x$pathways)
    }) %>% unlist() %>% unique()
    # group.max = 
    ht1 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[i]], pattern = "outgoing", 
                                            signaling = pathway.list,
                                            title = names(cellchat.objet.list.new)[i], 
                                            width = 10, height = length(pathway.list)/3)
    ht2 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[j]], pattern = "outgoing",
                                            signaling = pathway.list,
                                            title = names(cellchat.objet.list.new)[j],
                                            width = 10, height = length(pathway.list)/3)
    fig_save(paste0(group[i],"VS",group[j],'_Outgoing compare.heatmap'),plt = draw(ht1 + ht2, ht_gap = unit(0.5, "cm")),
             dpi = 300,w=10,h=length(pathway.list)/3)
    
    ht3 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[i]], pattern = "incoming", 
                                            signaling = pathway.list,
                                            title = names(cellchat.objet.list.new)[i], 
                                            width = 10, height = length(pathway.list)/3)
    ht4 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[j]], pattern = "incoming",
                                            signaling = pathway.list,
                                            title = names(cellchat.objet.list.new)[j], 
                                            width = 10, height = length(pathway.list)/3)
    fig_save(paste0(group[i],"VS",group[j],'_Incoming compare.heatmap'),plt = draw(ht3 + ht4, ht_gap = unit(0.5, "cm")),
             dpi = 300,w=10,h=length(pathway.list)/3)
    
    ht5 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[i]], pattern = "all", 
                                            signaling = pathway.list,
                                            title = names(cellchat.objet.list.new)[i], 
                                            width = 10, height = length(pathway.list)/3)
    ht6 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[j]], pattern = "all",
                                            signaling = pathway.list,
                                            title = names(cellchat.objet.list.new)[j], 
                                            width = 10, height = length(pathway.list)/3)
    fig_save(paste0(group[i],"VS",group[j],'_All compare.heatmap'),plt = draw(ht5 + ht6, ht_gap = unit(0.5, "cm")),
             dpi = 300,w=10,h=length(pathway.list)/3)
    
    
  }
}

# 比较群体间的传入传出区别
library(ComplexHeatmap)
pathway.list = lapply(cellchat@netP, function(x){
  unique(x$pathways)
}) %>% unlist() %>% unique()

ht1 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[1]], pattern = "outgoing", 
                                        signaling = pathway.list,
                                        title = names(cellchat.objet.list.new)[1], 
                                        width = length(group.max)/3, height = length(pathway.list)/3)
ht2 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[2]], pattern = "outgoing",
                                        signaling = pathway.list,
                                        title = names(cellchat.objet.list.new)[2],
                                        width = length(group.max)/3, height = length(pathway.list)/3)
fig_save('Outgoing compare.heatmap',plt = draw(ht1 + ht2, ht_gap = unit(0.5, "cm")),
         dpi = 300,w=length(group.max)/3,h=length(pathway.list)/3)

ht3 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[1]], pattern = "incoming", 
                                        signaling = pathway.list,
                                        title = names(cellchat.objet.list.new)[1], 
                                        width = length(group.max)/3, height = length(pathway.list)/3)
ht4 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[2]], pattern = "incoming",
                                        signaling = pathway.list,
                                        title = names(cellchat.objet.list.new)[2], 
                                        width = length(group.max)/3, height = length(pathway.list)/3)
fig_save('Incoming compare.heatmap',plt = draw(ht3 + ht4, ht_gap = unit(0.5, "cm")),
         dpi = 300,w=length(group.max)/3,h=length(pathway.list)/3)

ht5 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[1]], pattern = "all", 
                                        signaling = pathway.list,
                                        title = names(cellchat.objet.list.new)[1], 
                                        width = length(group.max)/3, height = length(pathway.list)/3)
ht6 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[2]], pattern = "all",
                                        signaling = pathway.list,
                                        title = names(cellchat.objet.list.new)[2], 
                                        width = length(group.max)/3, height = length(pathway.list)/3)
fig_save('All compare.heatmap',plt = draw(ht5 + ht6, ht_gap = unit(0.5, "cm")),
         dpi = 300,w=length(group.max)/3,h=length(pathway.list)/3)
dir.name = dir(pattern = "net_lr")
out.tmp = NULL
for (file in dir.name) {
  tmp = read.csv(file,sep = "\t")
  tmp$cluster = unlist(strsplit(file,".net_lr"))[1]
  out.tmp = rbind(out.tmp,tmp)
}
write.table(out.tmp,'All.net.lr.xls',sep = "\t",quote = F,row.names = F)


# 暂时先不跑
gg8$data$name
pathway.use = c('CCL','IFN-II')
# pathway.use = NULL
dir.create("compare_bubble")
for (j in pathway.use) {
  

for (i in celltype) {
  # i = 'M2_Mac_3'
  if (is.null(pathway.use)) {
    pathway.select = dif_pathway
  }else{
    # pathway.select = pathway.use
    pathway.select = j
  }
  try({
  gg9 = netVisual_bubble(cellchat,  
                   sources.use = i,
                   remove.isolate = T,
                   comparison = compare.list,
                   signaling = pathway.select,
                   title.name = paste0('singnal communication'),
                   angle.x = 90)
  gg9
  plot_save(paste0("compare_bubble/","Path_",j,"_celltype_",i," singnal communication"),plt = gg9,w=10,h=10)
  })
}
}

load('/lvdata/wzb/scRNA/FW2024-213/2024.04.23/Cellchat_sub/cellchat_object.RData')
i = 'NKT'
pairLR.use = extractEnrichedLR(cellchat, signaling = c("IFN-II","CCL"),enriched.only	 =T)
gg9 = netVisual_bubble(cellchat,  
                       sources.use = c("NKT",
                                       "KLRC2 NK",
                                       "GZMK NK",
                                       "CX3CR1 NK",
                                       "CD8TEX"
                       ),
                       targets.use = c("PADI Monocytes",
                                       "Monocytes/Macrophages",
                                       "IL7R Monocytes",
                                       "CCL Monocytes",
                                       "IFI Monocytes",
                                       "cDC",
                                       "pDC"
                       ),
                       remove.isolate = T,
                       comparison = compare.list,
                       pairLR.use = pairLR.use,
                       title.name = paste0('singnal communication'),
                       angle.x = 90)
ggsave(filename = "IFNG_CCL_bubble.pdf",gg9)
ggsave(filename = "IFNG_CCL_bubble.png",gg9)
cellchat@netP$iMCD$pathways
plot_save(paste0("compare_bubble/","Path_",j,"_celltype_",i," singnal communication"),plt = gg9,w=10,h=10)

rm(list = ls())
load('cellchat_object.RData')
cellchat@meta$datasets = factor(cellchat@meta$datasets)
extractEnrichedLR(cellchat, signaling = "BMP")
dev.off()
plotGeneExpression(cellchat, signaling = "BMP", split.by = "datasets", colors.ggplot = T, type = "violin")
# 特异上调的下调的配受体关系对

features.name = pos.dataset
cellchat_dif = identifyOverExpressedGenes(cellchat,
                                       group.dataset = "datasets",
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name,
                                       only.pos = FALSE,
                                       thresh.pc = 0.1, 
                                       thresh.fc = 0.1, 
                                       thresh.p = 0.05)
net = netMappingDEG(cellchat_dif, features.name = features.name)
net.up = subsetCommunication(cellchat_dif,
                              net = net,
                              datasets = pos.dataset,
                              # sources.use = 3,
                              ligand.logFC = 0.2,
                              receptor.logFC = NULL)
net.down = subsetCommunication(cellchat_dif,
                                net = net, 
                                datasets = ref.dataset,
                                # sources.use = 3,
                                ligand.logFC = -0.1,
                                receptor.logFC = -0.1)
gene.up = extractGeneSubsetFromPair(net.up, cellchat_dif)
gene.down = extractGeneSubsetFromPair(net.down, cellchat_dif)
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat_dif, pairLR.use = pairLR.use.up,
                        # sources.use = 3, targets.use = c(6,11),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = F,
                        title.name = paste0("Up-regulated signaling in ", pos.dataset))
gg1
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down,
                        sources.use = 3, targets.use = c(6,11),
                        comparison = c(1, 2), 
                        angle.x = 90, 
                        remove.isolate = T,
                        title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
gg1 + gg2

load('../cellchat_object_list.RData')
load('../cellchat_object.RData')



# load('./cellchat_object_list.RData')
# load('./cellchat_object.RData')

levels(cellchat@meta$cell_type)
cellchat_imcd
vertex.receiver = 25
pathways.show<- c("IFN-II",'CCL','BAFF')
# pathways.show<- c(
#                   "CXCL",
#                   "IFN−II",
#                   "CD40",
#                   "CD6",
#                   "GP1BA",
#                   "BAG",
#                   "SN",
#                   "CCL",
#                   "IL1",
#                   "VISFATIN",
#                   "PARS",
#                   "ANNEXIN",
#                   "ADGRE5",
#                   "ALCAM"
# )
pathways.show<- c('BAFF')
cellchat_iMCD = cellchat.objet.list.new$iMCD
cellchat_NC = cellchat.objet.list.new$NC
setwd('/lvdata/wzb/scRNA/FW2024-225_03/cellchat/BAFF')
netAnalysis_contribution(cellchat_iMCD, signaling = pathways.show)
ggsave('iMCD_barplot_BAFF.pdf')

netAnalysis_contribution(cellchat_NC, signaling = pathways.show)
ggsave('NC_barplot_BAFF.pdf')


cellchat_iMCD <- netAnalysis_computeCentrality(cellchat_iMCD, slot.name = "netP")

pdf('iMCD_BAFF.pdf')
netAnalysis_signalingRole_network(cellchat_iMCD, signaling = pathways.show, width = 14, height = 2.5, font.size = 10)
dev.off()

cellchat_NC <- netAnalysis_computeCentrality(cellchat_NC, slot.name = "netP")
pdf('NC_BAFF.pdf')
netAnalysis_signalingRole_network(cellchat_NC, signaling = pathways.show, width = 14, height = 2.5, font.size = 10)
dev.off()

pdf('scatter_iMCD_BAFF.pdf')
netAnalysis_signalingRole_scatter(cellchat_iMCD, signaling = pathways.show)
dev.off()
pdf('scatter_NC_BAFF.pdf')
netAnalysis_signalingRole_scatter(cellchat_NC, signaling = pathways.show)
dev.off()

gg2 = plotGeneExpression(cellchat_iMCD, signaling = pathways.show, enriched.only = TRUE, type = "violin")
ggsave('iMCD_Violin.pdf')
gg3 =plotGeneExpression(cellchat_NC, signaling = pathways.show, enriched.only = TRUE, type = "violin")
ggsave('NC_Violin.pdf')
num.link <- sapply(cellchat.objet.list.new, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link))
gg <- list()
for (i in 1:length(cellchat.objet.list.new)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(cellchat.objet.list.new[[i]], title = names(cellchat.objet.list.new)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)



setwd('/lvdata/wzb/scRNA/FW2024-225_03/cellchat/IFNG')
vertex.receiver = seq(1,15)
gg1 = netVisual_aggregate(cellchat_iMCD, signaling = pathways.show,  
                    vertex.receiver = vertex.receiver,
                    pt.title= 10,
                    layout = "hierarchy")
recor_save('IFNG_iMCD_Hierarchy',gg1,w = 22,h = 8)
netVisual_aggregate(cellchat_NC, signaling = pathways.show,  
                    vertex.receiver = vertex.receiver,
                    layout = "hierarchy")


netAnalysis_contribution(cellchat_iMCD, signaling = pathways.show)
netAnalysis_contribution(cellchat_NC, signaling = pathways.show)


cellchat_iMCD <- netAnalysis_computeCentrality(cellchat_iMCD, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat_iMCD, signaling = pathways.show, width = 14, height = 2.5, font.size = 10)
cellchat_NC <- netAnalysis_computeCentrality(cellchat_NC, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat_NC, signaling = pathways.show, width = 14, height = 2.5, font.size = 10)


netAnalysis_signalingRole_scatter(cellchat_iMCD, signaling = pathways.show)


gg2 = plotGeneExpression(cellchat_iMCD, signaling = pathways.show, enriched.only = TRUE, type = "violin")
ggsave('iMCD_Violin.pdf')
gg3 =plotGeneExpression(cellchat_NC, signaling = pathways.show, enriched.only = TRUE, type = "violin")
ggsave('NC_Violin.pdf')



num.link <- sapply(cellchat.objet.list.new, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link))
gg <- list()
for (i in 1:length(cellchat.objet.list.new)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(cellchat.objet.list.new[[i]], title = names(cellchat.objet.list.new)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)











levels(cellchat@meta$cell_type)
library(ggalluvial)
library(NMF)
# NMF分析输出
selectK(cellchat_iMCD, pattern = "outgoing",k.range = seq(2,6),nrun=64)
nPatterns  = 4
gg5 <- identifyCommunicationPatterns(cellchat_iMCD, pattern = "outgoing", k = nPatterns,width=10,height = 10)
ggsave('iMCD_Pattterns_NMF_outgoing.pdf',width = 10,height = 5)

gg4 = netAnalysis_river(cellchat_iMCD, pattern = "outgoing",font.size = 2)
ggsave('iMCD_NMF_outgoing.pdf',width = 14,height = 5)

netAnalysis_dot(cellchat_iMCD, pattern = "outgoing")
ggsave('iMCD_NMF_dotplot_outgoing.pdf',width = 10,height = 5)


# NMF分析输入
gg6
gg6 = selectK(cellchat_iMCD, pattern = "incoming",k.range = seq(2,6),nrun=64)
ggsave('iMCD_K_incoming.pdf',width = 10,height = 5)
nPatterns  = 3
cellchat_iMCD = identifyCommunicationPatterns(cellchat_iMCD, pattern = "incoming", k = nPatterns,width=10,height = 10)
ggsave('iMCD_Pattterns_NMF_incoming.pdf',width = 10,height = 5)
gg7
gg7 = netAnalysis_river(cellchat_iMCD, pattern = "incoming",font.size = 2)
ggsave('iMCD_NMF_incoming.pdf',width = 14,height = 5)

netAnalysis_dot(cellchat_iMCD, pattern = "incoming")
ggsave('iMCD_NMF_dotplot_incoming.pdf',width = 10,height = 5)


gg8 = selectK(cellchat_NC, pattern = "outgoing",k.range = seq(2,6),nrun=64)
gg8
ggsave('NC_K_outgoing.pdf',width = 10,height = 5)

nPatterns  = 5
cellchat_NC <- identifyCommunicationPatterns(cellchat_NC, pattern = "outgoing", k = nPatterns,width=10,height = 10)
ggsave('NC_Pattterns_NMF_outgoing.pdf',width = 10,height = 5)

gg10 = netAnalysis_river(cellchat_NC, pattern = "outgoing",font.size = 2)
ggsave('NC_NMF_outgoing.pdf',width = 14,height = 5)

netAnalysis_dot(cellchat_NC, pattern = "outgoing")
ggsave('NC_NMF_dotplot_outgoing.pdf',width = 10,height = 5)

gg11
gg11 = selectK(cellchat_NC, pattern = "incoming",k.range = seq(2,6),nrun=64)
ggsave('NC_K_incoming.pdf',width = 10,height = 5)
nPatterns  = 5
cellchat_NC = identifyCommunicationPatterns(cellchat_NC, pattern = "incoming", k = nPatterns,width=10,height = 10)
ggsave('NC_Pattterns_NMF_incoming.pdf',width = 10,height = 5)

gg12 = netAnalysis_river(cellchat_NC, pattern = "incoming",font.size = 2)
ggsave('NC_NMF_incoming.pdf',width = 14,height = 5)


netAnalysis_dot(cellchat_NC, pattern = "incoming")
ggsave('NC_NMF_dotplot_incoming.pdf',width = 10,height = 5)
