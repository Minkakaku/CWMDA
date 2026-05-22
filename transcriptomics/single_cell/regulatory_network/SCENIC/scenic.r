# rm(list = ls())
library(SCopeLoomR)
library(AUCell)
library(SCENIC)

# 绘图用包
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)

packageVersion("SCENIC")
setwd("/home/wzb/SCENIC/")
rm(list = ls())
seuset = readRDS("/home/wzb/GSE176078/output/myeloid.rmrb.rds")
grn_re = read.csv("/home/wzb/SCENIC/grn_result.csv")
ctx_re  = read.csv("/home/wzb/SCENIC/ctx_result.csv")

loom <- open_loom("/home/wzb/SCENIC/output/auc_result.loom")
  # Read information from loom file:
  exprMat <- get_dgem(loom)
  exprMat_log <- log2(exprMat + 1) # Better if it is logged/normalized
  regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
  regulons <- regulonsToGeneLists(regulons_incidMat)
  regulonAUC <- get_regulons_AUC(loom,column.attr.name = "RegulonsAUC")
  regulonAucThresholds <- get_regulon_thresholds(loom)
  embeddings <- get_embeddings(loom)
  cellClusters = seuset@meta.data[, ncol(seuset@meta.data), drop = FALSE]
  cellClusters$seurat_clusters = as.character(cellClusters$seurat_clusters)
close_loom(loom)
# ACU值文件
regulonAUC
selectedResolution = "seurat_clusters"
cellsPerCluster  <-  split ( rownames ( cellClusters ),  cellClusters [, selectedResolution ])  
regulonAUC  <-  regulonAUC [ onlyNonDuplicatedExtended ( rownames ( regulonAUC )),] 
# 计算平均表达式: 
regulonActivity_byCellType  <-  sapply ( cellsPerCluster , 
                                         function ( cells )  rowMeans (getAUC ( regulonAUC )[, cells ])) 
# 尺度表达式：
regulonActivity_byCellType_Scaled  <-  t ( scale ( t ( regulonActivity_byCellType ),  center  =  T ,  scale = T ))

hm = Heatmap(regulonActivity_byCellType_Scaled ,  name = "Regulon Activity",
                        row_names_gp = grid :: gpar ( fontsize = 6 ))  # 行字体大小
regulonOrder  <-  rownames ( regulonActivity_byCellType_Scaled )[row_order( hm )]  # 保存聚集的调节子供以后使用

topRegulators <- melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]

# 为了识别簇特异性调节子（特别是对于多种细胞类型的分析，其中一些调节子对于多种细胞类型来说是通用的），
# 我们发现调节子特异性评分（RSS）特别有用（由 Suo 等人在 2018 年为小鼠细胞图谱提出） ）。
rss  <-  calcRSS ( AUC = getAUC ( regulonAUC ),
                   cellAnnotation = cellClusters [ colnames ( regulonAUC ),  selectedResolution ])

# library(plotly)
rssPlot <- plotRSS(rss,cluster_columns = T)

gg = rssPlot$df
gg = gg + theme_bw()+theme(panel.grid.major=element_line(colour=NA),axis.text.x=element_text(angle=45,size=8))
ggsave("RSS.bubble.pdf",gg)
ggsave("RSS.bubble.png",gg)

# 散点图绘制强的转录因子
for (i in unique(colnames(rss))) {
  gg = plotRSS_oneSet(rss, setName = i)  
  ggsave(paste0(i,".Scatterplot.pdf"),gg)
}


umapinfo = seuset@reductions$umap@cell.embeddings
AUCell_plotTSNE(umapinfo,cellsAUC = regulonAUC,plots="AUC",asPNG = T)
binarizeAUC <- function(auc, thresholds)
{
  thresholds <- thresholds[intersect(names(thresholds), rownames(auc))]
  regulonsCells <- setNames(lapply(names(thresholds), 
                                   function(x) {
                                     trh <- thresholds[x]
                                     names(which(getAUC(auc)[x,]>trh))
                                   }),names(thresholds))
  
  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  class(binaryRegulonActivity) <- "matrix"  
  
  return(binaryRegulonActivity)
}
revnum <- function(numthe){
  val = as.numeric(names(numthe))
  name = as.character(numthe)
  numthe = val
  names(numthe) = name
  return(numthe)
}
binaryRegulonActivity = binarizeAUC(regulonAUC,revnum(regulonAucThresholds))
gg = Heatmap(binaryRegulonActivity, name="Binarized activity", 
        col = c("white", "black"),
        cluster_rows = TRUE, cluster_columns=TRUE,
        show_column_names = F,
        row_names_gp=grid::gpar(fontsize=3),use_raster = TRUE) # row font size
getwd()
pdf("Binarizedheatmap.pdf")
gg
dev.off()
