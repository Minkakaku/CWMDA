rm(list = ls())
.libPaths("/home/wzb/R/x86_64-pc-linux-gnu-library/4.2")
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
options(stringsAsFactors = FALSE,warn = -1)
command = matrix(c('outpath','o',1,"character",   #输出地址
                   'seusetinput','s',1,"character", #rds输入地址
                   'grninput','g',1,"character",  #grn结果文件输入
                   'ctxinput','c',1,"character",   #ctx结果文件输入
                   'loominput','l',2,"character")  #loom结果输入
,byrow=TRUE, ncol=4)
args=getopt(command)

tmpdir = "/lvdata/wzb/scRNA/FWSC20240471/2024.11.12/"
args$outpath = paste0(tmpdir,"/SCENIC/SCENIC_res")
args$seusetinput ="/lvdata/wzb/scRNA/FWSC20240471/2024.11.12/h5_trans.rds"
args$grninput = paste0(tmpdir,"/SCENIC/grn_result.csv")
args$ctxinput = paste0(tmpdir,"/SCENIC/ctx_result.csv")
args$loominput = paste0(tmpdir,"/SCENIC/aue_result.loom")
args$cellClusters = paste0(tmpdir,"/SCENIC/meta.xls")
args$umap = paste0(tmpdir,"SCENIC/umap.xls")

packageVersion("SCENIC")
setwd(args$outpath)
# rm(list = ls())
seuset = readRDS(args$seusetinput)
grn_re = read.csv(args$grninput)
ctx_re  = read.csv(args$ctxinput)

loom = open_loom(args$loominput)
  # Read information from loom file:
  exprMat = get_dgem(loom)
  exprMat_log = log2(exprMat + 1) # Better if it is logged/normalized
  regulons_incidMat = get_regulons(loom, column.attr.name="Regulons")
  regulons = regulonsToGeneLists(regulons_incidMat)
  regulonAUC = get_regulons_AUC(loom,column.attr.name = "RegulonsAUC")
  regulonAucThresholds = get_regulon_thresholds(loom)
  embeddings = get_embeddings(loom)
  # cellClusters = seuset@meta.data[, ncol(seuset@meta.data), drop = FALSE]
  # cellClusters$seurat_clusters = as.character(cellClusters$seurat_clusters)
  cellClusters = read.csv(args$cellClusters,row.names = 1,sep = "\t")
  cellClusters$seurat_clusters = paste0("Cluster",cellClusters$leiden)
close_loom(loom)
# ACU值文件
regulonAUC
selectedResolution = "cell_type"
cellsPerCluster  =  split ( rownames ( cellClusters ),  cellClusters [, selectedResolution ])  
regulonAUC  =  regulonAUC [ onlyNonDuplicatedExtended ( rownames ( regulonAUC )),] 
# 计算平均表达式: 
regulonActivity_byCellType  =  sapply ( cellsPerCluster , function ( cells )  rowMeans (getAUC ( regulonAUC )[, cells ])) 
# 尺度表达式：
regulonActivity_byCellType_Scaled  =  t ( scale ( t ( regulonActivity_byCellType ),  center  =  T ,  scale = T ))
regulonActivity_byCellType_Scaled[is.na(regulonActivity_byCellType_Scaled)] = 0
hm = Heatmap(regulonActivity_byCellType_Scaled ,  name = "Regulon Activity",
             row_names_gp = grid :: gpar ( fontsize = 6 ))  # 行字体大小
pdf("Regulon Activity.pdf",width = 20,height = 30)
draw(hm)
dev.off()
regulonOrder  <-  rownames ( regulonActivity_byCellType_Scaled )[row_order( hm )]  # 保存聚集的调节子供以后使用

topRegulators <- melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]

# 为了识别簇特异性调节子（特别是对于多种细胞类型的分析，其中一些调节子对于多种细胞类型来说是通用的），
# 我们发现调节子特异性评分（RSS）特别有用（由 Suo 等人在 2018 年为小鼠细胞图谱提出） ）。
rss  <-  calcRSS ( AUC = getAUC ( regulonAUC ),
                   cellAnnotation = cellClusters [ colnames ( regulonAUC ),  selectedResolution ])
rss[is.na(rss)] = 0
# library(plotly)
rssPlot <- plotRSS(rss,cluster_columns = T)

gg = rssPlot$plot
gg
gg = gg + theme_bw()+theme(panel.grid.major=element_line(colour=NA),axis.text.x=element_text(angle=45,size=8,hjust = 1))
ggsave("RSS.bubble.pdf",gg,height = 30)
ggsave("RSS.bubble.png",gg,height = 30)

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
