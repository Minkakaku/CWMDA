rm(list=ls())
# 获取当前工作目录
getwd()

# 设置工作目录，若需修改可更改路径
workingDir <- "."
setwd(workingDir)

# 加载 WGCNA 包
library(WGCNA)

# 重要设置，请勿省略
options(stringsAsFactors = FALSE)

# 加载第一部分保存的表达和性状数据
lnames <- load(file = "WGCNA\\opt_inputdata\\opt_inputdata.RData")
# 输出加载变量的名称
lnames

# 加载第二部分保存的网络数据
lnames <- load(file = "WGCNA\\results\\networkConstruction_auto.RData")
lnames

# 获取基因数量和样本数量
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

# 重新计算拓扑重叠矩阵
dissTOM <- 1 - TOMsimilarityFromExpr(datExpr, power = 6)

# 转换 dissTOM 以便在热图中更清晰显示中等强度连接
plotTOM <- dissTOM^7
# 设置对角线为 NA 使绘图更美观
diag(plotTOM) <- NA

# 打开绘图窗口并绘制全基因网络热图
sizeGrWindow(9, 9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

# 设置随机选择的基因数量
nSelect <- 400
# 设置随机种子以保证结果可重复
set.seed(10)
select <- sample(nGenes, size = nSelect)
selectTOM <- dissTOM[select, select]

# 重新聚类选择的基因
selectTree <- hclust(as.dist(selectTOM), method = "average")
selectColors <- moduleColors[select]

# 转换选择基因的 TOM 矩阵并设置对角线为 NA
plotDiss <- selectTOM^7
diag(plotDiss) <- NA

# 打开绘图窗口并绘制选择基因的网络热图
sizeGrWindow(9, 9)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

# 重新计算模块特征基因
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes

# 从临床性状中提取体重数据
weight <- as.data.frame(datTraits$weight_g)
names(weight) <- "weight"

# 将体重数据添加到现有模块特征基因中
MET <- orderMEs(cbind(MEs, weight))

# 绘制特征基因与性状的关系图
sizeGrWindow(5, 7.5)
par(cex = 0.9)
plotEigengeneNetworks(
  MET, "", 
  marDendro = c(0, 4, 1, 2), 
  marHeatmap = c(3, 4, 1, 2), 
  cex.lab = 0.8, 
  xLabelsAngle = 90
)

# 绘制特征基因树状图
sizeGrWindow(6, 6)
par(cex = 1.0)
plotEigengeneNetworks(
  MET, "Eigengene dendrogram",
  marDendro = c(0, 4, 2, 0),
  plotHeatmaps = FALSE
)

# 绘制特征基因邻接热图（会覆盖树状图）
par(cex = 1.0)
plotEigengeneNetworks(
  MET, "Eigengene adjacency heatmap", 
  marHeatmap = c(3, 4, 2, 2),
  plotDendrograms = FALSE, 
  xLabelsAngle = 90
)
