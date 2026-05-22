rm(list=ls())
# 加载WGCNA包并进行必要设置
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
# 加载数据
lnames <- load("WGCNA\\opt_inputdata\\opt_inputdata.RData")
lnames

# 选择一组软阈值幂次
powers <- c(1:10, seq(12, 20, by = 2))

# 调用网络拓扑分析函数
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# 绘制结果
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 <- 0.9

# 绘制无标度拓扑拟合指数
colnames(sft$fitIndices)
with(sft$fitIndices, {
  png("WGCNA\\figures/scaleIndependence.png")
  plot(Power, -sign(slope) * SFT.R.sq,
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit,signed R^2",
       type = "n",
       main = "Scale independence"
  )
  text(Power, -sign(slope) * SFT.R.sq, labels = Power, cex = cex1, col = "red")
  dev.off()
})
abline(h = 0.90, col = "red")

# 绘制平均连接性
with(sft$fitIndices, {
  png("WGCNA\\figures/meanConnectivity.png")
  plot(Power, mean.k.,
       xlab = "Soft Threshold (power)",
       ylab = "Mean Connectivity",
       type = "n",
       main = "Mean connectivity"
  )
  text(Power, mean.k., labels = Power, cex = cex1, col = "red")
  dev.off()
})

# 选择软阈值并计算邻接矩阵
softPower <- 6
adjacency <- adjacency(datExpr, power = softPower)

# 将邻接矩阵转换为拓扑重叠矩阵
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# 进行层次聚类
geneTree <- hclust(as.dist(dissTOM), method = "average")

# 绘制聚类树
sizeGrWindow(12, 9)
png("WGCNA\\figures/geneClustering.png")
plot(geneTree, xlab = "", sub = "", 
     main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

# 设置最小模块大小
minModuleSize <- 30

# 使用动态树切割识别模块
dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = minModuleSize
)
table(dynamicMods)

# 将数值标签转换为颜色
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

# 绘制树状图和模块颜色
sizeGrWindow(8, 6)
png("WGCNA\\figures/dynamicTreeCut.png")
plotDendroAndColors(
  geneTree, dynamicColors, "Dynamic Tree Cut",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)
dev.off()

# 计算模块特征基因
MEs <- moduleEigengenes(datExpr, colors = dynamicColors)$eigengenes

# 计算模块特征基因的相异度
MEDiss <- 1 - cor(MEs)

# 对模块特征基因进行聚类
METree <- hclust(as.dist(MEDiss), method = "average")

# 绘制聚类结果
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

# 设置模块合并阈值
MEDissThres <- 0.25

# 绘制切割线
abline(h = MEDissThres, col = "red")

# 调用自动合并函数
merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

# 绘制合并前后的模块颜色
sizeGrWindow(12, 9)
png("figures/mergedTreeCut.png")
plotDendroAndColors(
  geneTree, cbind(dynamicColors, mergedColors),
  c("Dynamic Tree Cut", "Merged dynamic"),
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05
)
dev.off()

# 重命名模块颜色
moduleColors <- mergedColors

# 构建与颜色对应的数值标签
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder) - 1
MEs <- mergedMEs

# 保存模块颜色和标签供后续使用
save(MEs, moduleLabels, moduleColors, geneTree, file = "WGCNA\\results/networkConstruction_man.RData")
