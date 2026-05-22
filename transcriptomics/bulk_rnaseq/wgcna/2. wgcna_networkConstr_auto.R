rm(list=ls())
# 显示当前工作目录
cat("当前工作目录: ", getwd(), "\n")

# 设置工作目录，可根据需要修改路径。"." 表示当前目录，Windows 下使用正斜杠 /
workingDir = "."
setwd(workingDir)

# 加载 WGCNA 包
library(WGCNA)

options(stringsAsFactors = FALSE)

# 启用 WGCNA 多线程，加速部分计算
# 注意：若使用 RStudio 或其他第三方 R 环境，请跳过此命令
# 当前此命令为代码运行所必需，若出现错误可忽略，但建议更新 WGCNA 包
if (!Sys.getenv("RSTUDIO") == "1") {
  enableWGCNAThreads()
}

# 加载第一部分保存的数据
lnames <- load(file = "V:\\repo\\bio-coding\\WGCNA\\opt_inputdata\\opt_inputdata.RData")
# 输出加载的变量名
cat("已加载的变量名: ", paste(lnames, collapse = ", "), "\n")

# 检查 datExpr 是否加载成功
if (!"datExpr" %in% lnames) {
  stop("未成功加载 datExpr 数据，请检查 'opt_inputdata.RData' 文件")
}

# 选择一组软阈值幂次
powers <- c(1:20)
powers
# 调用网络拓扑分析函数
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# 打开绘图窗口并设置布局
sizeGrWindow(9, 10)
par(mfrow = c(1, 2))
cex1 <- 0.9

# 绘制无标度拓扑拟合指数随软阈值幂次的变化
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
# 绘制 R^2 阈值线
abline(h = 0.80, col = "red")

# 绘制平均连接度随软阈值幂次的变化
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
abline(h = 20, col = "red")

optimal_power <- sft$powerEstimate
cat("选择的最优软阈值幂次为: ", optimal_power, "\n")

output_dir <- "results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
# 网络模块构造
net <- blockwiseModules(datExpr,
                         # == Adjacency Function ==
                         power = optimal_power, networkType = "signed",
                         # == Tree and Block Options ==
                         deepSplit = 2, minModuleSize = 30,maxBlockSize = 4000,
                         # == Module Adjustments ==
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         # == TOM == Archive the run results in TOM file (saves time)
                         saveTOMs = TRUE,TOMType = "signed", 
                         saveTOMFileBase = "WGCNA\\results/femaleMouse_auto_singleTOM", 
                         verbose = 3)

# 打开绘图窗口
sizeGrWindow(12, 9)

# 将标签转换为颜色以便绘图
mergedColors <- labels2colors(net$colors)

# 绘制树状图并在下方显示模块颜色
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

# 保存分析结果
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "WGCNA\\results/networkConstruction_auto.RData")
