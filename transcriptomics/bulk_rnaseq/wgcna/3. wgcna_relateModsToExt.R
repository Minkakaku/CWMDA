# 显示当前工作目录
getwd()

# 设置工作目录
# 如果需要，修改下面的路径为数据文件存储的目录
# "." 表示当前目录。在 Windows 系统中使用正斜杠 / 而不是反斜杠 \
working_dir = "."
# 原代码中变量名拼写错误，修正为 working_dir
setwd(working_dir)

# 加载 WGCNA 包
library(WGCNA)

# 设置重要参数，请勿省略
options(stringsAsFactors = FALSE)

# 加载第一部分保存的表达数据和表型数据
lnames <- load(file = "WGCNA/opt_inputdata/opt_inputdata.RData")
# 变量 lnames 包含已加载变量的名称
lnames

# 加载第二部分保存的网络数据
lnames <- load(file = "WGCNA\\results\\networkConstruction_auto.RData")
lnames

# 定义基因和样本的数量
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# ========== 模块与表型关联分析 ==========
# 重新计算带有颜色标签的模块特征基因 (MEs)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
MEs
# 计算模块特征基因与表型的相关性
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# 设置绘图窗口大小
sizeGrWindow(10, 6)

# 准备要显示的相关性和 p 值文本矩阵
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

# 设置绘图边距
par(mar = c(6, 8.5, 3, 3))

# 绘制热图展示模块与表型的相关性
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships"))

# ========== 基因模块成员关系与基因表型显著性分析 ==========
# 定义 weight 变量，包含 datTraits 中的 weight_g 列
weight = as.data.frame(datTraits$weight_g)
names(weight) = "weight"

# 获取模块名称（颜色）
modNames = substring(names(MEs), 3)

# 计算基因与模块的相关性（模块成员关系）
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

# 设置基因模块成员关系数据框的列名
names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")

# 计算基因与体重表型的相关性（基因表型显著性）
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

# 设置基因表型显著性数据框的列名
names(geneTraitSignificance) = paste("GS.", names(weight), sep = "")
names(GSPvalue) = paste("p.GS.", names(weight), sep = "")

# 选择棕色模块进行分析
module = "brown"
column = match(module, modNames)
moduleGenes = moduleColors == module

# 设置绘图窗口大小
sizeGrWindow(7, 7)
par(mfrow = c(1, 1))

# 绘制散点图展示棕色模块的模块成员关系与基因表型显著性
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# 查看基因表达数据的列名
names(datExpr)

# 查看棕色模块的基因名称
names(datExpr)[moduleColors == "brown"]