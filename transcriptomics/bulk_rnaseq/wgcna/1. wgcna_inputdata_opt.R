# 设置工作环境
rm(list = ls())
options(stringsAsFactors = FALSE)

# 加载库
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("WGCNA", quietly = TRUE))
  BiocManager::install("WGCNA")
if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr")

# 读取表达数据，明确指定参数
femData <- read.csv("LiverFemale3600.csv", header = TRUE, stringsAsFactors = FALSE)

# 提取表达矩阵，使用dplyr和tidyr进行更清晰的数据处理
datExpr0 <- femData %>% 
  select(-(1:8)) %>%  # 去除不包含表达量的列
  t() %>%  # 转置
  as.data.frame() # 转化为数据框

colnames(datExpr0) <- femData$substanceBXH # 基因名
rownames(datExpr0) <- colnames(femData)[-(1:8)] # 样本名

# 检查数据质量
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
# 计算每个基因的方差（按列计算，因为列是基因）
gene_vars <- apply(datExpr0, 2, var, na.rm = TRUE)

# 按方差从小到大排序并保留后75%基因（列）
filtered_genes <- names(sort(gene_vars))[(ceiling(length(gene_vars) * 0.75)):length(gene_vars)]

good_genes <- colnames(datExpr0) %in% filtered_genes
bad_genes <- sum(!good_genes)
# 计算每个样本的缺失值比例（按行计算，因为行是样本）
sample_na_ratio <- rowMeans(is.na(datExpr0))

# 判断样本是否合格（缺失值 <= 50%）
good_samples <- sample_na_ratio <= 0.5
bad_samples <- sum(!good_samples)

# 统计数据问题
initial_dims <- dim(datExpr0)

# 移除不良数据
if (!gsg$allOK) {
  if (bad_genes > 0) {
    cat(paste0("移除 ", bad_genes, " 个低质量基因\n"))
    datExpr0 <- datExpr0[, good_genes]
  }
  
  if (bad_samples > 0) {
    cat(paste0("移除 ", bad_samples, " 个低质量样本\n"))
    datExpr0 <- datExpr0[good_samples, ]
  }
  
  cat(paste0("数据清理完成: 从 ", initial_dims[1], "x", initial_dims[2], 
             " 减少到 ", dim(datExpr0)[1], "x", dim(datExpr0)[2], "\n"))
}

# 异常样本检测与移除
sampleTree <- hclust(dist(datExpr0), method = "average")

# 创建图形目录
if (!dir.exists("figures")) {
  dir.create("figures")
}

# 保存样本聚类图
png("figures/sample_clustering.png", width = 1500, height = 900, res = 100)
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "样本聚类以检测异常值", sub = "", xlab = "", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = 15, col = "red")
dev.off()

# 确定要保留的样本
clust <- cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
print(table(clust))

keepSamples <- (clust == 1)
datExpr <- datExpr0[keepSamples, ]
cat(paste0("保留 ", sum(keepSamples), " 个样本，移除 ", sum(!keepSamples), " 个异常样本\n"))

# 读取临床特征数据
traitData <- read.csv("ClinicalTraits.csv", header = TRUE, stringsAsFactors = FALSE)
allTraits <- traitData %>% 
  select(-c(31, 16)) %>% 
  select(c(2, 11:36))
# 检查数据
head(allTraits)
# 匹配样本与特征
femaleSamples <- rownames(datExpr)
traitRows <- match(femaleSamples, allTraits$Mice)

# 检查匹配情况
if (any(is.na(traitRows))) {
  cat(paste0("警告: 无法匹配 ", sum(is.na(traitRows)), " 个样本\n"))
  femaleSamples[is.na(traitRows)]
}

datTraits <- allTraits[traitRows, -1]
rownames(datTraits) <- allTraits[traitRows, 1]

# 清理内存
gc()

# 重新聚类样本
sampleTree2 <- hclust(dist(datExpr), method = "average")
traitColors <- numbers2colors(datTraits, signed = FALSE)

# 保存样本树与特征热图
png("figures/sample_dendrogram_traits.png", width = 1500, height = 1200, res = 100)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "样本树状图与特征热图")
dev.off()

# 保存处理后的数据
output_dir <- "opt_inputdata"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

save(datExpr, datTraits, file = file.path(output_dir, "opt_inputdata.RData"))
cat("数据处理完成！结果保存在", output_dir, "目录下\n")    
