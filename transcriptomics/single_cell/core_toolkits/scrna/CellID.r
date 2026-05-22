library(CellID)
library(Seurat)
library(tidyverse)

#获取原始表达矩阵
BaronMatrix   <- readRDS(url("https://storage.googleapis.com/cellid-cbl/BaronMatrix.rds"))
#仅考虑编码蛋白质的基因
data("HgProteinCodingGenes")
BaronMatrixProt <- BaronMatrix[rownames(BaronMatrix) %in% HgProteinCodingGenes,]

#这几步类似Seurat的标准流程
Baron <- CreateSeuratObject(counts = BaronMatrixProt, project = "Baron", min.cells = 5)
Baron <- NormalizeData(Baron)
Baron <- ScaleData(Baron, features = rownames(Baron)) 
Baron <- RunMCA(Baron) #该软件的主要分析函数，将细胞和基因同时降维到一个空间，离细胞近的基因被定义为细胞的特征基因——令人窒息的操作

#下载参考基因集
panglao <- read_tsv("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz")
#根据器官过滤，示例数据就是胰腺的
panglao_pancreas <- panglao %>% filter(organ == "Pancreas")
#物种包含人
panglao_pancreas <- panglao_pancreas %>%  filter(str_detect(species,"Hs"))
#下面两步会得到一个列表，列表的每一个元素是由基因集构成的字符向量，且每个元素被命名为对应的细胞类型的名字
panglao_pancreas <- panglao_pancreas %>%  
  group_by(`cell type`) %>%  
  summarise(geneset = list(`official gene symbol`))
pancreas_gs <- setNames(panglao_pancreas$geneset, panglao_pancreas$`cell type`)

#富集分析，用的超几何检验
HGT_pancreas_gs <- RunCellHGT(Baron, pathways = pancreas_gs, dims = 1:50, n.features = 200) #每个细胞选200个特征基因
HGT_pancreas_gs[1:4,1:4]
pancreas_gs_prediction <- rownames(HGT_pancreas_gs)[apply(HGT_pancreas_gs, 2, which.max)]
#矩阵的每一列依次：判断是否是最大值，得到一列布尔值，结合矩阵的行名会返回行名中的一个元素（也就是最有可能的细胞类型）。
#所有列运行完了之后会得到所有细胞最可能的注释
pancreas_gs_prediction_signif <- ifelse(apply(HGT_pancreas_gs, 2, max)>2, yes = pancreas_gs_prediction, "unassigned")
#如果`-log10 corrected p-value`的值小于等于2，则认为不显著，注释更正为"unassigned"
Baron$pancreas_gs_prediction <- pancreas_gs_prediction_signif
head(Baron@meta.data,2)

#加载预先知道的注释信息
BaronMetaData <- readRDS(url("https://storage.googleapis.com/cellid-cbl/BaronMetaData.rds"))
Baron@meta.data=merge(Baron@meta.data,BaronMetaData,by="row.names")
rownames(Baron@meta.data)=Baron@meta.data$Row.names
Baron@meta.data$Row.names=NULL
head(Baron@meta.data,2)

Baron <- FindVariableFeatures(Baron, selection.method = "vst", nfeatures = 2000)
Baron <- ScaleData(Baron)
Baron <- RunPCA(Baron, npcs = 50, verbose = FALSE)
Baron <- FindNeighbors(Baron, dims = 1:30)
Baron <- FindClusters(Baron, resolution = 0.5)
Baron <- RunUMAP(Baron, dims = 1:30)
Baron <- RunTSNE(Baron, dims = 1:30)

library(cowplot)
p1 <- DimPlot(Baron, reduction = "tsne", group.by = "cell.type", pt.size=0.5, label = TRUE,repel = TRUE)
p2 <- DimPlot(Baron, reduction = "tsne", group.by = "pancreas_gs_prediction", pt.size=0.5, label = TRUE,repel = TRUE)
fig_tsne <- plot_grid(p1, p2, labels = c('cell.type','pancreas_gs_prediction'),align = "v",ncol = 2)
ggsave(filename = "tsne.pdf", plot = fig_tsne, device = 'pdf', width = 40, height = 15, units = 'cm')