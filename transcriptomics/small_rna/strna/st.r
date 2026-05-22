# 载入所需的R包；
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)

# 可以使用Read10X_h5()或Read10X()读入表达矩阵后创建Seurat对象, 然后，使用Read10X_Image()函数读入空间图像信息文件，整合到Seurat对象上
# 比较快速的方式是将表达量数据和图像信息文件置于同一文件夹，用Load10X_Spatial()函数一次性读入，一步到位
data_dir <- "F:\\spatial_data"
file_name <- "V1_Mouse_Brain_Sagittal_Anterior_filtered_feature_bc_matrix"
# 读入数据并创建Seurat对象
brain <- Load10X_Spatial(data.dir = data_dir, filename = file_name,slice = "anterior1")
brain@project.name <- "anterior1"
Idents(brain) <- "anterior1"
brain$orig.ident <- "anterior1"

# 数据预处理
p1 <- VlnPlot(brain, features = "nCount_Spatial",
pt.size = 0,cols = "tomato") + NoLegend()
p2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") +
theme(legend.position = "right")
p1 | p2

p3 <- VlnPlot(brain, features = "nFeature_Spatial",
pt.size = 0,cols = "tomato") + NoLegend()
p4 <- FeatureScatter(brain, feature1 ="nCount_Spatial",
feature2 = "nFeature_Spatial")+ NoLegend()
p3 | p4

# 由于不同位置的细胞密度可能不一致所以需要标准化
# 推荐采用sctransform方法
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
#提取表达量变变化最高的10个基因;
top10 <- head(VariableFeatures(brain),10)
top10
p5 <- VariableFeaturePlot(brain,cols = c("gray60", "red"))+NoLegend()
p6 <- LabelPoints(plot = p5,points = top10, repel = TRUE,xnudge=0,ynudge=0)
p6

#展示特定基因的表达水平；
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))

#也可以对“点”进行个性化设置：点的大小和透明度设置；
p7 <- SpatialFeaturePlot(brain, features = "Ttr", pt.size.factor = 1.5,stroke = 0.0)
p8 <- SpatialFeaturePlot(brain, features = "Ttr", alpha = c(0.4, 1))
p7 + p8
# pt.size.factor：设置spots的大小，默认为1.6；
# alpha：设置spots的最小和最大透明度值，默认alpha = c(1, 1),表达量越低越透明；
# stroke：设置spots的描边粗细,默认stroke = 0.25

# 降维聚类分析
#接下来是常规单细胞转录组分析流程：PCA降维、UMAP聚类；
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
#绘制UMAP分群图；
p9 <- DimPlot(brain, reduction = "umap", label = TRUE)
p9

#在切片图像上映射分群信息；
p10 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p10

#使用cells.highlight参数单独展示单个cluster的映射情况,这里仅展示1~4个cluster；
cluster <- CellsByIdentities(object = brain, idents = c(1, 2, 3, 4))
SpatialDimPlot(brain, cells.highlight = cluster ,
facet.highlight = TRUE,
cols.highlight = c("green", "grey70"),
ncol = 2)

# 差异分析
#对5、6两个cluster进行差异分析；
de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
head(de_markers)
#绘制前3个差异基因的表达映射图；
SpatialFeaturePlot(object = brain,
features = rownames(de_markers)[1:3],
alpha = c(0.1, 1), ncol = 3)
genes <- VariableFeatures(brain)
#从前100个中查找差异基因；
brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT",
features = genes[1:100],
selection.method = "markvariogram")
DEGs <- SpatiallyVariableFeatures(brain, selection.method = "markvariogram")
top.features <- head(DEGs, 6)
#对差异基因进行可视化；
SpatialFeaturePlot(brain, features = top.features, ncol = 2, alpha = c(0.1, 1))

# 提取子亚群
# 粗略提取额叶皮层的spots数据，然后根据坐标逐步“精修”；
cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
p_img <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE, label.size = 3)
p_img
#获取切片图像坐标和Key：
Key(object =cortex@images$anterior1)
img<- GetTissueCoordinates(cortex)
head(img)
#增加分群信息；
ident <- Idents(cortex)
imgId <- data.frame(img,ident)
head(imgId)

#spots坐标信息的探索；
p_plot <- ggplot(img, aes(x=imagecol, y=600-imagerow,color = ident)) +
geom_point(size = 0.6)
p_img+p_plot

#根据切片图像坐标（前缀需和Key保持一致）进一步去除多余的spots数据；
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
SpatialDimPlot(cortex, crop = TRUE, label = TRUE, label.size = 3)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
SpatialDimPlot(cortex, crop = TRUE, label = TRUE, label.size = 3)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)
SpatialDimPlot(cortex, crop = TRUE, label = TRUE, label.size = 3)


#对切片特定区域的数据进行局部和整体可视化；
p11 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE, label.size = 3)
p12 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE,
pt.size.factor = 1, label.size = 3)
p12 + p11

# 多切片数据整合
data_dir <- "C:/Users/MHY/Desktop/brain2/"
#查看目录中的文件；
list.files(data_dir)
file_name <- "V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.h5"
#读入数据并创建Seurat对象；
brain2 <- Load10X_Spatial(data.dir = data_dir, filename = file_name,slice="Posterior1")
brain2@project.name <- "Posterior1"
Idents(brain2) <- "Posterior1"
brain2$orig.ident <- "Posterior1"

#读入数据并对表达量数据标准化；
brain2 <- SCTransform(brain2, assay = "Spatial", verbose = FALSE)
#合并两个对象；
brain.merge <- merge(brain, brain2)
save(brain.merge,file = "brain.merge_sct.Rdata")

#合并高变基因集，执行常规的PCA、分群聚类流程；
DefaultAssay(brain.merge) <- "SCT"
VariableFeatures(brain.merge) <- c(VariableFeatures(brain), VariableFeatures(brain2))
brain.merge <- RunPCA(brain.merge, verbose = FALSE)
brain.merge <- FindNeighbors(brain.merge, dims = 1:30)
brain.merge <- FindClusters(brain.merge, verbose = FALSE)
brain.merge <- RunUMAP(brain.merge, dims = 1:30)
#对合并后的分群数据进行可视化；
DimPlot(brain.merge, reduction = "umap", group.by = c("ident", "orig.ident"))

#比较两个切片的分群情况；
#默认使用两个slice的名称分组；
SpatialDimPlot(brain.merge)
#如果创建对象时没有对slice重命名，也可使用ident进行分组，绘图效果一样；
SpatialDimPlot(brain.merge,group.by = c("ident"))


#比较两个切片的特定基因表达水平；
SpatialFeaturePlot(brain.merge, features = c("Hpca", "Plp1"))