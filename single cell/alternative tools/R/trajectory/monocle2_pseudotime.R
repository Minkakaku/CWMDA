# Seurat4Pseudotime
library(Seurat)
library(monocle)
library(dplyr)
library(ggplot2)
library(getopt)
library(reshape2)
library(cowplot)
# rm(list = ls())
# seu_T =readRDS("Neutrophil.rds")
# geneFile = NULL                        #/t分割符的文本
# lambda = 8                             #数值
# RootState = NULL                       #起点位置state数值
# Reverse = F                            #逻辑值
# celllist = NULL
# cluster = NULL
# RandomCells = T
# cellsPercent = 0.1
# minCellsPerCluster = 1000
# seuset = seu_T
options(stringsAsFactors = FALSE,warn = -1)
command = matrix(c('rds','r',1,"character",
                   'gene','g',2,"character",
                   'lambda','l',1,"numeric",
                   'outpath','o',1,"character",
                   'RootState','S',2,"numeric",
                   'Reverse','v',2,"logical",
                   'celllist','C',2,"character",
                   'clusters','c',2,"character",
                   'RandomCells','R',2,"logical",
                   'CellsPercent','p',2,"numeric",
                   'MinCellsPerCluster','M',2,"numeric"
                   )
                 ,byrow=TRUE, ncol=4)
args=getopt(command)
if(is.null(args$gene)){args$gene=NULL}
if(is.null(args$RootState)){args$RootState = NULL}
if(is.null(args$CellsPercent)){args$CellsPercent=0.1}
if(is.null(args$MinCellsPerCluster)){args$MinCellsPerCluster=1000}
if(is.null(args$RandomCells)){args$RandomCells= FALSE}
if(is.null(args$Reverse)){args$Reverse = NULL}
rdsFile = args$rds
geneFile = args$gene;print(paste0("gene use:",head(geneFile)))
outpath = args$outpath;print(paste0("out path:",outpath))
lambda = args$lambda;print(paste0("lambda use:",lambda))
RootState = args$RootState
Reverse = args$Reverse
celllist = args$celllist
cluster = args$clusters;if(!is.null(cluster)){cluster = unlist(strsplit(cluster,","))}
print(paste0("cluster use:",cluster))
RandomCells = args$RandomCells;print(paste0("Random Cell:",RandomCells))
cellsPercent = args$CellsPercent;print(paste0("cells Percent:",cellsPercent))
minCellsPerCluster = args$MinCellsPerCluster;print(paste0("min Cells PerCluster:",minCellsPerCluster))
dir.create(paste0(outpath,"Pseudotime"))
setwd(paste0(outpath,"Pseudotime"))
# rm(list = ls())
# 函数名: Psepre
# 参数: seuset - Seurat对象 celllists -抽取的细胞列表 cluster -选择细胞类型 RandomCells -随机抽取
# 返回值: seurat - Seurat对象
# 功能: 对传入的seuset进行处理并返回
Psepre= function(seuset,celllist,cluster,RandomCells,cellsPercent,minCellsPerCluster){
  # if(!is.null(celllist)){
  #   newcluster = read.table(celllist,sep="\t",header=F)
  #   realcluster = newcluster[is.finite(match(newcluster[,1],colnames(seuset))),]
  #   seuset = subset(seuset,cells=as.character(realcluster[,1]))
  #   seuset = SetIdent(object = seuset, cells= as.character(realcluster[,1]),
  #                     as.character(realcluster[,2]))
  # }else if(!is.null(cluster)){
  #   clusters = unlist(strsplit(cluster, "[,]"))
  #   celllist = names(seuset@active.ident[is.finite(match(seuset@active.ident,clusters))])
  #   seuset = subset(seuset,cells=celllist)
  # }
  if(!is.null(RandomCells)){
    newcelllist = randomSample(seuset@active.ident,cellsPercent = 0.1,minCellsPerCluster = 1000)
    seuset = subset(seuset,cells=as.character(newcelllist))
    write.table(data.frame(Cell = colnames(seuset),Cluster=seuset@active.ident),
                file=paste0("RandomCells.txt"),sep="\t",row.names=F,quot=F)
  }
  seuset <- AddMetaData(object = seuset, metadata = seuset@active.ident, col.name = "Cluster")
  return(seuset)
}
# 函数名: setNum
# 参数: number - 传入细胞数 cellsPercent -抽取的比例 minCellsPerCluster -最小细胞数
# 返回值: finalNum - 需要抽取的细胞数量
# 功能: 对传入的细胞数以及参数进行处理并返回最终抽取细胞数量
setNum = function(number,cellsPercent,minCellsPerCluster){
  selectNum = round(number*cellsPercent)
  if(selectNum>=minCellsPerCluster){
    finalNum = selectNum
  }else{
    finalNum = min(number,minCellsPerCluster)
  }
  return(finalNum)
}
# 函数名: setNum
# 参数: cluster - 细胞类型定义 cellsPercent -抽取的比例 minCellsPerCluster -最小细胞数
# 返回值: finalNum - 需要抽取的细胞数量
# 功能: 对传入的细胞数以及参数进行处理并返回最终抽取细胞数量
randomSample = function(cluster,size=NULL,cellsPercent,minCellsPerCluster){
  samples = c()
  id = unique(cluster)
  for(i in id){
    clusterInfo = cluster[cluster==i]
    if(is.null(size)){
      finalNum = setNum(length(clusterInfo),cellsPercent,minCellsPerCluster)
    }else{
      finalNum = size
    }
    insample = sample(names(clusterInfo),finalNum,replace = FALSE)
    samples = c(samples,insample)
  }
  return(samples)
}
# 函数名: ImportSeurat4CDS
# 参数: otherCDS - Seurat对象
# 返回值: monocle_cds - 转换后的Monocle对象
# 功能: 将Seurat对象转换为Monocle对象，以便在Monocle中进行进一步的分析和可视化。
ImportSeurat4CDS = function(otherCDS){
  # 获取原始数据，并判断是否为稀疏矩阵，不是则进行转换
  # 官方文档推荐加快计算速度
  # otherCDS = seuset
  data = GetAssayData(otherCDS, slot = "counts", assay = "RNA")
  if (class(data) == "data.frame") {
    data = as(as.matrix(data),"sparseMatrix")
  }
  # 获取meta.data,用于monocle定义数据
  pd = tryCatch(
    {
      pd = new("AnnotatedDataFrame", data = otherCDS@meta.data)
      pd
    },error = function(e){
      pData = data.frame(cell_id = colnames(data), row.names = colnames(data))
      pd = new("AnnotatedDataFrame", data = pData)
      message("This Seurat object doesn't provide any meta data")
      pd
    }
  )
  # raw meta不一致以meta为主，用于monocle的featureData
  if (length(setdiff(colnames(data),rownames(pd))) >0) {
    data = data[,rownames(pd)]
  }
  fd = data.frame(gene_short_name = row.names(data),row.names = row.names(data))
  fd = new("AnnotatedDataFrame",data = fd)
  # 这里使用VGAM包的函数定义传入数据的分布类型，其中negbinomial.size()以及negbinomial()适用于UMI数据前者速度更快
  # 通常情况下将会运行第一种，FPKM值适用tobit
  # if (all(data == floor(data))) {
    expressionFamily <- negbinomial.size()
  # }else if(any(data < 0)){
  #   expressionFamily <- uninormal()
  # }else{
  #   expressionFamily <- tobit()
  # }
  # 创建monocle使用object，expressionFamily表达相应的VGAM函数
  monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd,
                                lowerDetectionLimit = 0, expressionFamily = expressionFamily)
  # 如果传入的Seurat包含monocle的数据进行清除
  if ("Monocle" %in% names(otherCDS@misc)) {
    otherCDS@misc$Monocle@auxClusteringData$seurat = NULL
    otherCDS@misc$Monocle@auxClusteringData$scran = NULL
    monocle_cds = otherCDS@misc$Monocle
    mist_list = otherCDS
  }else{
    mist_list = otherCDS
  }
  # 其他信息存储
  monocle_cds@auxClusteringData$seurat <- mist_list
  return(monocle_cds)
}
# 函数名: PsePlot
# 参数: CDS - monocle_cds对象
# 返回值: 画图
# 功能: 画各种图
PsePlot = function(seuCDS){
  gg =  plot_cell_trajectory(seuCDS, color_by = "State")+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.State.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.State.pdf"),limitsize = FALSE)
  gg =  plot_cell_trajectory(seuCDS, color_by = "Cluster")+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.Cluster.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.Cluster.pdf"),limitsize = FALSE)
  gg =  plot_cell_trajectory(seuCDS, color_by = "Pseudotime")+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.Pseudotime.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.Pseudotime.pdf"),limitsize = FALSE)
  gg =  plot_cell_trajectory(seuCDS, color_by = "Cluster") + 
    facet_wrap(~Cluster, nrow = floor(sqrt(length(unique(seuCDS$Cluster)))))+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.ClusterEach.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.ClusterEach.pdf"),limitsize = FALSE)
  gg =  plot_cell_trajectory(seuCDS, color_by = "State") +
    facet_wrap(~State, nrow = floor(sqrt(length(unique(seuCDS$State)))))+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.StateEach.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.StateEach.pdf"),limitsize = FALSE)
  gg = plot_complex_cell_trajectory(seuCDS, color_by = 'Cluster', 
                                    show_branch_points = T, cell_size = 0.5, 
                                    cell_link_size = 0.3)+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.tree.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.tree.pdf"),limitsize = FALSE)
}
# 函数名: PsePlotSample
# 参数: CDS - monocle_cds对象
# 返回值: 画样本图
# 功能: 画样本图
PsePlotSample = function(seuCDS){
  if(("sample"%in% colnames(seuCDS@phenoData@data))|("orig.ident"%in% colnames(seuCDS@phenoData@data))){
    if("orig.ident"%in% colnames(seuCDS@phenoData@data)){
      seuCDS$Sample = as.factor(seuCDS$orig.ident)
    }else {
      seuCDS$Sample = as.factor(seuCDS$sample)     
    }
  gg =  plot_cell_trajectory(seuCDS, color_by = "Sample")+theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.Sample.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.Sample.pdf"),limitsize = FALSE)
  }
  if("group" %in% colnames(seuCDS@phenoData@data)){
  gg =  plot_cell_trajectory(seuCDS, color_by = "group")+theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.group.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.group.pdf"),limitsize = FALSE)
  }
  if("Sample"%in% colnames(seuCDS@phenoData@data)){
  gg =  plot_cell_trajectory(seuCDS, color_by = "Cluster") + facet_wrap(~Sample, nrow = floor(sqrt(length(unique(seuCDS$Sample)))))+theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.SampleEach.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.SampleEach.pdf"),limitsize = FALSE)
  }
  if("group" %in% colnames(seuCDS@phenoData@data)){
  gg =  plot_cell_trajectory(seuCDS, color_by = "Cluster") + facet_wrap(~group, nrow = floor(sqrt(length(unique(seuCDS$group)))))+theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.groupEach.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.groupEach.pdf"),limitsize = FALSE)
  }
  if("Sample"%in% colnames(seuCDS@phenoData@data)){
  gg = plot_complex_cell_trajectory(seuCDS, color_by = 'Cluster', show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~Sample, nrow = floor(sqrt(length(unique(seuCDS$Sample))))) + scale_size(range = c(0.2, 0.2))+theme(legend.position="right",axis.text.x = element_text(angle = 30, hjust = 1))+theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("SampleEach.tree.png"),limitsize = FALSE)
  ggsave(paste0("SampleEach.tree.pdf"),limitsize = FALSE)
  }
  if("group" %in% colnames(seuCDS@phenoData@data)){
  gg = plot_complex_cell_trajectory(seuCDS, color_by = 'Cluster', show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~group, nrow = floor(sqrt(length(unique(seuCDS$Sample))))) + scale_size(range = c(0.2, 0.2))+theme(legend.position="right",axis.text.x = element_text(angle = 30, hjust = 1))+theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("groupEach.tree.png"),limitsize = FALSE)
  ggsave(paste0("groupEach.tree.pdf"),limitsize = FALSE)
  }
}
PseDenPlot = function(forplot){
  if("Sample" %in%colnames(forplot)){
    if(length(unique(forplot$Sample))>1){
    nrow = floor(sqrt(length(unique(forplot$Sample))))
    ncol = ceiling(length(unique(forplot$Sample))/nrow)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=Sample))+geom_density(alpha=0.5)+theme_cowplot()+facet_wrap(~Sample,nrow = nrow)
    ggsave("PseudotimeDensity.SampleEach.png",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE)
    ggsave("PseudotimeDensity.SampleEach.pdf",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=Sample))+geom_density(alpha=0.5)+theme_cowplot()
    ggsave("PseudotimeDensity.Sample.png",gg,width= 8,height=5,limitsize = FALSE)
    ggsave("PseudotimeDensity.Sample.pdf",gg,width= 8,height=5,limitsize = FALSE)
    }else{
    # plotSampledata=ggplot_build(gg)$data[[1]]
    # write.table(plotSampledata,file="../DCPlotForVis/orig.ident.DCPlot.txt",sep="\t",row.names=F,quot=F)
    nrow = floor(sqrt(length(unique(forplot$Cluster))))
    ncol = ceiling(length(unique(forplot$Cluster))/nrow)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=Cluster))+geom_density(alpha=0.5)+theme_cowplot()+facet_wrap(~Cluster,nrow = nrow)
    ggsave("PseudotimeDensity.ClusterEach.png",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE)
    ggsave("PseudotimeDensity.ClusterEach.pdf",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=Cluster))+geom_density(alpha=0.5)+theme_cowplot()
    ggsave("PseudotimeDensity.Cluster.png",gg,width= 8,height=5,limitsize = FALSE)
    ggsave("PseudotimeDensity.Cluster.pdf",gg,width= 8,height=5,limitsize = FALSE)
    # plotdata=ggplot_build(gg)$data[[1]]
    # write.table(plotdata,file="../DCPlotForVis/Cluster.DCPlot.txt",sep="\t",row.names=F,quot=F)   
    }
  }
  if("group" %in%colnames(forplot)){
    nrow = floor(sqrt(length(unique(forplot$group))))
    ncol = ceiling(length(unique(forplot$group))/nrow)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=group))+geom_density(alpha=0.5)+theme_cowplot()+facet_wrap(~group,nrow = nrow)
    ggsave("PseudotimeDensity.groupEach.png",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE)
    ggsave("PseudotimeDensity.groupEach.pdf",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=group))+geom_density(alpha=0.5)+theme_cowplot()
    ggsave("PseudotimeDensity.group.png",gg,width= 8,height=5,limitsize = FALSE)
    ggsave("PseudotimeDensity.group.pdf",gg,width= 8,height=5,limitsize = FALSE)
    
    gg = ggplot(forplot,aes(x=Pseudotime))+geom_density(alpha=0.5,fill="blue")+theme_cowplot()
    ggsave("PseudotimeDensity.png",gg,width= 8,height=5,limitsize = FALSE)
    ggsave("PseudotimeDensity.pdf",gg,width= 8,height=5,limitsize = FALSE)
  }
  
}

mainfun = function(seuset,geneFile,lambda,RootState,Reverse,celllist,cluster,RandomCells,
                   cellsPercent,minCellsPerCluster){
  
  seuset = readRDS(seuset)
  seuset@active.ident = seuset$cell_type_sub
  # seuset@active.ident = seuset$leiden
  seuset = Psepre(seuset,celllist,cluster,RandomCells,cellsPercent,minCellsPerCluster)
  print("Psepre done")
  # 若没有传入基因选择列表则这里获取seuset对象里面的vargene进行使用
  if(!is.null(geneFile)){
    genelist = unique(read.table(geneFile,sep="\t",header=T)[,1])
  }else{
    seuset=FindVariableFeatures(seuset,selection.method ='vst',nfeatures =2000)
    genelist = VariableFeatures(seuset)
    # print(length(genelist))
    # genelist = genelist[1:100]
  }
  # 获取lambda参数
  lambda = lambda*ncol(seuset)
  # 创建monocle对象,标记用于聚类的基因，标化数据，并进行降维,排序，计算拟时间
  # estimateSizeFactors 用于裱花细胞间的差异 lamda越大branch数越小
  # reduceDimension max_components降维空间维数
  seuCDS = ImportSeurat4CDS(seuset) %>% setOrderingFilter(genelist) %>%
    estimateSizeFactors() %>% 
    reduceDimension(reduction_method = "DDRTree",
    lambda=lambda,
    # auto_param_selection = F,    #细胞数量较低的时候需要有这行代码 Error in kmeans(t(Z), K, centers = centers) : initial centers are not distinct
    # ncenter = 2
    ) %>% orderCells()
  print("orderCells done")
  
  # 修改起点
  if(!is.null(RootState) | !is.null(Reverse)){
    seuCDS = orderCells(seuCDS,root_state=RootState,reverse=Reverse)
  }
  table = table(seuCDS$Cluster,seuCDS$State) %>% as.data.frame() %>% 
    dcast(Var1~Var2)
  colnames(table)[1] = "Cluster"
  table %>% write.table("Pseudotime.Statistics.txt",sep = "\t",quote = F)
  seuCDS@phenoData@data %>% data.frame(CellName=row.names(seuCDS@phenoData@data)) %>%
    write.table("Pseudotime.Summary.txt",sep = "\t",quote = F,row.names = F)
  PsePlot(seuCDS)
  if (length(unique(seuCDS$sample))>1|length(unique(seuCDS$orig.ident))>1) {
    PsePlotSample(seuCDS)
  }
  if("sample"%in% colnames(seuCDS@phenoData@data)){
  forplot = data.frame(Sample=seuCDS$sample,Cluster=seuCDS$Cluster,Pseudotime=seuCDS$Pseudotime)
  PseDenPlot(forplot)
  }
  if("orig.ident"%in% colnames(seuCDS@phenoData@data)){
  forplot = data.frame(Sample=seuCDS$orig.ident,Cluster=seuCDS$Cluster,Pseudotime=seuCDS$Pseudotime)
  PseDenPlot(forplot)
  }
  if("group"%in% colnames(seuCDS@phenoData@data)){
  forplot = data.frame(group=seuCDS$group,Cluster=seuCDS$Cluster,Pseudotime=seuCDS$Pseudotime)
  PseDenPlot(forplot)
  }
  # 储存monocle rds
  saveRDS(seuCDS, file = paste0("Pseudotime.monocle.rds"))
}

# rdsFile = "T.rds"
mainfun(rdsFile,geneFile = geneFile,lambda = lambda,RootState = RootState,
        Reverse = Reverse,celllist = celllist,cluster = cluster,
        RandomCells = RandomCells,cellsPercent = cellsPercent,
        minCellsPerCluster = minCellsPerCluster)
# markerlist = FindAllMarkers(Neurds,only.pos = T)
# rm(list = ls())
# library(dplyr)
# marker = markerlist %>% group_by(cluster) %>% top_n(n=-10,wt = avg_log2FC)
# diff_res = differentialGeneTest(cds_subset[marker$gene,],fullModelFormulaStr = "~sm.ns(Pseudotime)")
# sig_gene_name = row.names(subset(diff_res,qval<0.1))
# pdf("heatmap.pdf",width = 12,height = 7)
# heatmap = plot_pseudotime_heatmap(cds_subset[sig_gene_name,],
#                         num_clusters = 3,
#                         cores = 64,
#                         show_rownames = T,
#                         use_gene_short_name = T,
#                         return_heatmap = T)
# dev.off()
