library(monocle)
library(getopt)
library(reshape2)
library(dplyr)
library(ggplot2)
library(cowplot)
command = matrix(c('rds','r',1,"character",
                'outpath','o',1,"character",
                'genelist','g',2,"character",
                'ReOrderCell','R',2,"logical",
                'Reverse','s',2,"logical",
                'RootState','S',2,"numeric")
                ,byrow=TRUE, ncol=4)
args=getopt(command)
options(stringsAsFactors = FALSE)
useGenelist = FALSE
if(is.null(args$RootState)){args$RootState = NULL}
if(is.null(args$Reverse)){args$Reverse = FALSE}
if(is.null(args$ReOrderCell)){args$ReOrderCell = FALSE}
if(!is.null(args$genelist)){useGenelist = TRUE}
fplotcolor = unlist(strsplit(("grey,red"),"[,]"))

rdsFile = args$rds;
outpath = args$outpath;
print("Start")
setwd(outpath)
seuCDS=readRDS(rdsFile);
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

if(args$ReOrderCell){
    print("ReOrderCell")
    RootState = args$RootState
    print(paste0('Use start',RootState))
    Reverse = args$Reverse
    seuCDS = orderCells(seuCDS,root_state=RootState,reverse=Reverse)
    saveRDS(seuCDS, file = paste0("ReOrderCells.monocle.pseudotime.rds"))
    phenoData = data.frame(CellName=rownames(seuCDS@phenoData@data),seuCDS@phenoData@data)
    write.table(phenoData,file=paste0("ReOrderCells.Summary_Cell.txt"),sep="\t",row.names=F,quot=F)
    # table = table(seuCDS$FinalCluster,seuCDS$State)
    table = table(seuCDS$Cluster,seuCDS$State)
    print(table)
    table = dcast(data.frame(table),Var1~Var2)
    colnames(table)[1]="Cluster"
    write.table(table,file=paste0("ReOrderCells.Statistics.txt"),sep="\t",row.names=F,quot=F)
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
    
}
