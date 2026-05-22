#PCA plot
library("ggplot2")
library(stats)
library(ggrepel)
library(RColorBrewer)


## yxx ph
directory="D:/Users/24432/Desktop/����/pheatmap/yxx/ph"
setwd(directory)
list.files()

countdata <- read.csv("all_compare.csv", header = T, skip = 0, stringsAsFactors = F)
countdata1 <- countdata[,1:6]
data1 <- as.matrix(sapply(countdata1, as.numeric))
data2 <- countdata1[which(rowSums(countdata1) > 0),]
#pca analysis
data.pca <- princomp(data2, cor = T, scores = T)
summary(data.pca)
proportion = data.pca$sdev^2 / sum(data.pca$sdev^2)
##ggplot2
prinl = loadings(data.pca)
dd = as.data.frame.matrix(prinl)
dd$Sample = unlist(lapply(rownames(dd), FUN = function(x){
  strsplit(x, split = ".FP")[[1]][1]
}))
xlab <- paste("PC1","(",round(proportion[[1]]*100,1),"%)",sep="")
ylab <- paste("PC2","(",round(proportion[[2]]*100,1),"%)",sep="")

p <- ggplot(dd,aes(x=-dd$Comp.1,y=-dd$Comp.2,color=dd$Sample)) +  #
  geom_point(aes(color = dd$Sample), size = 2, alpha = 0.8) + #
  #ggrepel::geom_text_repel(mapping = aes(label=dd$Sample)) +
  scale_color_manual(guide=guide_legend(title="Sample"), 
                     breaks = dd$Sample,
                     values = c("KO_D0_fpkm" = "green","KO_D2_fpkm" = "gray","KO_D7_fpkm" = "black",
                                "WT_D0_fpkm" = "blue","WT_D2_fpkm" = "yellow","WT_D7_fpkm" = "red"))+
  #scale_fill_manual(values = colorRampPalette(brewer.pal(14, "Paired"))(14))+
  theme_bw()+
  #geom_label_repel(aes(x=dd$Comp.1, y=dd$Comp.2, fill = factor(dd$Sample),label=dd$Sample),
  #                fontface="bold", color="white", box.padding=unit(0.15, "lines"), 
  #               point.padding=unit(0.5, "lines"), segment.colour = "grey50", force = 2)+
  theme(plot.background=element_blank(),
        panel.background=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.title=element_text(color="black",size=11),
        axis.text=element_text(size=9),
        legend.title = element_text(face = "bold"),
        legend.key.height = unit(1,"line"))+
  labs(x = xlab, y = ylab)
p
pdf("PCA_yxx.pdf",width=5.5,height=2.5)
print(p)
dev.off()




