#MFuzz

## fig4
library(Mfuzz)

directory = "D:/Users/24432/Desktop/生信/mfuzz时间序列/csx/fig4"
setwd(directory)
list.files()
countdata <- read.csv("sheet3.csv",header = T,stringsAsFactors = F) ##sheet4
countdata[countdata == "--"] <- NA
countdata[countdata == ""] <- NA
data1 <- countdata[complete.cases(countdata),]
data2 <- data1[-which(data1$Other.Gene.ID == "1-Mar"),]
data2 <- data2[-which(data2$Other.Gene.ID == "2-Mar"),]
data4 <- as.data.frame(sapply(data2[,2:20],as.numeric))  #2:15
## 
rownames(data4) <- data2$Other.Gene.ID
data4[data4 == "0"] <- NA
data4 <- data4[complete.cases(data4),]

data5 <- new("ExpressionSet",exprs = as.matrix(data4))
data.r <- filter.NA(data5, thres = 0.25)
data.y <- fill.NA(data.r, mode = "mean")
tmp <- filter.std(data.y, min.std = 0)
data.s <- standardise(data.y)
#m1 <- mestimate(data.s)
#算聚类个数c1（两个方法选一个）
#cselection(data.s, m1, crange = seq(4,40,4),repeats = 3, visu = T)
#Dmin(data.s,m1,seq(3,10,1))
c <- 12
#聚类
cl <- mfuzz(data.s, c = c, m = 1.25)
O = overlap(cl)
overlap.plot(cl, O, thres = 0.05)
pdf("mfyzz(1-60).pdf",width = 9,height = 6) #1-60
mfuzz.plot2(data.s, cl=cl, mfrow=c(3,4), xlab = "", ylab = "",x11 = F,
           time.labels = c("P1","P4","P7","P10","P13","P16","P30","P60")) ##
dev.off()

library(stringr)
tf <- read.csv("mouse_TF.csv",header = T,stringsAsFactors = F)
tf$Gene <- paste0(str_sub(tf$Gene_symbol,1,1), 
                         str_to_lower(str_sub(tf$Gene_symbol,2,str_length(tf$Gene_symbol))))
cl.class <- as.data.frame(cl$cluster)
cl.class$gene <- rownames(cl.class)
all <- merge(x = cl.class, y = tf, by.x = "gene", by.y = "Gene", all.x = T)
all.exp <- merge(x = all, y = data1, by.x = "gene", by.y = "Other.Gene.ID", all.x = T)
write.csv(all.exp,file = "tf_1-60.csv") ##1-16






## test
library(Mfuzz)
directory = 'D:/Users/24432/Desktop/生信/基因表达时间序列图'
setwd(directory)
data1 <- read.table("kegg_cellcycle_allgene_BGfpkm.txt",sep = "\t",header = T)
row.names(data1) <- data1$Other.Gene.ID
data2 <- data1[,-1]
data3 <- data.frame(P1BGwt = rowMeans(data2[,1:2]),P4BGwt = rowMeans(data2[,3:4]),
                    P7BGwt = rowMeans(data2[,5:8]),P10BGwt = rowMeans(data2[,9:10]),
                    P13BGwt = rowMeans(data2[,11:12]),P16BGwt = rowMeans(data2[,13:14]),
                    X1MBGwt = rowMeans(data2[,15:16]),X2MBGwt = rowMeans(data2[,17:19]))
data3[data3==0] <- NA
data4 <- new("ExpressionSet",exprs = as.matrix(data3))
data.r <- filter.NA(data4, thres = 0.25)
data.y <- fill.NA(data.r, mode = "mean")
tmp <- filter.std(data.y, min.std = 0)
data.s <- standardise(data.y)
m1 <- mestimate(data.s)
#算聚类个数c1（两个方法选一个）
#cselection(data.s, m1, crange = seq(4,40,4),repeats = 3, visu = T)
Dmin(data.s,m1,seq(3,10,1))
c <- 6
#聚类
c1 <- mfuzz(data.s, c = c, m = m1)
O = overlap(c1)
overlap.plot(c1, O, thres = 0.05)
mfuzz.plot(data.s,cl=c1,mfrow=c(3,2),new.window = F)
## heatmap
data5 <- data3
data5[is.na(data5)] <- 0
b = as.data.frame(c1$cluster)
b$Gene = rownames(b)
colnames(b) = c("Cluster","Gene")
b = b[,c(2,1)]
data5$Gene = rownames(data5)
result = merge(b,data5,by="Gene",all=F)
result = result[order(result$Cluster,decreasing = F),]
annRow = paste0("cluster",result$Cluster)
annRow = as.data.frame(annRow)
rownames(annRow) = rownames(result)
colnames(annRow) = "Cluster"
library(pheatmap)
pheatmap(result[,c(3:ncol(result))],scale = "row", 
         cluster_rows = F,cluster_cols = F,
         gaps_row = cumsum(c1$size),
         annotation_row = annRow,
         color = colorRampPalette(colors = c("#0569cd","white","#ff0200"))(100),
         show_rownames = F#,annotation_col = annCol
)
ord <- c(4,6,5,1,3,2)
result$Cluster <- factor(result$Cluster, levels = ord)
result <- result[order(result$Cluster),]
pheatmap(result[,c(3:ncol(result))],scale = "row", 
         cluster_rows = F,cluster_cols = F,
         gaps_row = c(95,110,122),
         annotation_row = annRow,
         color = colorRampPalette(colors = c("#0569cd","white","#ff0200"))(100),
         show_rownames = F#,annotation_col = annCol
)
dev.off()
save.image(file='env_cellcycleheatmap1.RData')

