library(ggplot2)
library(pheatmap)
library(reshape2)

d <- "/home/molab/Share/hf/Project/ZHR"
setwd(d)
list.files()

mscs <- read.delim("MSCs2.csv",
    header = TRUE,
    sep = ",",
    row.names = 1
)
colnames(mscs)
colnames(mscs) <- c(1)

coln <- colnames(mscs)
rown1 <- rownames(mscs[1:36, ])
rown2 <- rownames(mscs[37:47, ])
rownames(mscs) <- seq_len(nrow(mscs))
colnames(mscs) <- seq_len(ncol(mscs))
mscs <- apply(mscs, 1, scale)
mscs <- t(mscs)


bk <- c(seq(-2, -0.01, by = 0.001), seq(0, 1.5, by = 0.001))

mscs1 <- mscs[1:36, ]
p <- pheatmap(mscs1,
    border = "white", # 设置边框为白色
    breaks = bk,
    color = c(
        colorRampPalette(colors = c("#076AC1", "white"))(length(bk) / 2),
        colorRampPalette(colors = c("white", "#F6507A"))(length(bk) / 2)
    ),
    cluster_cols = F, # 去掉横向、纵向聚类
    cluster_rows = T,
    treeheight_col = 50, # 分别设置横、纵向聚类树高
    treeheight_row = 45,
    labels_row = rown1,
    legend = T, # 添加图例
    legend_breaks = c(-2, -1, 0, 1),
    fontsize_row = 8,
    fontsize_col = 6,
    cellwidth = 20, cellheight = 10, # 设置热图方块宽度和高度
    filename = "P111.jpeg",
    width = 12, height = 8
)


mscs2 <- mscs[37:47, ]
p2 <- pheatmap(mscs2,
    border = "white", # 设置边框为白色
    cluster_cols = F, # 去掉横向、纵向聚类
    cluster_rows = F,
    breaks = bk,
    color = c(
        colorRampPalette(colors = c("#076AC1", "white"))(length(bk) / 2),
        colorRampPalette(colors = c("white", "#F6507A"))(length(bk) / 2)
    ),
    treeheight_col = 50, # 分别设置横、纵向聚类树高
    treeheight_row = 45,
    labels_col = coln,
    labels_row = rown2,
    # legend = F, # 添加图例
    fontsize_row = 8,
    fontsize_col = 6,
    angle_col = 45,
    legend = F,
    cellwidth = 20, cellheight = 10, # 设置热图方块宽度和高度
    filename = "P222.jpeg",
    width = 12, height = 8
)