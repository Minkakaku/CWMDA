rm(list = ls())
setwd("/data/gpfs02/wmo/work/HF/06.ST-seq/project/wk2")

mypvals <- read.table("./out/pvalues.txt", header = T, sep = "\t", stringsAsFactors = F)
mymeans <- read.table("./out/means.txt", header = T, sep = "\t", stringsAsFactors = F)

kp <- grepl(pattern = "Microglial", colnames(mypvals))
table(kp)
pos <- (seq_len(ncol(mypvals)))[kp]
choose_pvalues <- mypvals[, c(c(1, 5, 6, 8, 9), pos)]
choose_means <- mymeans[, c(c(1, 5, 6, 8, 9), pos)]

logi <- apply(choose_pvalues[, 5:ncol(choose_pvalues)] < 0.001, 1, sum)
# 只保留具有细胞特异性的一些相互作用对
choose_pvalues <- choose_pvalues[logi >= 10, ]

# 去掉空值
logi1 <- choose_pvalues$gene_a != ""
logi2 <- choose_pvalues$gene_b != ""
logi <- logi1 & logi2
choose_pvalues <- choose_pvalues[logi, ]

# 同样的条件保留choose_means
choose_means <- choose_means[choose_means$id_cp_interaction %in% choose_pvalues$id_cp_interaction, ]

# 将choose_pvalues和choose_means数据宽转长
library(tidyverse)
meansdf <- choose_means %>% reshape2::melt()
meansdf <- data.frame(
    interacting_pair = paste0(meansdf$gene_a, "_", meansdf$gene_b),
    CC = meansdf$variable,
    means = meansdf$value
)
pvalsdf <- choose_pvalues %>% reshape2::melt()
pvalsdf <- data.frame(
    interacting_pair = paste0(pvalsdf$gene_a, "_", pvalsdf$gene_b),
    CC = pvalsdf$variable,
    pvals = pvalsdf$value
)

# 合并p值和mean文件
pvalsdf$joinlab <- paste0(pvalsdf$interacting_pair, "_", pvalsdf$CC)
meansdf$joinlab <- paste0(meansdf$interacting_pair, "_", meansdf$CC)
pldf <- merge(pvalsdf, meansdf, by = "joinlab")
unique(pldf$interacting_pair.x)
pairs <- unique(unlist(str_split(unique(pldf[, 2]), "_")))
library(Seurat)
wk <- readRDS("finally_wk2_te.rds")
for (i in pairs) {
    p <- SpatialFeaturePlot(wk, features = i, alpha = c(0.1, 1))
    ggsave(paste0(i, "wk2.FeaturePlot.pdf"), plot = p)
}
# dotplot可视化
summary((filter(pldf, means > 0))$means)
head(pldf)
pcc <- pldf %>%
    filter(means > 0) %>%
    ggplot(aes(CC.x, interacting_pair.x)) +
    geom_point(aes(color = means, size = -log10(pvals + 0.0001))) +
    scale_size_continuous(range = c(1, 3)) +
    scale_color_gradient2(high = "red", mid = "yellow", low = "darkblue", midpoint = 1) +
    theme_bw() +
    # scale_color_manual(values = rainbow(100))+
    theme(axis.text.x = element_text(angle = -45, hjust = -0.1, vjust = 0.8))
pcc
ggsave("wk2.interreact.pdf", width = 7, height = 15, plot = pcc)
ggsave("wk2.interreact.png", width = 7, height = 15, plot = pcc)



library(Seurat)
sce <- wk
Idents(sce) <- sce$seurat_clusters
table(Idents(sce))
# sce <- sce[, Idents(sce) %in% "Microglial cell"]
DimPlot(sce, label = T, repel = T)

cg <- unique(unlist(str_split(unique(pldf[, 2]), "_")))

library(ggplot2)
th <- theme(axis.text.x = element_text(angle = -45, hjust = -0.1, vjust = 0.8))
pdopt <- DotPlot(sce,
    features = cg,
    assay = "RNA"
) + coord_flip() + th
ggsave("wk2.features.png", width = 7, height = 23, plot = pdopt)
