# http://www.ehbio.com/ImageGP/index.php/Home/Index/Volcanoplot.html
rm(list = ls())
suppressPackageStartupMessages({
  library(ggpubr)
  library(ggrepel)
  library(gg.gap)
  library(ggplot2)
  library(dplyr)
})

directory <- "/home/hongfan/todolist/todo/220906_CQ_graduate"
setwd(directory)
list.files()

res <- read.csv("siSetdb1_1vssiNC/siSetdb1_1vssiNC_deg.csv",
  header = T,
  stringsAsFactors = F
)
res1 <- res[, c(7, 4, 5)]
colnames(res1) <- c("gene_id", "log2FC", "Pvalue")
res2 <- res1[complete.cases(res1), ]
head(res2)

res3 <- as.data.frame(sapply(res2[, 2:3], as.numeric))
res3$gene_id <- res2$gene_id
# res3$label <- ifelse(res3$Padj < 1 & abs(res3$log2FC) > 1, res3$gene_id, "")
res3$label <- ifelse(res3$gene_id %in% c("Ucp1", "Elovl3", "Dio2", "Ppargc1a", "Elovl6"), res3$gene_id, "")

res3$change <- as.factor(ifelse(res3$Pvalue < 0.05 & abs(res3$log2FC) > 1,
  ifelse(res3$log2FC > 1, "UP", "DOWN"), "NOT"
))

res3$geneid <- res2$gene_id
res3$minuslog10Pvalue <- -log10(res3$Pvalue)
count <- table(res3$change)
a <- ggscatter(res3,
  y = "minuslog10Pvalue",
  x = "log2FC",
  ylab = "-log10(Pvalue)",
  color = "change",
  size = 2,
  # label = "label",
  # font.label = c(12, "#000000"),
  # repel = FALSE,
  palette = c("#00AFBB", "#999999", "#FC4E07"),
  alpha = 0.6
)
a
a2 <- a + annotate("label", x = c(-5, 5), y = 10, label = count[c(1, 3)], col = c("#00AFBB", "#FC4E07")) +
  annotate("text", x = c(-5, 5), y = 10.5, label = c("DOWN", "UP")) +
  geom_text_repel(aes(label = label),
    size = 4, color = "black",
    point.padding = 0.3,
    fontface = "bold",
    hjust = 0.5,
    arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
    segment.size = 0.1, segment.alpha = 0.8, nudge_y = 0.7,
    max.overlaps = 10000000 # , xlim = c(0, 3.5)
  )
a3 <- a2 + ylab("-log10 (pvalue)") + # 修改y轴名称
  xlab("log2 (FoldChange)") + # 修改x轴名称
  geom_vline(xintercept = c(-1, 1), lty = 3, col = "black", lwd = 0.5) + # 添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(0.05), lty = 3, col = "black", lwd = 0.5) +
  theme(
    legend.title = element_blank() # 不显示图例标题
  )
a3
ggsave("siSetdb1_1vssiNC.pdf", width = 10, height = 10)
