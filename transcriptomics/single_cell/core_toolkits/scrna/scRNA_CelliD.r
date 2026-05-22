# PLOT
rm(list = ls())
library(Seurat)
library(CelliD)
library(Seurat)
library(tidyverse)
library(homologene)
setwd("E:\\CODE\\PJ_For_PR\\hep")
options(radian.auto_match = FALSE)

# load data
df <- readRDS("hep_integrated_sctv2.rds")
data("MgProteinCodingGenes")

# subset ProteinCodingGenes
df <- df[rownames(df) %in% MgProteinCodingGenes, ]
# df <- NormalizeData(df)
# df <- ScaleData(df, features = rownames(df))
df <- RunMCA(df) # 该软件的主要分析函数，将细胞和基因同时降维到一个空间，离细胞近的基因被定义为细胞的特征基因——令人窒息的操作

# 下载参考基因集
panglao <- read_tsv("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz")

# filter based on organ
organs <- c("Liver", "Immune system", "Vasculature")
panglao_liver <- panglao %>% filter(organ == organs)
panglao_liver <- panglao_liver %>% filter(str_detect(species, "Mm"))
m_gene <- homologene(panglao_liver$`official gene symbol`,
    inTax = 9606, outTax = 10090
)
panglao_liver_m_gene <- merge(panglao_liver, m_gene,
    by.x = "official gene symbol", by.y = "9606"
)

# 下面两步会得到一个列表，列表的每一个元素是由基因集构成的字符向量，且每个元素被命名为对应的细胞类型的名字
panglao_liver_m_gene <- panglao_liver_m_gene %>%
    group_by(`cell type`) %>%
    summarise(geneset = list(`10090`))
# summarise( )几乎适用于任何聚合函数，并允许进行额外的算术运算：
# n() - 给出观察次数
# n_distinct(var) - 给出唯一值的数量 var
# sum(var), max(var), min(var), ...
# mean(var), median(var), sd(var), IQR(var)

liver_gs <- setNames(
    panglao_liver_m_gene$geneset,
    panglao_liver_m_gene$`cell type`
)

liver_gs <- liver_gs[sapply(liver_gs, length) >= 5]

# Performing per-cell hypergeometric tests against the gene signature collection
HGT_liver_gs <- RunCellHGT(df, pathways = liver_gs, dims = 1:50, n.features = 200)
liver_gs_prediction <- rownames(HGT_liver_gs)[apply(HGT_liver_gs, 2, which.max)]
# 矩阵的每一列依次：判断是否是最大值，得到一列布尔值，结合矩阵的行名会返回行名中的一个元素（也就是最有可能的细胞类型）。
# 所有列运行完了之后会得到所有细胞最可能的注释
liver_gs_prediction_signif <- ifelse(apply(HGT_liver_gs, 2, max) > 2,
    yes = liver_gs_prediction, "unassigned"
)
# 如果`-log10 corrected p-value`的值小于等于2，则认为不显著，注释更正为"unassigned"
df$liver_gs_prediction <- liver_gs_prediction_signif
df$treat <- df$orig.ident
Idents(df) = df$liver_gs_prediction
DimPlot(df, reduction = "umap", group.by = "liver_gs_prediction", pt.size = 0.5, label = TRUE, repel = TRUE)

levels(Idents(df)) # 查看细胞亚群，与上述结果一致
## 根据levels(Idents(df)) 顺序重新赋予对应的 细胞亚群名称，顺序不能乱
new.cluster.ids <- c(
    "Endothelial cells", "Endothelial cells", "Dendritic cells", "Kupffer cells", "Macrophages", "B cells", "B cells memory", "B cells naive", "T memory cells",
    "NK cells", "T cells", "Monocytes", "Plasma cells", "Pericytes", "Hepatocytes", "Neutrophils"
)
names(new.cluster.ids) <- levels(df)
df <- RenameIdents(df, new.cluster.ids)
levels(df) # 查看是否已改名
DimPlot(df,
    reduction = "umap",
    label = TRUE, pt.size = 0.5
)
ggsave("final_obj_heptocyte.png")
saveRDS(df, "final_obj_heptocyteV2.rds")
# deg_list <- FindAllMarkers(
#     df,
#     assay = "RNA",
#     logfc.threshold = 0.7,
#     test.use = "wilcox",
#     group.by = "liver_gs_prediction",
#     min.pct = 0.1,
#     min.diff.pct = -Inf,
#     verbose = TRUE,
#     only.pos = TRUE,
# )
# write.csv(deg_list, "deg_list.csv")

only_hep <- df[, Idents(df) %in% c("Hepatocytes")]
Idents(only_hep) <- only_hep$treat
levels(only_hep)

new_id <- c("chow", "chow", "nash", "nash", "nash")
names(new_id) <- levels(only_hep)
only_hep <- RenameIdents(only_hep, new_id)
hep_CHow_vs_Nash_Deg <- FindMarkers(only_hep,
    ident.1 = "chow",
    ident.2 = "nash"
)
write.csv(hep_CHow_vs_Nash_Deg, "hep_CHow_vs_Nash_Deg.csv")
# saveRDS(df, "only_hep_obj_V2.rds")
