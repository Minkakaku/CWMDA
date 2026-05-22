library(Seurat)
library(SeuratDisk)
sce <- LoadH5Seurat("Microglia.h5seurat")

table(Idents(sce))

sce$celltype <- sce$sublabel
table(Idents(sce))

library(monocle)
sample_ann <- sce@meta.data
head(sample_ann)
# rownames(sample_ann)=sample_ann[,1]
gene_ann <- data.frame(
    gene_short_name = rownames(sce@assays$RNA),
    row.names = rownames(sce@assays$RNA)
)
head(gene_ann)

pd <- new("AnnotatedDataFrame",
    data = sample_ann
)
fd <- new("AnnotatedDataFrame",
    data = gene_ann
)
matrix <- as.sparse(sce@assays$RNA@counts)
ct[1:4, 1:4]

sc_cds <- newCellDataSet(
    matrix,
    phenoData = pd,
    featureData = fd,
    expressionFamily = negbinomial.size(),
    lowerDetectionLimit = 1
)
sc_cds
sc_cds <- estimateSizeFactors(sc_cds)
sc_cds <- estimateDispersions(sc_cds)

disp_table <- dispersionTable(sc_cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(sc_cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(sc_cds)
plot_pc_variance_explained(sc_cds, return_all = F) # norm_method='log'
# 其中 num_dim 参数选择基于上面的PCA图
sc_cds <- reduceDimension(sc_cds,
    max_components = 2, num_dim = 6,
    reduction_method = "tSNE", verbose = T
)
sc_cds <- clusterCells(sc_cds, num_clusters = 6)
plot_cell_clusters(sc_cds, 1, 2)
table(pData(sc_cds)$Cluster)
colnames(pData(sc_cds))

table(pData(sc_cds)$Cluster, pData(sc_cds)$new)
plot_cell_clusters(sc_cds, 1, 2)
library(monocle)
library(Seurat)
load(file = "input_cds.Rdata")

# 接下来很重要，到底是看哪个性状的轨迹
colnames(pData(sc_cds))
table(pData(sc_cds)$Cluster)
table(pData(sc_cds)$Cluster, pData(sc_cds)$celltype)
plot_cell_clusters(sc_cds, 1, 2)

## 我们这里并不能使用 monocle的分群
# 还是依据前面的 seurat分群, 其实取决于自己真实的生物学意图
pData(sc_cds)$Cluster <- pData(sc_cds)$celltype
table(pData(sc_cds)$Cluster)

Sys.time()
diff_test_res <- differentialGeneTest(sc_cds,
    fullModelFormulaStr = "~Cluster"
)
Sys.time()
# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval < 0.1)
sig_genes <- sig_genes[order(sig_genes$pval), ]
head(sig_genes[, c("gene_short_name", "pval", "qval")])
cg <- as.character(head(sig_genes$gene_short_name))
# 挑选差异最显著的基因可视化
plot_genes_jitter(sc_cds[cg, ],
    grouping = "Cluster",
    color_by = "Cluster",
    nrow = 3,
    ncol = NULL
)
cg2 <- as.character(tail(sig_genes$gene_short_name))
plot_genes_jitter(sc_cds[cg2, ],
    grouping = "Cluster",
    color_by = "Cluster",
    nrow = 3,
    ncol = NULL
)

# 前面是找差异基因，后面是做拟时序分析

# 第一步: 挑选合适的基因. 有多个方法，例如提供已知的基因集，这里选取统计学显著的差异基因列表
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
ordering_genes
sc_cds <- setOrderingFilter(sc_cds, ordering_genes)
plot_ordering_genes(sc_cds)
sc_cds <- reduceDimension(sc_cds,
    max_components = 2,
    method = "DDRTree"
)

# 第二步: 降维。降维的目的是为了更好的展示数据。函数里提供了很多种方法,
# 不同方法的最后展示的图都不太一样, 其中“DDRTree”是Monocle2使用的默认方法
sc_cds <- reduceDimension(sc_cds,
    max_components = 2,
    method = "DDRTree"
)
# 第三步: 对细胞进行排序
sc_cds <- orderCells(sc_cds)
# 最后两个可视化函数，对结果进行可视化
plot_cell_trajectory(sc_cds, color_by = "Cluster")
ggsave("monocle_cell_trajectory_for_seurat.pdf")

length(cg)
plot_genes_in_pseudotime(sc_cds[cg, ],
    color_by = "Cluster"
)
ggsave("monocle_plot_genes_in_pseudotime_for_seurat.pdf")

# https://davetang.org/muse/2017/10/01/getting-started-monocle/

my_sc_cds_subset <- sc_cds
# pseudotime is now a column in the phenotypic data as well as the cell state
head(pData(my_sc_cds_subset))
# 这个differentialGeneTest会比较耗费时间
my_pseudotime_de <- differentialGeneTest(my_sc_cds_subset,
    fullModelFormulaStr = "~sm.ns(Pseudotime)",
    cores = 1
)
# 不知道为什么无法开启并行计算了

head(my_pseudotime_de)
save(my_sc_cds_subset, my_pseudotime_de,
    file = "output_of_phe2_monocle.Rdata"
)
