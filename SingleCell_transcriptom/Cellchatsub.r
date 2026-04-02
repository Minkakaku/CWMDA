rm(list=ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
  library(CellChat)
  library(ComplexHeatmap)
  library(RColorBrewer)
})

set.seed(123)

epi <- readRDS("03_Subcluster\\Epithelial\\subcluster_annotated.rds")
caf = readRDS("03_Subcluster\\Fibroblast\\subcluster_annotated.rds")  
seu <- merge(epi, caf)
# 定义需要剔除的 VSMC 名称（可按你自己的 annotation 微调）
# vsmc_pattern <- c("Epithelial")

# vsmc_cells <- grep(
#   paste(vsmc_pattern, collapse = "|"),
#   seu$sub_cell_type,
#   ignore.case = TRUE
# )

# cat("Remove Epithelial cells:", length(vsmc_cells), "\n")

# # 真正剔除
# seu <- subset(
#   seu,
#   cells = colnames(seu)[-vsmc_cells]
# )

# # 清理 factor，避免残留 level
# seu$sub_cell_type <- droplevels(factor(seu$sub_cell_type))


seu <- NormalizeData(seu, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)

outdir <- "06_epi_caf_cellchat"
dir.create(outdir, showWarnings = FALSE)
DefaultAssay(seu) <- "RNA"

# 看一眼 sub_cell_type 都有什么（第一次跑建议保留）
print(table(seu$cell_type))

## group
seu$group <- factor(seu$group, levels = c("CON", "T1","T2"))
print(table(seu$group, useNA = "ifany"))

Idents(seu) <- "group"

# CON    <- subset(seu, idents = "CON")
T1 <- subset(seu, idents = "T1")
T2<- subset(seu, idents = "T2")

## Seurat v5 兼容
# try(CON<- JoinLayers(CON),    silent = TRUE)
try(T1 <- JoinLayers(T1), silent = TRUE)
try(T2 <- JoinLayers(T2), silent = TRUE)
## =========================
## 3) CellChat function
## =========================
run_cellchat_one <- function(
  seu = T1,
  species = c("mouse", "human"),
  group_col = "sub_cell_type",
  min_cells = 10,
  out_rds = NULL
) {
  species <- match.arg(species)

  DefaultAssay(seu) <- "RNA"
  data.input <- GetAssayData(seu, layer = "data")

  meta <- seu@meta.data[, c("group", group_col)]
  colnames(meta) <- c("Group", "CellTypes")

  cellchat <- createCellChat(
    object = data.input,
    meta = meta,
    group.by = "CellTypes"
  )

  ## 正确的人/鼠数据库
  if (species == "mouse") {
    cellchat@DB <- CellChatDB.mouse
  } else {
    cellchat@DB <- CellChatDB.human
  }
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  if (species == "mouse") {
    cellchat <- smoothData(cellchat, adj = PPI.mouse)
  } else {
    cellchat <- smoothData(cellchat, adj = PPI.human)
  }
  cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  cellchat <- computeNetSimilarity(cellchat, type = "functional")

  saveRDS(cellchat, out_rds)
  cellchat
}
## =========================
## 4) Run
## =========================

# cco.CON <- run_cellchat_one(
#   CON,
#   species = "mouse",
#   out_rds = file.path(outdir, "cco_CON.rds")
# )
cco.T1 <- run_cellchat_one(
  T1,
  species = "mouse",
  out_rds = file.path(outdir, "cco_T1.rds")
)

cco.T2 <- run_cellchat_one(
  T2,
  species = "mouse",
  out_rds = file.path(outdir, "cco_T2.rds")
)

# cco.CON = readRDS("05_cellchat\\cco_CON.rds")
cco.T1 = readRDS("06_epi_caf_cellchat\\cco_T1.rds")
cco.T2 = readRDS("06_epi_caf_cellchat\\cco_T2.rds")
cco.list <- list(T1 = cco.T1,T2 = cco.T2)

cellchat_merged <- mergeCellChat(
  cco.list,
  add.names = names(cco.list),
  cell.prefix = TRUE
)

joint_levels <- levels(cellchat_merged@idents$joint)
print(joint_levels)

names(cco.list)

plotdir <- file.path(outdir, "plots")
tbldir  <- file.path(outdir, "tables")
dir.create(plotdir, showWarnings = FALSE, recursive = TRUE)
dir.create(tbldir,  showWarnings = FALSE, recursive = TRUE)

object.list <- cco.list   # list(T1=..., T2=...)
cellchat <- cellchat_merged

## optional: update (safe)
try(cellchat <- updateCellChat(cellchat), silent = TRUE)

## dataset names / indices
dataset.names <- names(object.list)
stopifnot(length(dataset.names) == 2)
cmp <- c(1, 2)  # (T1 vs T2) in merged object

## -------------------------
## 5.1) Compare total interactions (count / weight)
## -------------------------
pdf(file.path(plotdir, "01_compareInteractions_total.pdf"), width = 8, height = 4)
gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = cmp)
gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = cmp, measure = "weight")
print(gg1 + gg2)
dev.off()

## -------------------------
## 5.2) Differential interactions (circle + heatmap)
## -------------------------
pdf(file.path(plotdir, "02_diffInteraction_circle_count_weight.pdf"), width = 12, height = 6)
par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE)                    # count diff
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")# weight diff
dev.off()

pdf(file.path(plotdir, "03_diffInteraction_heatmap_count_weight.pdf"), width = 10, height = 5)
gg1 <- netVisual_heatmap(cellchat)                      # count diff heatmap
gg2 <- netVisual_heatmap(cellchat, measure = "weight")  # weight diff heatmap
print(gg1 + gg2)
dev.off()

## -------------------------
## 5.3) Circle plots per dataset with unified edge scaling
## -------------------------
weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))
pdf(file.path(plotdir, "04_circle_eachDataset_count_unifiedScale.pdf"), width = 12, height = 6)
par(mfrow = c(1,2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(
    object.list[[i]]@net$count,
    weight.scale   = TRUE,
    label.edge     = FALSE,
    edge.weight.max= weight.max[2],
    edge.width.max = 12,
    title.name     = paste0("Number of interactions - ", names(object.list)[i])
  )
}
dev.off()

weight.max.w <- getMaxWeight(object.list, attribute = c("idents", "weight"))
pdf(file.path(plotdir, "05_circle_eachDataset_weight_unifiedScale.pdf"), width = 12, height = 6)
par(mfrow = c(1,2), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(
    object.list[[i]]@net$weight,
    weight.scale   = TRUE,
    label.edge     = FALSE,
    edge.weight.max= weight.max.w[2],
    edge.width.max = 12,
    title.name     = paste0("Interaction weights - ", names(object.list)[i])
  )
}
dev.off()

## 5.4.0) 你给的顺序
celltype_order <- joint_levels

## 5.4.1) 固定颜色（按 celltype_order）
color.use <- rainbow(length(celltype_order))
names(color.use) <- celltype_order

## 5.4.2) 统一 edge 上限（跨 dataset，使用 weight）
weight.max.w <- getMaxWeight(object.list, attribute = c("idents", "weight"))
edge_max_weight <- weight.max.w[2]  # 用于 edge.weight.max 统一尺度

## 5.4.3) 每个 dataset 循环输出（T1 / T2）
for (k in seq_along(object.list)) {

  obj      <- object.list[[k]]
  obj_name <- names(object.list)[k]

  ## groupSize：节点大小（按 idents 统计）
  groupSize <- as.numeric(table(obj@idents))
  names(groupSize) <- names(table(obj@idents))

  ## 取 weight 矩阵并按 celltype_order 排序
  mat <- as.matrix(obj@net$weight)

  ## 防御：确保 order 都在矩阵里
  miss_row <- setdiff(celltype_order, rownames(mat))
  miss_col <- setdiff(celltype_order, colnames(mat))
  if (length(miss_row) > 0 || length(miss_col) > 0) {
    stop(
      "对象 ", obj_name, " 的 net$weight 缺少这些 cell type：\n",
      "missing rows: ", paste(miss_row, collapse = ", "), "\n",
      "missing cols: ", paste(miss_col, collapse = ", ")
    )
  }

  mat <- mat[celltype_order, celltype_order, drop = FALSE]

  ## groupSize 也按同样顺序排列（没有的会 NA，这里不会，因为上面已经检查）
  groupSize <- groupSize[celltype_order]
  groupSize[is.na(groupSize)] <- 0

  ## 输出 PDF：每个 dataset 一份
  pdf(file.path(plotdir, paste0("06_circle_senderOnly_eachCelltype_", obj_name, "_weight_unifiedScale.pdf")),
      width = 14, height = 18)

  ## 如果显示不全：调 mfrow / mar
  nplot <- nrow(mat)
  ncolp <- 4
  nrowp <- ceiling(nplot / ncolp)
  par(mfrow = c(nrowp, ncolp), xpd = TRUE, mar = c(1, 1, 2, 1))

  for (i in seq_len(nrow(mat))) {

    mat2 <- matrix(
      0,
      nrow = nrow(mat), ncol = ncol(mat),
      dimnames = dimnames(mat)
    )
    ## 只保留 sender = 第 i 行
    mat2[i, ] <- mat[i, ]

    netVisual_circle(
      mat2,
      vertex.weight   = groupSize,
      weight.scale    = TRUE,
      label.edge      = FALSE,
      arrow.size      = 0.05,
      arrow.width     = 1,
      edge.weight.max = edge_max_weight,  # 跨 dataset 统一尺度
      edge.width.max  = 12,
      title.name      = paste0(obj_name, " | Sender: ", rownames(mat)[i]),
      color.use       = color.use
    )
  }

  dev.off()
}

# -------------------------
# count
# -------------------------
weight.max.c <- getMaxWeight(object.list, attribute = c("idents", "count"))
edge_max_count <- weight.max.c[2]
for (k in seq_along(object.list)) {
  obj <- object.list[[k]]; obj_name <- names(object.list)[k]
  groupSize <- as.numeric(table(obj@idents)); names(groupSize) <- names(table(obj@idents))
  mat <- as.matrix(obj@net$count)
  miss_row <- setdiff(celltype_order, rownames(mat)); miss_col <- setdiff(celltype_order, colnames(mat))
  if (length(miss_row) > 0 || length(miss_col) > 0) stop("对象 ", obj_name, " 的 net$count 缺少 cell type")
  mat <- mat[celltype_order, celltype_order, drop = FALSE]
  groupSize <- groupSize[celltype_order]; groupSize[is.na(groupSize)] <- 0
  pdf(file.path(plotdir, paste0("07_circle_senderOnly_eachCelltype_", obj_name, "_count_unifiedScale.pdf")),
      width = 14, height = 18)
  nplot <- nrow(mat); ncolp <- 4; nrowp <- ceiling(nplot / ncolp)
  par(mfrow = c(nrowp, ncolp), xpd = TRUE, mar = c(1, 1, 2, 1))
  for (i in seq_len(nrow(mat))) {
    mat2 <- matrix(0, nrow=nrow(mat), ncol=ncol(mat), dimnames=dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight=groupSize, weight.scale=TRUE, label.edge=FALSE,
                     arrow.size=0.05, arrow.width=1, edge.weight.max=edge_max_count,
                     edge.width.max=12, title.name=paste0(obj_name," | Sender: ",rownames(mat)[i]),
                     color.use=color.use)
  }
  dev.off()
}


## -------------------------
## 5.4) Compare major sources/targets (scatter) for each dataset
## -------------------------
num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})
weight.MinMax <- c(min(num.link), max(num.link))

pdf(file.path(plotdir, "08_signalingRole_scatter_eachDataset.pdf"), width = 10, height = 5)
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(
    object.list[[i]],
    title = names(object.list)[i],
    weight.MinMax = weight.MinMax
  )
}
print(patchwork::wrap_plots(plots = gg))
dev.off()

## -------------------------
## 5.6) Functional / structural similarity (pairwise) + embedding + clustering
## -------------------------
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
pdf(file.path(plotdir, "10_embedding_functional.pdf"), width = 10, height = 8)
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 1.5)
dev.off()

pdf(file.path(plotdir, "11_rankSimilarity_functional.pdf"), width = 7, height = 10)
rankSimilarity(cellchat, type = "functional")
dev.off()

cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
pdf(file.path(plotdir, "12_embedding_structural.pdf"), width = 10, height = 8)
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 1.5)
dev.off()

pdf(file.path(plotdir, "13_rankSimilarity_structural.pdf"), width = 7, height = 10)
rankSimilarity(cellchat, type = "structural")
dev.off()

saveRDS(cellchat, file.path(outdir, "cellchat_merged_T1_T2_afterComparisonAnalysis.rds"))

## -------------------------
## 5.7) Compare information flow (rankNet)
## -------------------------
pdf(file.path(plotdir, "14_rankNet_comparison.pdf"), width = 12, height = 12)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = TRUE,  do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE)
print(gg1 + gg2)
dev.off()

## -------------------------
## 5.8) Outgoing / incoming / all signaling heatmaps (per dataset; pathway union)
## -------------------------
i <- 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

pdf(
  file.path(plotdir, "15_signalingRole_heatmap_outgoing.pdf"),
  width = 12,
  height = 25
)
ht1 <- netAnalysis_signalingRole_heatmap(
  object.list[[i]],
  pattern = "outgoing",
  signaling = pathway.union,
  title = names(object.list)[i],
  width = 5,
  height = 25
)
ht2 <- netAnalysis_signalingRole_heatmap(
  object.list[[i + 1]],
  pattern = "outgoing",
  signaling = pathway.union,
  title = names(object.list)[i + 1],
  width = 5,
  height = 25
)
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

pdf(
  file.path(plotdir, "16_signalingRole_heatmap_incoming.pdf"),
  width = 12,
  height = 25
)
ht1 <- netAnalysis_signalingRole_heatmap(
  object.list[[i]],
  pattern = "incoming",
  signaling = pathway.union,
  title = names(object.list)[i],
  width = 5,
  height = 25,
  color.heatmap = "GnBu"
)
ht2 <- netAnalysis_signalingRole_heatmap(
  object.list[[i + 1]],
  pattern = "incoming",
  signaling = pathway.union,
  title = names(object.list)[i + 1],
  width = 5,
  height = 25,
  color.heatmap = "GnBu"
)
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

pdf(
  file.path(plotdir, "17_signalingRole_heatmap_all.pdf"),
  width = 12,
  height = 25
)
ht1 <- netAnalysis_signalingRole_heatmap(
  object.list[[i]],
  pattern = "all",
  signaling = pathway.union,
  title = names(object.list)[i],
  width = 5,
  height = 25,
  color.heatmap = "OrRd"
)
ht2 <- netAnalysis_signalingRole_heatmap(
  object.list[[i + 1]],
  pattern = "all",
  signaling = pathway.union,
  title = names(object.list)[i + 1],
  width = 5,
  height = 25,
  color.heatmap = "OrRd"
)
ComplexHeatmap::draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

## -------------------------
## 5.9) Bubble plot: compare L-R probabilities between T1 and T2 (src -> tgt)
##      (A) all L-R in these groups
##      (B) increased vs decreased (by comm. prob)
## -------------------------

pdf(file.path(plotdir, "18_bubble_compare_allLR_src_to_tgt.pdf"), width = 20, height = 50)
p <- netVisual_bubble(
  cellchat,
  comparison  = cmp,
  angle.x     = 45
)
print(p)
dev.off()

pdf(file.path(plotdir, "19_bubble_increased_decreased_byProb.pdf"), width = 50, height = 50)
gg.increase <- netVisual_bubble(
  cellchat,
  comparison     = cmp,
  max.dataset    = 2,  # increased in T2 vs T1
  title.name     = paste0("Increased signaling in ", dataset.names[2]),
  angle.x        = 45,
  remove.isolate = TRUE
)
gg.decrease <- netVisual_bubble(
  cellchat,
  comparison     = cmp,
  max.dataset    = 1,  # decreased in T2 vs T1 (i.e., higher in T1)
  title.name     = paste0("Decreased signaling in ", dataset.names[2]),
  angle.x        = 45,
  remove.isolate = TRUE
)
print(gg.increase + gg.decrease)
dev.off()

## 导出 bubble data（重要：后续你可以筛 L-R）
try(write.csv(gg.increase$data, file.path(tbldir, "bubble_increased_in_T2.csv"), row.names = FALSE), silent = TRUE)
try(write.csv(gg.decrease$data, file.path(tbldir, "bubble_decreased_in_T2.csv"), row.names = FALSE), silent = TRUE)

## -------------------------
## 5.10) DEG-based up/down L-R (within CellChat framework)
##       NOTE: 这里 group.dataset="datasets" 依赖 cellchat@meta$datasets
## -------------------------
pos.dataset   <- dataset.names[2]          # T2 as "positive" vs T1
features.name <- pos.dataset

cellchat <- identifyOverExpressedGenes(
  cellchat,
  group.dataset = "datasets",
  pos.dataset   = pos.dataset,
  features.name = features.name,
  only.pos      = FALSE,
  thresh.pc     = 0.1,
  thresh.fc     = 0.1,
  thresh.p      = 1
)

net.deg <- netMappingDEG(cellchat, features.name = features.name)

## up: ligands higher in T2 (logFC threshold adjustable)
net.up <- subsetCommunication(
  cellchat,
  net      = net.deg,
  datasets = pos.dataset,
  ligand.logFC   = 0.2,
  receptor.logFC = NULL
)

## down: ligands+receptors higher in T1 (i.e., negative in T2)
net.down <- subsetCommunication(
  cellchat,
  net      = net.deg,
  datasets = dataset.names[1],
  ligand.logFC   = -0.1,
  receptor.logFC = -0.1
)

write.csv(net.up,   file.path(tbldir, "DEG_based_up_LR_in_T2.csv"),   row.names = FALSE)
write.csv(net.down, file.path(tbldir, "DEG_based_down_LR_in_T2.csv"), row.names = FALSE)

gene.up   <- extractGeneSubsetFromPair(net.up,   cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
write.csv(gene.up,   file.path(tbldir, "DEG_based_up_genes_in_T2.csv"),   row.names = FALSE)
write.csv(gene.down, file.path(tbldir, "DEG_based_down_genes_in_T2.csv"), row.names = FALSE)

## bubble: only DEG-based up/down L-R, restricted to your src->tgt
pairLR.use.up   <- net.up[,   "interaction_name", drop = FALSE]
pairLR.use.down <- net.down[, "interaction_name", drop = FALSE]

pdf(file.path(plotdir, "20_bubble_DEG_based_up_down_src_to_tgt.pdf"), width = 25, height = 10)
gg1 <- netVisual_bubble(
  cellchat,
  pairLR.use    = pairLR.use.up,
  comparison    = cmp,
  angle.x       = 90,
  remove.isolate= TRUE,
  title.name    = paste0("Up-regulated (DEG-based) in ", dataset.names[2])
)
gg2 <- netVisual_bubble(
  cellchat,
  pairLR.use    = pairLR.use.down,
  comparison    = cmp,
  angle.x       = 90,
  remove.isolate= TRUE,
  title.name    = paste0("Down-regulated (DEG-based) in ", dataset.names[2])
)
print(gg1 + gg2)
dev.off()

## chord: up/down (DEG-based) — 注意：chord 基于单个 dataset 的 CellChat object
##       这里用 object.list[[2]](T2) 画 up，用 object.list[[1]](T1) 画 down
pdf(file.path(plotdir, "21_chord_gene_DEG_based_up_down.pdf"), width = 14, height = 7)
par(mfrow = c(1,2), xpd = TRUE)
netVisual_chord_gene(
  object.list[[2]],
  slot.name   = "net",
  net         = net.up,
  lab.cex     = 0.8,
  small.gap   = 3.5,
  title.name  = paste0("Up-regulated (DEG-based) in ", dataset.names[2])
)
netVisual_chord_gene(
  object.list[[1]],
  slot.name   = "net",
  net         = net.down,
  lab.cex     = 0.8,
  small.gap   = 3.5,
  title.name  = paste0("Down-regulated (DEG-based) in ", dataset.names[2])
)
dev.off()


# 看特异性的LR
pathways <- c("EPHA", "EPHB", "CXCL", "IL6","TGFb","VEGF")

for (grp in names(cco.list)) {
  
  obj <- cco.list[[grp]]
  
  ## 当前对象里真实存在的通路
  available_pathways <- obj@netP$pathways
  
  for (pw in pathways) {
    
    ## 如果不存在，直接跳过
    if (!pw %in% available_pathways) {
      message("Skip: ", pw, " not found in ", grp)
      next
    }
    
    outfile <- file.path(
      outdir,
      "/plots/signalingRole",
      paste0(pw, "_pathways2_signalingRole_heatmap_", tolower(grp), ".pdf")
    )
    
    pdf(outfile, width = 7, height = 6)
    
    netAnalysis_signalingRole_network(
      obj,
      signaling = pw
    )
    
    dev.off()
    
    message("Saved: ", outfile)
  }
}



for (grp in names(cco.list)) {
  obj <- cco.list[[grp]]
  available_pathways <- obj@netP$pathways

  for (pw in pathways) {
    if (!pw %in% available_pathways) {
      message("Skip: ", pw, " not found in ", grp)
      next
    }

    outfile <- file.path(
      outdir,
      "/plots/aggregate",
      paste0("6NetVisual_", pw, "_circle_", tolower(grp), ".pdf")
    )

    pdf(outfile, width = 6, height = 5)
    par(xpd = TRUE)

    netVisual_aggregate(
      obj,
      signaling = pw,
      layout = "circle"
    )

    dev.off()
    message("Saved: ", outfile)
  }
}



## 批量输出 chord 图
for (grp in names(cco.list)) {
  obj <- cco.list[[grp]]
  available_pathways <- obj@netP$pathways

  for (pw in pathways) {
    if (!pw %in% available_pathways) {
      message("Skip: ", pw, " not found in ", grp)
      next
    }

    outfile <- file.path(
      outdir,
      "/plots/chord",
      paste0("7NetVisual_", pw, "_chord_", tolower(grp), ".pdf")
    )

    pdf(outfile, width = 10, height = 8)

    netVisual_aggregate(
      obj,
      signaling = pw,
      layout = "chord"
    )

    dev.off()

    message("Saved: ", outfile)
  }
}



## save final merged object
saveRDS(cellchat, file.path(outdir, "cellchat_merged_T1_T2_final.rds"))