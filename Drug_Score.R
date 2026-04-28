rm(list = ls())
options(stringsAsFactors = FALSE)

## =========================
## 1) Install & load packages
## =========================
pkgs <- c(
  "data.table",
  "fgsea",
  "pheatmap",
  "ggplot2",
  "dplyr",
  "stringr",
  "tibble"
)
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(pheatmap)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(tibble)
})

set.seed(123)

## =========================
## 2) IO helper: fread smart
## =========================
fread_smart <- function(file, sep = NULL, header = TRUE, ...) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  ext <- tolower(tools::file_ext(file))
  if (ext %in% c("xls", "xlsx")) {
    message(
      "[fread_smart] Detected .",
      ext,
      " - will try fread as text. If it fails, this is a real Excel binary file."
    )
  }

  dt <- tryCatch(
    fread(file, sep = sep, header = header, data.table = TRUE, ...),
    error = function(e) {
      stop(
        "fread failed: ",
        conditionMessage(e),
        "\nIf this is a real .xls/.xlsx, use readxl::read_excel instead."
      )
    }
  )
  return(dt)
}

## =========================
## 3) Read GMT and keep 4 pathways
## =========================
read_gmt_fgsea <- function(gmt_file) {
  if (!file.exists(gmt_file)) {
    stop("GMT not found: ", gmt_file)
  }
  gmt <- fread_smart(gmt_file, sep = "\t", header = FALSE)
  # GMT format: pathway \t description \t gene1 \t gene2 ...
  pathways <- split(gmt[, -c(1, 2)], gmt[[1]])
  pathways <- lapply(pathways, function(x) {
    unique(na.omit(as.character(unlist(x))))
  })
  return(pathways)
}

## 你需要改：gmt 文件路径 + 你关心的 4 条通路名字（必须与 gmt 里的 pathway 名完全一致）
gmt_file <- "path/to/your.gmt"
pathways_all <- read_gmt_fgsea(gmt_file)

pathways_keep <- c(
  "ICD",
  "Antigen",
  "Glycolysis",
  "Lipogenesis"
)

missing_pw <- setdiff(pathways_keep, names(pathways_all))
if (length(missing_pw) > 0) {
  stop(
    "These pathways are not found in GMT: ",
    paste(missing_pw, collapse = ", ")
  )
}

pathways_4 <- pathways_all[pathways_keep]

## =========================
## 4) One-drug GSEA scoring
## =========================
run_fgsea_one <- function(
  dep_dt,
  pathways,
  gene_col = "gene",
  stat_col = "stat",
  minSize = 10,
  maxSize = 500,
  nperm = 10000
) {
  if (!all(c(gene_col, stat_col) %in% colnames(dep_dt))) {
    stop("dep_dt must contain columns: ", gene_col, " and ", stat_col)
  }
  x <- dep_dt[, .(
    gene = as.character(get(gene_col)),
    stat = as.numeric(get(stat_col))
  )]
  x <- x[!is.na(gene) & gene != "" & !is.na(stat)]
  x <- x[!duplicated(gene)] 

  stats <- x$stat
  names(stats) <- x$gene
  stats <- sort(stats, decreasing = TRUE)

  fg <- fgsea(
    pathways = pathways,
    stats = stats,
    minSize = minSize,
    maxSize = maxSize,
    nperm = nperm
  )

  # 输出 4 条通路的 NES（也可换成 ES）
  out <- fg[, .(pathway, NES, ES, pval, padj, size)]
  out <- out[match(names(pathways), pathway)] # 按 pathways 输入顺序对齐
  return(out)
}

## =========================
## 5) Batch: 30 drugs -> score matrix
## =========================
# 需要准备：30 个药物差异表文件路径（一个药物一个文件）
# drug_names 与 dep_files 一一对应
dep_files <- c(
  "drug01_dep.tsv",
  "drug02_dep.tsv"
  # ...补齐到 30 个
)

drug_names <- paste0("Drug", sprintf("%02d", seq_along(dep_files)))
if (length(drug_names) != length(dep_files)) {
  stop("drug_names length must match dep_files length.")
}

# 你需要根据你的差异表实际列名改：gene_col / stat_col
gene_col <- "gene"
stat_col <- "stat"

gsea_list <- vector("list", length(dep_files))
names(gsea_list) <- drug_names

for (i in seq_along(dep_files)) {
  dt <- fread_smart(dep_files[i])
  gsea_res <- run_fgsea_one(
    dep_dt = dt,
    pathways = pathways_4,
    gene_col = gene_col,
    stat_col = stat_col,
    nperm = 10000
  )
  gsea_list[[i]] <- gsea_res
}

# 组装 Drug x Pathway 矩阵（NES）
score_mat <- do.call(
  rbind,
  lapply(names(gsea_list), function(dn) {
    x <- gsea_list[[dn]]
    data.frame(
      Drug = dn,
      Pathway = x$pathway,
      NES = x$NES,
      stringsAsFactors = FALSE
    )
  })
)

score_wide <- tidyr::pivot_wider(
  score_mat,
  names_from = Pathway,
  values_from = NES
)
score_df <- as.data.frame(score_wide)
rownames(score_df) <- score_df$Drug
score_df$Drug <- NULL

## =========================
## 6) Heatmap: Drug x Pathway
## =========================
outdir <- "drug_scoring_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

pdf(file.path(outdir, "01_program_scores_heatmap.pdf"), width = 4.5, height = 8)
pheatmap(
  mat = as.matrix(score_df),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  border_color = NA,
  fontsize_row = 7,
  fontsize_col = 9,
  main = "Program scores (NES)"
)
dev.off()

fwrite(
  data.table::as.data.table(cbind(Drug = rownames(score_df), score_df)),
  file.path(outdir, "01_program_scores_matrix.csv")
)

## =========================
## 7) Drug single score + pairwise synergy (A+B)
## =========================
combine_mode <- "sum" # "sum" | "mean" | "maxabs"
calc_drug_score <- function(v, mode = "sum") {
  v <- as.numeric(v)
  v <- v[!is.na(v)]
  if (length(v) == 0) {
    return(NA_real_)
  }
  if (mode == "sum") {
    return(sum(v))
  }
  if (mode == "mean") {
    return(mean(v))
  }
  if (mode == "maxabs") {
    return(v[which.max(abs(v))])
  }
  stop("Unknown combine_mode: ", mode)
}

drug_score <- apply(score_df, 1, calc_drug_score, mode = combine_mode)
names(drug_score) <- rownames(score_df)

# 两两组合（不含自身）
pair_df <- expand.grid(
  A = seq_along(drug_names),
  B = seq_along(drug_names)
) |>
  dplyr::filter(A < B) |>
  dplyr::mutate(
    DrugA = drug_names[A],
    DrugB = drug_names[B],
    scoreA = drug_score[DrugA],
    scoreB = drug_score[DrugB],
    synergy_raw = scoreA + scoreB
  ) |>
  dplyr::filter(!is.na(synergy_raw))

# 只画 Top N（稀疏图）
top_n <- 120
pair_top <- pair_df |>
  dplyr::arrange(dplyr::desc(synergy_raw)) |>
  dplyr::slice_head(n = top_n)

# 为了像 demo 一样 0-10 显示：做 min-max scaling
scale_0_10 <- function(x) {
  if (length(unique(x)) == 1) {
    return(rep(10, length(x)))
  }
  (x - min(x)) / (max(x) - min(x)) * 10
}
pair_top$synergy_0_10 <- scale_0_10(pair_top$synergy_raw)

# 生成 Top pairs 文本（前 6 个）
top_k_label <- 6
top_pairs_txt <- pair_top |>
  dplyr::slice_head(n = top_k_label) |>
  dplyr::mutate(
    line = paste0(
      dplyr::row_number(),
      ". ",
      DrugA,
      " + ",
      DrugB,
      " = ",
      sprintf("%.2f", synergy_0_10)
    )
  ) |>
  dplyr::pull(line) |>
  paste(collapse = "\n")

# 作图
p_pair <- ggplot(pair_top, aes(x = A - 1, y = B - 1)) + 
  geom_point(aes(size = synergy_0_10, color = synergy_0_10), alpha = 0.9) +
  scale_size_continuous(range = c(2, 7)) +
  labs(
    title = "Top synergistic pairs as a sparse interaction map",
    x = "Drug index (A)",
    y = "Drug index (B)",
    color = "Integrated synergy score",
    size = "Integrated synergy score"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 12),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  ) +
  annotate(
    geom = "label",
    x = max(pair_top$A - 1) * 0.75,
    y = min(pair_top$B - 1) + 2,
    label = paste0("Top pairs\n(rank: A + B = score)\n\n", top_pairs_txt),
    hjust = 0,
    size = 3
  )

pdf(file.path(outdir, "02_top_pairs_sparse_map.pdf"), width = 12, height = 5)
print(p_pair)
dev.off()

fwrite(pair_top, file.path(outdir, "02_top_pairs_table.csv"))