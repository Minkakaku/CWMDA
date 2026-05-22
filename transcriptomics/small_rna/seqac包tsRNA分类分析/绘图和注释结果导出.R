out_pdf <- "PAC_trna_Group_type_class.pdf"

# 如果之前有同名损坏文件，先删掉
if (file.exists(out_pdf)) file.remove(out_pdf)

# 用 tryCatch 确保设备一定被关闭
tryCatch({
  pdf(out_pdf, width = 9, height = 6, useDingbats = FALSE)
  
  PAC_trna(
    PAC = pac_trna,
    norm = "cpm",
    filter = NULL,
    join = TRUE,
    top = 15,
    log2fc = TRUE,
    pheno_target = list(
      "Group",
      c("NC", "TNBCexo", "Her2exo", "HRexo")
    ),
    anno_target_1 = list("type"),
    anno_target_2 = list("class")
  )
  
}, finally = {
  dev.off()
})

cowplot::plot_grid(trna_result$plots$Expression_Anno_1$Grand_means,
                   trna_result$plots$Log2FC_Anno_1,
                   trna_result$plots$Percent_bars$Grand_means,
                   nrow=1, ncol=3)
anno_result <- anno(pac_trna)



saveRDS(anno_result, file = "pac_trna_annotation.rds")
write.csv(
  anno_result,
  file = "pac_trna_annotation.csv",
  row.names = TRUE
)



counts_mat <- pac_trna@Counts
anno_result <- anno(pac_trna)

trf_anno <- anno_result[grepl("tRF", anno_result$class), ]

trf_seqs <- rownames(trf_anno)

length(trf_seqs)
counts_mat <- pac_trna@Counts

counts_trf <- counts_mat[trf_seqs, , drop = FALSE]
dim(counts_trf)
counts_trf[1:5, 1:5]
trf_anno$tRF_type <- NA_character_

trf_anno$tRF_type[grepl("5'-tRF", trf_anno$class)] <- "tRF-5"
trf_anno$tRF_type[grepl("3'-tRF", trf_anno$class)] <- "tRF-3"
trf_anno$tRF_type[grepl("i'-tRF", trf_anno$class)] <- "i-tRF"
table(trf_anno$tRF_type)
counts_trf
trf_type_sample_mat <- rowsum(
  counts_trf,
  group = trf_anno$tRF_type
)
saveRDS(trf_anno,   file = "tRF_annotation.rds")

write.csv(
  counts_trf,
  file = "tRF_sequence_by_sample_counts.csv",
  row.names = TRUE
)



# 读 counts
counts_df <- read.csv(
  "tRF_sequence_by_sample_counts.csv",
  row.names = 1,
  check.names = FALSE
)

# 读 annotation
anno_df <- read.csv(
  "pac_trna_annotation.csv",
  row.names = 1,
  stringsAsFactors = FALSE
)
anno_trf <- anno_df[grepl("tRF", anno_df$class), , drop = FALSE]
# 取交集（防止顺序/多余行问题）
common_seqs <- intersect(rownames(counts_df), rownames(anno_trf))

length(common_seqs)   # sanity check
counts_trf_anno <- cbind(
  sequence = common_seqs,
  anno_trf[common_seqs, , drop = FALSE],
  counts_df[common_seqs, , drop = FALSE]
)
rownames(counts_trf_anno) <- counts_trf_anno$sequence
write.csv(
  counts_trf_anno,
  file = "tRF_sequence_by_sample_counts_annotated.csv",
  row.names = FALSE
)
