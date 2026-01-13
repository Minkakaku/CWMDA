library(ShortRead)

count_fastq <- function(fq_file) {
  fq <- readFastq(fq_file)
  seqs <- as.character(sread(fq))
  as.data.frame(table(seqs), stringsAsFactors = FALSE)
}
count_tables <- lapply(fq, count_fastq)
names(count_tables) <- basename(fq)


all_seqs <- Reduce(union, lapply(count_tables, function(x) x$seqs))

counts_mat <- matrix(
  0,
  nrow = length(all_seqs),
  ncol = length(count_tables),
  dimnames = list(all_seqs, names(count_tables))
)

for (i in seq_along(count_tables)) {
  ct <- count_tables[[i]]
  counts_mat[ct$seqs, i] <- ct$Freq
}

counts_mat[1:10, 1:5]


count_list <- list(
  counts = as.data.frame(counts_mat)
)
