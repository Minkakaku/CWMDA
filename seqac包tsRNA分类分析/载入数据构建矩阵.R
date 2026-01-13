library(seqpac)
setwd("C:/Users/Administrator/Desktop/WORK/蔡总20260106/data")
fq <- list.files(pattern = "fastq", all.files = FALSE, full.names = TRUE)
fq


count_list <- make_counts(input=fq, trimming="seqpac", threads=16, parse="default_neb")

# 查看count_list的结构
str(count_list, max.level = 1)