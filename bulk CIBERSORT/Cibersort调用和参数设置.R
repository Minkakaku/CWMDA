# 加载CIBERSORT源码
source("C:/Users/Administrator/Desktop/WORK/tgca/20251028/CIBERSORT.R")  # 改成你的实际路径

# 运行CIBERSORT
# 参数解释：
#   sig_matrix   = LM22特征矩阵（固定的免疫签名）
#   mixture_file = 你的表达矩阵
#   perm         = 100 或 1000 等，用于置换检验p值
#   QN           = TRUE 对 microarray 推荐；对RNA-seq通常建议 FALSE
#                  RNA-seq + FPKM/TPM 的常见做法是 QN=FALSE，以避免强行分位数标准化破坏生物差异
#
results <- CIBERSORT(
  sig_matrix    = "C:/Users/Administrator/Desktop/WORK/tgca/20251028/LM22.txt",
  mixture_file  = "for_cibersort_clean.txt",
  perm          = 100,
  QN            = FALSE
)

# 查看结果
head(results)

