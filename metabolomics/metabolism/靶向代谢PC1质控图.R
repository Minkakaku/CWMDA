library(readr)
library(dplyr)
library(ggplot2)

# 读取 PCA 结果（第一列为样本名）
pca_df <- read.csv(
  "PCA_Score.csv",
  row.names = 1,
  check.names = FALSE
)

# 转为数据框并显式保留 Sample
pca_df <- pca_df %>%
  tibble::rownames_to_column("Sample")

pca_df <- pca_df %>%
  mutate(
    Type = ifelse(grepl("^QC", Sample), "QC", "Sample")
  )


pc1_mean <- mean(pca_df$PC1)
pc1_sd   <- sd(pca_df$PC1)

pc1_mean
pc1_sd



# -------------------------------------------------
# 1. 构造一个带“左侧留白”的 X 轴
# -------------------------------------------------
pca_df <- pca_df %>%
  mutate(
    Sample_plot = factor(
      Sample,
      levels = c(" ", levels(Sample))  # 人为加一个空位
    )
  )
# -------------------------------------------------
# 2. 画图（完全示例对齐版）
# -------------------------------------------------
p <- ggplot(pca_df, aes(x = Sample_plot, y = PC1)) +
  
  # 样本点
  geom_point(
    aes(color = Type),
    size = 3
  ) +
  
  # 均值线
  geom_hline(
    yintercept = pc1_mean,
    color = "grey40",
    linewidth = 0.8
  ) +
  
  # ±2SD
  geom_hline(
    yintercept = c(pc1_mean + 2 * pc1_sd,
                   pc1_mean - 2 * pc1_sd),
    linetype = "dashed",
    color = "#E69F00",
    linewidth = 0.9
  ) +
  
  # ±3SD
  geom_hline(
    yintercept = c(pc1_mean + 3 * pc1_sd,
                   pc1_mean - 3 * pc1_sd),
    linetype = "dashed",
    color = "#D55E00",
    linewidth = 1
  ) +
  
  # -----------------------------
# SD 文本：画在“左侧留白区”
# -----------------------------
annotate("text",
         x = 1,
         y = pc1_mean + 3 * pc1_sd,
         label = "+3SD",
         color = "#D55E00",
         hjust = 1, vjust = -0.3, size = 4) +
  
  annotate("text",
           x = 1,
           y = pc1_mean + 2 * pc1_sd,
           label = "+2SD",
           color = "#E69F00",
           hjust = 1, vjust = -0.3, size = 4) +
  
  annotate("text",
           x = 1,
           y = pc1_mean - 2 * pc1_sd,
           label = "-2SD",
           color = "#E69F00",
           hjust = 1, vjust = 1.2, size = 4) +
  
  annotate("text",
           x = 1,
           y = pc1_mean - 3 * pc1_sd,
           label = "-3SD",
           color = "#D55E00",
           hjust = 1, vjust = 1.2, size = 4) +
  
  # 颜色
  scale_color_manual(
    values = c(
      "QC" = "#F0C808",
      "Sample" = "black"
    )
  ) +
  
  # 标签
  labs(
    title = "PC1 scores / Outlier Detection",
    y = "PC1 scores / Outlier Detection",
    x = NULL,
    color = NULL
  ) +
  
  # 关键：允许画到面板外
  coord_cartesian(clip = "off") +
  
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    
    # 给左侧真正留出空间
    plot.margin = margin(t = 10, r = 10, b = 10, l = 40)
  )

p

