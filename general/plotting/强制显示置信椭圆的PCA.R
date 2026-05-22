library(ellipse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(ggExtra)     # 新增：用于添加边际分布

expr <- read.csv("protein.csv")
group_info <- read.csv("group.csv")
rownames(expr) <- expr$Protein.accession
expr <- expr[ , -1]

common_samples <- intersect(colnames(expr), group_info$sample)
expr <- expr[ , common_samples]
group_info <- group_info[group_info$sample %in% common_samples, ]
group_info <- group_info[match(colnames(expr), group_info$sample), ]
expr <- expr[complete.cases(expr), ]

# ===== 3. PCA 计算 =====
pca_res <- prcomp(t(expr), scale. = TRUE)
pc_df <- as.data.frame(pca_res$x[, 1:2])
pc_df$group <- group_info$group
pc_df$sample <- rownames(pc_df)

expl <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

# ===== 4. 设置颜色 =====
color_map <- c("NG" = "#40d1d6", "HG" = "#ff6e38")
# 生成每个 group 的椭圆坐标
ellipses <- pc_df %>%
  group_by(group) %>%
  do({
    df <- .
    cov_mat <- cov(cbind(df$PC1, df$PC2))
    center <- colMeans(cbind(df$PC1, df$PC2))
    ell <- as.data.frame(ellipse(cov_mat, centre = center, level = 0.95))
    ell$group <- unique(df$group)
    ell
  })

p_main <- ggplot(pc_df, aes(PC1, PC2, color = group)) +
  geom_polygon(
    data = ellipses,
    aes(x = x, y = y, fill = group),
    alpha = 0.15, color = NA
  ) +
  geom_point(size = 4, alpha = 0.95) +
  scale_color_manual(values = color_map) +
  scale_fill_manual(values = color_map) +
  labs(
    x = paste0("PC1 (", expl[1], "%)"),
    y = paste0("PC2 (", expl[2], "%)"),
    color = NULL, fill = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right"
  )


x_lim <- range(ellipses$x)
y_lim <- range(ellipses$y)

expand_ratio <- 0.15
p_main <- p_main +
  scale_x_continuous(limits = c(x_lim[1] - diff(x_lim)*expand_ratio,
                                x_lim[2] + diff(x_lim)*expand_ratio)) +
  scale_y_continuous(limits = c(y_lim[1] - diff(y_lim)*expand_ratio,
                                y_lim[2] + diff(y_lim)*expand_ratio)) +
  coord_cartesian(clip = "off")

p_main


ggsave("ZQY PCA.pdf",p_main, width =8, height = 6)
