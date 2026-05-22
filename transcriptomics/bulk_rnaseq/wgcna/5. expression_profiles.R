# 清空工作区
rm(list = ls())

# 加载必要的包
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(purrr)

# 加载数据
paths <- c(
  "WGCNA\\opt_inputdata\\opt_inputdata.RData",
  "WGCNA\\results\\networkConstruction_auto.RData"
)
for (path in paths) {
  load(path)
}
# 创建模块数据框
module_df <- tibble(
  gene_id = names(moduleLabels),
  colors = moduleColors
)

# 定义感兴趣的颜色
colors_of_interest <- c("green", "turquoise", "tan")

# 筛选感兴趣模块并提取表达数据
subexpr <- datExpr[
  ,
  module_df %>% 
    filter(colors %in% colors_of_interest) %>% 
    pull(gene_id)
]
head(subexpr)

# 转换数据格式并添加模块信息
submod_df <- subexpr %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene_id") %>% 
  pivot_longer(-gene_id) %>% 
  left_join(module_df %>% rename(module = colors), by = "gene_id")

# 绘制表达谱图
submod_df %>% 
  ggplot(aes(x = name, y = value, group = gene_id, color = module)) +
  geom_line(alpha = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment", y = "normalized expression")
