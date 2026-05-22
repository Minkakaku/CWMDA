rm(list = ls())
library(vegan)
library(dplyr)
library(ggcor)
library(ggplot2)
# OTU表格
d = "/home/hongfan/Punchure/projects/ZJK/"
setwd(d)
df <- read.csv("/home/hongfan/Punchure/projects/ZJK/相关性分析数据zjk.csv",
    row.names = 1, check.names = FALSE
)

spec <- df[, c(1, 2, 3, 4)]
env <- df[, c(-1, -2, -3, -4)]

library(dplyr)
mantel <- mantel_test(spec, env,
    spec.dist.method = "euclidean",
    spec.select = list(
        Bacillariophyta = 1,
        Dinophyceae = 2,
        Chlorophyta = 3,
        Haptophyceae = 4
    )
) %>%
    mutate(
        rd = cut(r,
            breaks = c(-Inf, 0.2, 0.4, Inf),
            labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")
        ),
        pd = cut(p.value,
            breaks = c(-Inf, 0.01, 0.05, Inf),
            labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")
        )
    )
write.csv(mantel, "mantel,p,r.csv")
pdf("mantel_test.pdf")
quickcor(env, type = "upper") +
    geom_square() +
    anno_link(aes(colour = pd, size = rd), data = mantel) +
    scale_size_manual(values = c(0.5, 1, 2)) +
    scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
    guides(
        size = guide_legend(
            title = "Mantel's r",
            override.aes = list(colour = "grey35"),
            order = 2
        ),
        colour = guide_legend(
            title = "Mantel's p",
            override.aes = list(size = 3),
            order = 1
        ),
        fill = guide_colorbar(title = "Pearson's r", order = 3)
    )
dev.off()
