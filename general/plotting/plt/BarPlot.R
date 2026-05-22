rm(list = ls())
# hongfan
library(ggplot2)
library(grid)
library(gridExtra)
library(ggrepel)

# data prepare
directory <- "pathtoyourdata"
setwd(directory)
list.files()
DATA <- read.csv("pathtoyourdata", header = T)
data <- DATA[1:15, ] # 1:10

g2 <- ggplot(data = data) +
      geom_bar(aes(x = reorder(Pathway.Name, -P.value), y = -log2(P.value)),
            stat = "identity", fill = "#ACD5DF", width = 0.9,
            position = position_dodge(width = 0.8)
      ) +
      theme(
            axis.title.x = element_text(face = "bold", size = 10),
            axis.title.x.top = element_text(face = "bold"),
            axis.text.x = element_text(size = 9),
            axis.title.y = element_text(face = "bold", size = 11, vjust = 3),
            axis.text = element_text(face = "bold"),
            axis.text.y = element_text(size = 11),
            axis.ticks.y = element_blank(),
            axis.line = element_line(lineend = "round"),
            panel.background = element_blank(),
            plot.margin = unit(c(1, 1, 1, 2), "mm")
      ) +
      labs(x = "Number of genes", y = "-log2(P)") +
      geom_text_repel(aes(
            x = reorder(data$Pathway.Name, -P.value), y = 0,
            label = reorder(data$Pathway.Name, -P.value), fontface = "bold"
      ),
      size = 3.5, inherit.aes = T, hjust = 0.7,
      position = position_dodge(width = 0.8)
      ) +
      scale_x_discrete(
            breaks = reorder(data$Pathway.Name, -data$P.value),
            labels = as.character(data$Term.Candidate.Gene.Num)
      ) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_flip()
g2
pdf("p10 cm total-Pathway(ǰ15).pdf", width = 4.42, height = 5) ## 4.42  3.4
print(g2)
dev.off()
