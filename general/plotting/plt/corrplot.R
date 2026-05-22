rm(list = ls())
# hongfan
library(corrplot)

directory = "pathtoyourdata"
setwd(directory)
list.files()

data <- read.csv("pathtoyourdata", header = TRUE)
data1 <- -log10(data[, 2:11])
rownames(data1) <- data$GO.Term

data2 <- as.matrix(data1)
color <- colorRampPalette(c("blue", "#faf16e", "#fc4b4b"))(100)


pdf("corrplot_GO.pdf", width = 7, height = 4)
corrplot(data2,
    is.corr = FALSE, col = color, bg = "#f8fcdb", outline = T,
    tl.col = "#4e4d4d", tl.srt = 45
)

dev.off()
