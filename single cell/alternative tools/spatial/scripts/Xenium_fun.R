install.packages("Seurat")
install.packages("harmony")
remotes::install_github('prabhakarlab/Banksy')
remotes::install_github('satijalab/seurat-wrappers')
devtools::install_github('immunogenomics/presto')
remotes::install_github("bnprks/BPCells/r")
install.packages("hdf5r")
install.packages("sf")
install.packages("arrow")

library(Seurat)
library(sf)
library(ggplot2)


# directory "Xenium_Prime_Breast_Cancer_FFPE_outs" contains those 3 required files and 1 folder

# for this big dataset, loading time is ~16mins
xenium.obj <- LoadXenium("./one-sample-analysis/Xenium_Prime_Breast_Cancer_FFPE_outs",
                         fov = "fov", segmentations="cell", flip.xy=TRUE)

metadata <- xenium.obj@meta.data
head(metadata, n=3)

options(repr.plot.width=10, repr.plot.height=16)

# max.cutoff="q95" - Color is capped at 95th percentile
ImageFeaturePlot(xenium.obj, fov = "fov",
                 features = c("nCount_Xenium"),
                 max.cutoff="q95")+
  scale_y_reverse()