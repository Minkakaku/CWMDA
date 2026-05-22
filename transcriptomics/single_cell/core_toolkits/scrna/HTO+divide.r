# hongfan
# singlecellexpriement
rm(list = ls())
# droplet filtered
suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(scuttle)
    library(Biobase)
    library(DropletUtils)
    library(dplyr)
})

d <- "D:\\xinqinlaozuo\\220820_hep_scRNAseq(cell hashing)"
setwd(d)
fs <- list.files(path = "h5", pattern = ".h5")
lapply(fs, function(x) {
    df <- Read10X_h5("x")
    p <- str_split(x, "_", simplify = TRUE)[, 1]
    # bcrank <- barcodeRanks(df$"Gene Expression")
    # # Only showing unique points for plotting speed.
    # uniq <- !duplicated(bcrank$rank)
    # plot(bcrank$rank[uniq], bcrank$total[uniq],
    #     log = "xy",
    #     xlab = "Rank", ylab = "Total UMI count", cex.lab = 1.2
    # )
    # abline(h = metadata(bcrank)$inflection, col = "darkgreen", lty = 2)
    # abline(h = metadata(bcrank)$knee, col = "dodgerblue", lty = 2)
    # legend("bottomleft",
    #     legend = c("Inflection", "Knee"),
    #     col = c("darkgreen", "dodgerblue"), lty = 2, cex = 1.2
    # )
    set.seed(100)
    e.out <- emptyDrops(df$"Gene Expression", lower = 100)
    # summary(e.out$FDR <= 0.1)

    retained <- df$"Gene Expression"[, which(e.out$FDR <= 0.1)]
    # bcrank <- barcodeRanks(retained)
    # uniq <- !duplicated(bcrank$rank)
    # plot(bcrank$rank[uniq], bcrank$total[uniq],
    #     log = "xy",
    #     xlab = "Rank", ylab = "Total UMI count", cex.lab = 1.2
    # )

    # saveRDS(retained, "filtered_empty_droplets.rds")
    #######################################################################
    dd <- retained
    umi <- CreateSeuratObject(dd)
    # Add number of genes per UMI for each cell to metadata
    umi$log10GenesPerUMI <- log10(umi$nFeature_RNA) / log10(umi$nCount_RNA)
    # Compute percent mito ratio
    umi$mitoRatio <- PercentageFeatureSet(object = umi, pattern = "^mt-")
    umi$mitoRatio <- umi$mitoRatio / 100
    ff <- GetAssayData(umi, slot = "data")
    # write.csv(ff, "ff.csv")
    VlnPlot(umi,
        features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"),
        ncol = 3
    )
    # Create metadata dataframe
    metadata <- umi@meta.data
    # Add cell IDs to metadata
    metadata$cells <- rownames(metadata)
    # Rename columns
    metadata <- metadata %>%
        dplyr::rename(
            seq_folder = orig.ident,
            nUMI = nCount_RNA,
            nGene = nFeature_RNA
        )
    # Create sample column
    metadata$sample <- NA
    metadata$sample <- "one"
    # Add metadata back to Seurat object
    umi@meta.data <- metadata

    # Filter out low quality reads using selected thresholds
    filtered_seurat <- subset(
        x = umi,
        subset = (nUMI >= 3000) &
            (log10GenesPerUMI > 0.70) &
            (mitoRatio < 0.3)
    )

    # Output a logical vector for every gene on whether the more than zero counts per cell
    # Extract counts
    counts <- GetAssayData(object = filtered_seurat, slot = "counts")

    # Output a logical vector for every gene on whether the more than zero counts per cell
    nonzero <- counts > 0

    # Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
    keep_genes <- Matrix::rowSums(nonzero) >= 5

    # Only keeping those genes expressed in more than 10 cells
    filtered_counts <- counts[keep_genes, ]

    # Reassign to filtered Seurat object
    filtered_seurat <- CreateSeuratObject(filtered_counts,
        meta.data = filtered_seurat@meta.data, project = p
    )
    # saveRDS(filtered_seurat, "p48_filtered_seurat.RDS")

    htos <- df$"Antibody Capture"
    umis <- filtered_seurat
    # htos <- htos[, colSums(htos) > 0]
    # write.csv(rownames(htos), file = "htos.csv")
    # stop()
    bcs <- intersect(colnames(umis), colnames(htos))
    umis <- umis[, bcs]
    htos <- as.matrix(htos[, bcs])
    # p <- str_split(x, "_", simplify = TRUE)[, 1]
    # sce <- CreateSeuratObject(umis, project = p)
    sce <- umis
    rm(umis)
    sce <- NormalizeData(sce)
    sce <- FindVariableFeatures(sce, selection.method = "mean.var.plot")
    sce <- ScaleData(sce, features = VariableFeatures(sce))
    sce[["HTO"]] <- CreateAssayObject(counts = htos)
    sce <- NormalizeData(sce,
        assay = "HTO",
        normalization.method = "CLR"
    )
    sce <- HTODemux(sce,
        assay = "HTO",
        positive.quantile = 0.99
    )
    n2 <- names(table(Idents(sce)))
    for (i in seq_len(length(n2))) {
        c <- assign(n2[i], subset(sce, idents = n2[i]))
        saveRDS(c, file = paste(n2[i], ".rds", collapse = ""))
    }
})
