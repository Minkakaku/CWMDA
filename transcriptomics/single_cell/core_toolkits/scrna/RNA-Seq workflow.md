# RNA-Seq workflow: gene-level exploratory analysis and differential expression

## Align
**STAR**
If you have performed a single-end experiment, you would only have one file per ID
```bash
for f in 'cat files'; do STAR --genomeDir ../STAR/ENSEMBL.homo_sapiens.release-75 \
--readFilesIn fastq/$f\_1.fastq fastq/$f\_2.fastq \
--runThreadN 12 --outFileNamePrefix aligned/$f.; done
```
SAMtools10 was used to generate BAM files. The –@ flag can be used to allocate additional threads

```bash
for f in 'cat files'; do samtools view -bS aligned/$f.Aligned.out.sam \
-o aligned/$f.bam; done
```

## featurecounts
```r
sampleTable <- read.csv(csvfile,row.names=1)

## SampleName cell dex albut Run avgLength Experiment Sample BioSample
## SRR1039508 GSM1275862 N61311 untrt untrt SRR1039508 126 SRX384345 SRS508568 SAMN02422669
## SRR1039509 GSM1275863 N61311 trt untrt SRR1039509 126 SRX384346 SRS508567 SAMN02422675
## SRR1039512 GSM1275866 N052611 untrt untrt SRR1039512 126 SRX384349 SRS508571 SAMN02422678
## SRR1039513 GSM1275867 N052611 trt untrt SRR1039513 87 SRX384350 SRS508572 SAMN02422670
## SRR1039516 GSM1275870 N080611 untrt untrt SRR1039516 120 SRX384353 SRS508575 SAMN02422682
## SRR1039517 GSM1275871 N080611 trt untrt SRR1039517 126 SRX384354 SRS508576 SAMN02422673
## SRR1039520 GSM1275874 N061011 untrt untrt SRR1039520 101 SRX384357 SRS508579 SAMN02422683
## SRR1039521 GSM1275875 N061011 trt untrt SRR1039521 98 SRX384358 SRS508580 SAMN02422677
```

We indicate in Bioconductor that these files are BAM files using the BamFileList function from the Rsamtools package that provides an R interface to BAM files
```r
library("Rsamtools")
bamfiles <- BamFileList(filenames, yieldSize=2000000)
```
See the `seqlevelsStyle` function in the `GenomeInfoDb` package for solutions. 
```r
seqinfo(bamfiles[1])
```
```r
gtffile <- file.path(dir,"Homo_sapiens.GRCh37.75_subset.gtf")
(txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character()))
```
```r
(ebg <- exonsBy(txdb, by="gene"))
```

## Read counting step

```r
library("GenomicAlignments")
library("BiocParallel")

register(SerialParam())
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
 mode="Union",
 singleEnd=FALSE,
 ignore.strand=TRUE,
 fragments=TRUE )
```
Setting singleEnd to FALSE indicates that 
the experiment produced paired-end reads, and we want to count a pair of reads (a fragment) only once toward the count for a gene.

condition is a column in colData(dds) that specifies which of two (or more groups) the samples belong to. For the airway experiment, we will specify ~ cell + dex meaning that we want to test for the effect of dexamethasone (dex) controlling for the effect of different cell line (cell).

```r
Note: it is prefered in R that the first level of a factor be the reference level (e.g. control, or untreated samples), so we 
can relevel the dex factor like so:
se$dex <- relevel(se$dex, "untrt")
se$dex
```

## DESeq2
```r
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ cell + dex)
countdata <- assay(se)
coldata <- colData(se)
(ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
 colData = coldata,
 design = ~ cell + dex))
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)
```
## PCA
```r
sampleDists <- dist( t( assay(rld) ) )

library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
 clustering_distance_rows=sampleDists,
 clustering_distance_cols=sampleDists,
 col=colors)
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
 clustering_distance_rows=poisd$dd,
 clustering_distance_cols=poisd$dd,
 col=colors)

# PCA
plotPCA(rld, intgroup = c("dex", "cell"))
(data <- plotPCA(rld, intgroup = c( "dex", "cell"), returnData=TRUE))
percentVar <- round(100 * attr(data, "percentVar"))
library("ggplot2")
ggplot(data, aes(PC1, PC2, color=dex, shape=cell)) + geom_point(size=3) +
 xlab(paste0("PC1: ",percentVar[1],"% variance")) +
 ylab(paste0("PC2: ",percentVar[2],"% variance"))

# MDS plot
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld)))
ggplot(mds, aes(X1,X2,color=dex,shape=cell)) + geom_point(size=3)
mdsPoisData <- data.frame(cmdscale(samplePoisDistMatrix))
mdsPois <- cbind(mdsPoisData, as.data.frame(colData(dds)))
ggplot(mdsPois, aes(X1,X2,color=dex,shape=cell)) + geom_point(size=3)
```

## Differential expression analysis
```r
dds <- DESeq(dds)
(res <- results(dds))
mcols(res, use.names=TRUE)
res.05 <- results(dds, alpha=.05)
table(res.05$padj < .05)
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)
results(dds, contrast=c("cell", "N061011", "N61311"))
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])
# Plotting results
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene=topGene, intgroup=c("dex"))
data <- plotCounts(dds, gene=topGene, intgroup=c("dex","cell"), returnData=TRUE)
ggplot(data, aes(x=dex, y=count, color=cell)) +
 scale_y_log10() +
 geom_point(position=position_jitter(width=.1,height=0), size=3)
ggplot(data, aes(x=dex, y=count, fill=dex)) +
 scale_y_log10() +
 geom_dotplot(binaxis="y", stackdir="center")
ggplot(data, aes(x=dex, y=count, color=cell, group=cell)) +
 scale_y_log10() + geom_point(size=3) + geom_line()

plotMA(res, ylim=c(-5,5))
plotMA(resLFC1, ylim=c(-5,5))
topGene <- rownames(resLFC1)[which.min(resLFC1$padj)]
with(resLFC1[topGene, ], {
 points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
 text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
hist(res$pvalue[res$baseMean > 1], breaks=0:20/20, col="grey50", border="white")
# Gene clustering
library("genefilter")
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),20)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("cell","dex")])
pheatmap(mat, annotation_col=df)
# Independent filtering
qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
bins <- cut(resLFC1$baseMean, qs) 
levels(bins) <- paste0("~",round(signif(.5*qs[-1] + .5*qs[-length(qs)],2)))
ratios <- tapply(resLFC1$pvalue, bins, function(p) mean(p < .05, na.rm=TRUE))
barplot(ratios, xlab="mean normalized count", ylab="ratio of small p values")

# Annotating and exporting results
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
res$symbol <- mapIds(org.Hs.eg.db,
 keys=row.names(res),
column="SYMBOL",
 keytype="ENSEMBL",
 multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
 keys=row.names(res),
column="ENTREZID",
 keytype="ENSEMBL",
 multiVals="first")
resOrdered <- res[order(res$padj),]
head(resOrdered)
# Exporting results
resOrderedDF <- as.data.frame(resOrdered)[1:100,]
write.csv(resOrderedDF, file="results.csv")
library("ReportingTools")
htmlRep <- HTMLReport(shortName="report", title="My report",
 reportDirectory="./report")
publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
browseURL(url)

# Plotting fold changes in genomic space
(resGR <- results(dds, lfcThreshold=1, format="GRanges"))
resGR$symbol <- mapIds(org.Hs.eg.db, names(resGR), "SYMBOL", "ENSEMBL")
library("Gviz")
window <- resGR[topGene] + 1e6
strand(window) <- "*"
resGRsub <- resGR[resGR %over% window]
naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)
sig <- factor(ifelse(resGRsub$padj < .1 & !is.na(resGRsub$padj),"sig","notsig"))

options(ucscChromosomeNames=FALSE)
g <- GenomeAxisTrack()
a <- AnnotationTrack(resGRsub, name="gene ranges", feature=sig)
d <- DataTrack(resGRsub, data="log2FoldChange", baseline=0,
 type="h", name="log2 fold change", strand="+")
plotTracks(list(g,d,a), groupAnnotation="group", notsig="grey", sig="hotpink")
# Removing hidden batch effects
library("sva")
dat <- counts(dds, normalized=TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~ dex, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv=2)
par(mfrow=c(2,1),mar=c(3,5,3,1))
stripchart(svseq$sv[,1] ~ dds$cell,vertical=TRUE,main="SV1")
abline(h=0)
stripchart(svseq$sv[,2] ~ dds$cell,vertical=TRUE,main="SV2")
abline(h=0)
ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + dex
ddssva <- DESeq(ddssva)

# Time course experiments
library("fission")
data("fission")
ddsTC <- DESeqDataSet(fission, ~ strain + minute + strain:minute)
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ strain + minute)
resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$symbol
head(resTC[order(resTC$padj),],4)
data <- plotCounts(ddsTC, which.min(resTC$padj),
 intgroup=c("minute","strain"), returnData=TRUE)
ggplot(data, aes(x=minute, y=count, color=strain, group=strain)) +
 geom_point() + stat_smooth(se=FALSE,method="loess") + scale_y_log10()
resultsNames(ddsTC)
res30 <- results(ddsTC, name="strainmut.minute30", test="Wald")
res30[which.min(resTC$padj),]
betas <- coef(ddsTC)
colnames(betas)
library("pheatmap")
topGenes <- head(order(resTC$padj),20)
mat <- betas[topGenes, -c(1,2)]
thr <- 3
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
 cluster_col=FALSE)
```


