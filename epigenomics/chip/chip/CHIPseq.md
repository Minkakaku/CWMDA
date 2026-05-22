# Methods for ChIP-seq analysis: A practical workflow and advanced applications 

**H3K4me1 H3K27ac** associated with enhancer regions  
**H3K4me3** associated with promoter regions  
**H3K36me3** associated with transcribed regions in gene bodies  
**H3K27me3** associated with Polycomb repression  
**H3K9me3** associated with heterochromatin  

## Technical considerations of ChIP-seq analysis for histone modifications
The reliability of a ChIP analysis is governed by antibody quality, including specificity and signal-to-noise ratio (S/N)

**H3K27ac** and **H3K4me1**, produce sharp peaks, but sometimes construct broadly enriched regions called “super enhancers”  

## Read mapping

CHIP-seq: Bowtie2、BWA用的比较多  
RNA-seq: Tophat、Bsmap  
甲基化：BS-seeker  

Bowtie2
如果目的是对齐两个非常大的序列（例如两个基因组），请考虑使用MUMmer

自建索引

```bash
{
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz 
tar -zxvf chromFa.tar.gz 
cat *.fa > mm10.fa
bowtie2-build mm10.fa mm10
}
```
单末端
```bash
{
bowtie2 -p 10 -x genome_index -U input.fq | samtools sort -O bam -@ 10 -o
}
```
双末端
```bash
{
bowtie2 -p 10 -x genome_index -1 input_1.fq -2 input_2.fq | samtools sort -O bam -@ 10 -o
}
```

## Peak calling
MASC2
```bash
{
# regular peak calling：
macs2 callpeak -t ChIP.bam -c Control.bam -f BAM -g hs  --shift -100 --extsize 200 -n test -B -q 0.01
# broad peak calling:
macs2 callpeak -t ChIP.bam -c Control.bam -g hs --shift -100 --extsize 200 --broad --broad-cutoff 0.1
}
```

## 使用 trim_galore 进行数据清洗
```bash
{
ls ~/chipseq/raw/*.gz | while read id;
do
nohup trim_galore -j 10 -q 20 --phred33 --length 25 -e 0.1 --stringency 4 -o ~/chipseq/clean $id 1>trim_galore.log & 
done
# 出错时查看 trim_galore.log 可能是参数等输错了。
# -j 线程
# -q -quality<int> 设定phred quality阈值。默认20（99%的read质量），如果测序深度较深，可以设定25
# --phred33 可以选择-phred33或者-phred64，表示测序平台使用的Phred quality score
# --length 设定输出reads长度阈值，小于设定值会被抛弃
# --fastqc 同时做质控
# --fastqc_args 可以设置 fastqc 的参数
}
```
#### 输入文件参数：
`-t`:实验组，IP的数据文件
`-c`: 对照组, control 或 mock(非特异性抗体，如IgG)组
`-f`：指定输入文件的格式，默认是自动检测输入数据是什么格式，支持bam,sam,bed等
`-g`:有效基因组大小，由于基因组序列的重复性，基因组实际可以mapping的大小小于原始的基因组。这个参数要根据实际物种计算基因组的有效大小。软件里也给出了几个默认的-g 值：hs -- 2.7e9表示人类的基因组有效大小(UCSC human hg18 assembly).
hs:	2.7e9
mm:	1.87e9
ce:	9e7
dm:	1.2e8
`-S/–TSIZE`:测序读长；如果不设定，MACS 利用输入的前10个序列自动检测；
`–BW`：湿实验中，声波打断基因组的片段长度，用来建立模型；

#### 输出文件参数：
`--outdir`:输出结果的存储路径  
`-n`:输出文件名的前缀
`-B/--bdg`:输出bedgraph格式的文件，输出文件以NAME+'_treat_pileup.bdg' for treatment data, NAME+'_control_lambda.bdg' for local lambda values from control显示。  

#### peak calling 参数:
`--Q/–QVALUE`：qvalue (minimum FDR)设定call significant regions的阈值；默认，0.01，对于 broad marks(组蛋白修饰的chipseq)，可以使用0.05；Q-values are calculated from p-values using Benjamini-Hochberg procedure.
`-M/–MFOLD`：构建模型时，enrichment regions 选用标准（MFOLD range of high-confidence enrichment ratio against background to build model）;DEFAULT:5,50 means using all regions not too low (>5) and not too high (<50) to build paired-peaks model. MACS 无法找到超过100 regions 用来构建模型时，只有设定–fix-bimodal情况下，MACS 会调用参数–extsize。
`--broad`:peak有narrow peak和broad peak, 设置时可以call broad peak 的结果文件。  
`--broad-cutoff`:和pvalue、以及qvalue相似  
`--nolambda`: 不要考虑在峰值候选区域的局部偏差/λ  
`–SLOCAL, –LLOCAL`：设定两个水平检测peak 区域，从而计算最大λ作为local λ。默认，MACS 采用1000bp为small local region(–slocal)，10000bps为large local region(–llocal)计算开放染色体区域的偏差。区域设置的太小，尖峰会掩盖掉旁边显著性的峰。

#### Shift 模型参数：
`--nomodel`:这个参数和extsize、shift是配套使用的，有这个参数才可以设置extsize和shift。对`双端测序`而言，它本身测的就是文库的两端，因此不用建立模型和偏倚。我们只需要对 MACS 设置参数 --nomodel 就能略过双峰模型建立的过程。
`--extsize`:当设置了nomodel时，MACS会用--extsize这个参数从5'->3'方向扩展reads修复fragments。比如说你的转录因子结合范围200bp，就设置这个参数是200。
`--shift`:当设置了--nomodel，MACS用这个参数从5' 端移动剪切，然后用--extsize延伸，如果--shift是负值表示从3'端方向移动。建议ChIP-seq数据集这个值保持默认值为0，对于检测富集剪切位点如DNAsel数据集设置为EXTSIZE的一半。
示例：
想找富集剪切位点，如DNAse-seq，所有5'端的序列reads应该从两个方向延伸，如果想设置移动的窗口是200bp，参数设置如下：
--nomodel --shift -100 --extsize 200
对nucleosome-seq数据，用核小体大小的一半进行extsize,所以参数设置如下：
--nomodel --shift 37 --extsize 73
`--call-summits`:MACS利用此参数重新分析信号谱，解析每个peak中包含的subpeak。对相似的结合图谱，推荐使用此参数，当使用此参数时，输出的subpeak会有相同的peak边界，不同的绩点和peak summit poisitions.

#### ATAC-Seq call peaks示例
ATAC-seq关心的是在哪切断，断点才是peak的中心，所以使用shift模型，--shift -75或-100
```bash
{
macs2 callpeak -t H1hesc.final.bam -n sample --shift -100 --extsize 200 --nomodel -B --SPMR -g hs --outdir Macs2_out 2> sample.macs2.log
}
```

### bdgdiff
通过bdgdiff子命令来进行差异peak分析， 该命令不需要基于已有的peak calling结果，只需要输入每个样本对应的bedGraph格式的文件。需要注意的是，该命令只针对两个样本间的差异peak进行设计，适用于没有生物学重复的情况。对于使用macs2来进行差异peak的完整流程，官方给出了详细的说明文档，链接如下
`https://github.com/taoliu/MACS/wiki/Call-differential-binding-events`

可以分为以下3步
1. 预测插入片段长度
通过predictd子命令可以预测样本的fragment size，命令如下
```bash
{
macs2 predictd -i input.bam
}
```
2. peak  calling
在peak calling时，需要添加-B参数，这样才可以输出样本对应的bedgraph文件，同时需要保证peak  calling时采用一致的--extsize的值，就是第一步预测出来的数值，取多个样本的均值即可。官方也给出了推荐值，对于大多数的转录因子chip_seq数据，推荐值为200， 对于大部分组蛋白修饰的chip_seq数据，推荐值为147，命令如下
```bash
{
# condition1
macs2 callpeak -B -t cond1_ChIP.bam -c cond1_Control.bam -n cond1 --nomodel --extsize 120
# condition2
macs2 callpeak -B -t cond1_ChIP.bam -c cond1_Control.bam -n cond1 --nomodel --extsize 120
}
```
在运行这一步的时候，会输出每个样本过滤之后的reads数目，示意如下
```bash
{
# tags after filtering in treatment: 19291269
# tags after filtering in control: 12914669
}
```
这个数值在差异分析中会用到，所以要记录下来, 选择最小的tags。
3. 差异peak分析
命令如下
```bash
{
macs2 bdgdiff --t1 cond1_treat_pileup.bdg --c1 cond1_control_lambda.bdg --t2 cond2_treat_pileup.bdg \
--c2 cond2_control_lambda.bdg --d1 12914669 --d2 14444786 -g 60 -l 120 --o-prefix diff_c1_vs_c2
}
```

## BAM 转为UCSC可视化需要的 BIGWIG 格式
deeptools bamCoverage
```bash
{
bamCoverage -e 170 -bs 10 -b test.bam -o test.bw 
}
```
这里的参数-e/--extendReads拓展了原来的read长度，-bs/--binSize 设置分箱的大小。

## MACS2 输出的 BedGraph 文件转为 Bigwig 文件，然后通过 UCSC 查看
首先，使用 macs2 bdgcmp 得到 FE 或者 logLR 转化后的文件 ：
```bash
{
# 使用 macs2 bdgcmp 得到 FE 转化后的文件
# -m FE 计算富集倍数降低数据噪音
macs2 bdgcmp -t NAME_treat_pileup.bdg -c NAME_control_lambda.bdg --outdir ../test -o NAME_FE.bdg -m FE

# 使用 macs2 bdgcmp 得到 logLR 转化后的文件
# -p 为避免log的时候input值为0时出现error，赋值为0.00001
macs2 bdgcmp -t NAME_treat_pileup.bdg -c NAME_control_lambda.bdg --outdir ../test -o NAME_logLR.bdg -m logLR -p 0.00001
}
```
**bdg2bw**
```{
# bedgraph to bigwig 
sh /scripts/bdg2bw.sh NAME_FE.bdg /index/hg38_UCSC.chromInfo.txt > NAME_FE.bw
sh /scripts/bdg2bw.sh NAME_logLR.bdg /index/hg38_UCSC.chromInfo.txt > NAME_logLR.bw
}
```
```bash
{
#!/bin/bash

# check commands: slopBed, bedGraphToBigWig and bedClip

which bedtools &>/dev/null || { echo "bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
which bedGraphToBigWig &>/dev/null || { echo "bedGraphToBigWig not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }
which bedClip &>/dev/null || { echo "bedClip not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }

# end of checking

if [ $# -lt 2 ];then
    echo "Need 2 parameters! <bedgraph> <chrom info>"
    exit
fi

F=$1
G=$2

bedtools slop -i ${F} -g ${G} -b 0 | bedClip stdin ${G} ${F}.clip

LC_COLLATE=C sort -k1,1 -k2,2n ${F}.clip > ${F}.sort.clip

bedGraphToBigWig ${F}.sort.clip ${G} ${F/bdg/bw}

rm -f ${F}.clip ${F}.sort.clip
}
```
染色体长度文件:http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/chromInfo.txt.gz


## 鉴定可重复的Peaks
1. 用Bedtools进行简单的overlap合并重复样本
2. IDR（Irreproducibility Discovery Rate）的方法获得高重复性的peaks

### 1. Overlapping peaks using bedtools
```bash
{
bedtools intersect [OPTIONS] -a <bed/gff/vcf/bam> -b <bed/gff/vcf/bam>
}
```
`-a`: 参数后加重复样本1（A）
`-b`：参数后加重复样本2（B），也可以加多个样本
其他常用参数解释和图解如下：
`-wo`：Write the original A and B entries plus the number of base pairs of overlap between the two features.
`-wa`：Write the original entry in A for each overlap.
`-v`：Only report those entries in A that have no overlaps with B
如果只有-a和-b参数，返回的是相对于A的overlaps。加上参数-wo返回A和B原始的记录加上overlap的记录。参数-wa返回每个overlap中A的原始记录。

### 2. Irreproducibility Discovery Rate (IDR)

> 建议使用IDR时，MACS2 call peaks的步骤参数设置不要过于严格，以便鉴定出更多的peaks。
> 使用IDR需要先对MACS2的结果文件narrowPeak根据-log10(p-value)进行排序。

```bash
{
# Call peaks
macs2 callpeak -t  sample.final.bam -n sample --shift -100 --extsize 200 --nomodel -B --SPMR -g hs --outdir Macs2_out 2> sample.macs2.log
#Sort peak by -log10(p-value)
sort -k8,8nr NAME_OF_INPUT_peaks.narrowPeak > macs/NAME_FOR_OUPUT_peaks.narrowPeak
}
```
**使用IDR示例**
```bash
{
idr --samples sample_Rep1_sorted_peaks.narrowPeak sample_Rep2_sorted_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file sample-idr \
--plot \
--log-output-file sample.idr.log
}
```
参数说明:
`--samples`:narrowPeak的输入文件（重复样本）
`--input-file-type`：输入文件格式包括narrowPeak,broadPeak,bed
`--rank p.value`：以p-value排序
`--output-file`: 输出文件路径
`--plot`：输出IDR度量值的结果

wc -l *-idr 计算下common peaks的个数，接着可再计算下与总peaks的比率。

## 用ChIPseeker对peaks进行注释和可视化

### 1. 制作TxDb方法示例：
```R
{
# makeTxDbFromUCSC
library(GenomicFeatures)
hg19.refseq.db <- makeTxDbFromUCSC(genome="hg19", table="refGene")
# makeTxDbFromGFF
download.file("ftp://ftp.ebi.ac.uk/pub/databases/pombase/pombe/Chromosome_Dumps/gff3/schizosaccharomyces_pombe.chr.gff3", "schizosaccharomyces_pombe.chr.gff3")require(GenomicFeatures)
spombe <- makeTxDbFromGFF("schizosaccharomyces_pombe.chr.gff3")
}
```
### 2. 加载包：
```R
{
# 下载人的基因和lincRNA的TxDb对象
biocLite("org.Hs.eg.db")
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
biocLite("TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts")
biocLite("clusterProfiler")
#载入各种包
library("ChIPseeker")
library(clusterProfiler)
library("org.Hs.eg.db")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library("TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts")
lincRNA_txdb=TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts
}
```

### 3. 读取文件
```R
{
nanog <- readPeakFile("./idr_out.bed/nanog_idr-bed")
pou5f1 <- readPeakFile("./idr_out.bed/pou5f1_idr-bed")
}
```

### 4. 注释peaks
```R
{
# 制作多个样本比较的list
peaks <- list(Nanog=nanog,Pou5f1=pou5f1)
# promotor区间范围可以自己设定
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peaks, getTagMatrix, windows=promoter)
#annotatePeak传入annoDb参数,可进行基因ID转换（Entrez，ENSEMBL，SYMBOL，GENENAME）
peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Hs.eg.db")
}
```

### 5. visualize
```R
{
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList,title="Distribution of transcription factor-binding loci n relative to TSS")
}
```
### 6. enrichment
```R
{
# Create a list with genes from each sample
gene = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
# Run GO enrichment analysis 
ego <- enrichGO(gene = entrez, 
                    keytype = "ENTREZID", 
                    OrgDb = org.Hs.eg.db, 
                    ont = "BP", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, 
                    readable = TRUE)
# Dotplot visualization
dotplot(ego, showCategory=50)
# Multiple samples KEGG analysis
compKEGG <- compareCluster(geneCluster = gene, 
                         fun = "enrichKEGG",
                         organism = "human",
                         pvalueCutoff  = 0.05, 
                         pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 20, title = "KEGG Pathway Enrichment Analysis")
}
```

```R
{
# Output peakAnnolist file
save(peakAnnoList,file="peakAnnolist.rda")
write.table(as.data.frame(peakAnnoList$Nanog),file="Nanog.PeakAnno",sep='t',quote = F)
# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
write.csv(cluster_summary, "results/clusterProfiler_Nanog.csv")
}
```

## other ways like DROMPAplus: a pipeline tool for ChIP-seq analysis


