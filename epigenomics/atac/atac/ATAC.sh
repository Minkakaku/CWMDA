#!/bin/bash
#PBS -l nodes=1:ppn=20,walltime=100:00:00
#PBS -j oe
#PBS -N zl

#atac peak calling progress

cd /data/gpfs02/wmo/work/DNA-seq/ATAC/JMY

#skewer -f sanger -t 20 -z -m pe -x /share/apps/pepatac/tools/NexteraPE-PE.fa \
#--quiet -o $id/$id $id/$id\_1.clean.fq.gz $id/$id\_2.clean.fq.gz

trim_galore -j 4 -q 10 --phred33 --length 35 -e 0.1 --stringency 3 --paired --fastqc --gzip \
-o 02.trim/ 01.cleandata/$id/$id\_1.clean.fq.gz 01.cleandata/$id/$id\_2.clean.fq.gz


## QC
cat name |while read id;do
bsub -q normal_1week -n 10 -e $id.err -o $id.log -J $id \
#fastqc -t 10 01.cleandata/$id/$id\_1.clean.fq.gz -o 03.QC/
#fastqc -t 10 01.cleandata/$id/$id\_2.clean.fq.gz -o 03.QC/
## align
bowtie2 -p 10 --very-sensitive -X 500 --rg-id $id --rg SM:$id \
-x /data/gpfs01/wmo/genome/reference/bowtie2_mm10/mm10 \
-1 02.trim/$id\_1.clean_val_1.fq.gz -2 02.trim/$id\_2.clean_val_2.fq.gz \
-S 03.align/$id.sam

samtools view -@ 10 -bS 03.align/$id.sam -o 03.align/$id.bam 
samtools sort -@ 10 -o 03.align/$id.sorted.bam 03.align/$id.bam
samtools index -@ 10 03.align/$id.sorted.bam

## picard duplicates
# 这里选用sambamba来去重复
cat name |while read id;do
bsub -q normal -n 10 -e $id.err -o $id.log -J $id \

bedtools bamtobed -i 03.align/$id.sorted.bam > 03.align/$id.sorted.bed
samtools flagstat 03.align/$id.sorted.bam > 03.align/$id.sorted.stat

sambamba markdup -t 10 --overflow-list-size 200000 --tmpdir='04.sambamba/' \
-r 03.align/$id.sorted.bam 04.sambamba/$id.rmdup.bam

samtools flagstat 04.sambamba/$id.rmdup.bam >04.sambamba/$id.rmdup.stat
## 接下来只保留两条reads要比对到同一条染色体(Proper paired) ，还有高质量的比对结果(Mapping quality>=30)
## 顺便过滤 线粒体reads
cat name |while read id;do
bsub -q normal_1week -n 10 -e $id.err -o $id.log -J $id \

samtools view -@ 10 -h 04.sambamba/$id.rmdup.bam|grep -v chrM >04.sambamba/$id.rmMT.bam
samtools sort -O bam -@ 10 -o 04.sambamba/$id.last.bam 04.sambamba/$id.rmMT.bam
samtools index -@ 10 04.sambamba/$id.last.bam
bedtools bamtobed -i 04.sambamba/$id.last.bam > 04.sambamba/$id.last.bed

## 原来的去除线粒体染色质和转换成bed
#samtools view -h 04.sambamba/$id.rmdup.bam|removeChrom.py - - chrM|samtools view -b - > 04.sambamba/$id.rmdup.rmMT.bam
#samtools view -h /home/casual/ATAC_seq_out/zhangliang/$id/$id.sorted.dedup.rmMT.bam | \
#SAMtoBED.py -i - -o /home/casual/ATAC_seq_out/zhangliang/$id/$id.bed -x -v
source activate python3
 
macs2 callpeak -t 04.sambamba/$id.last.bed -f BED -g mm --outdir 05.callpeak/ -n $id --shift -75 --extsize 150 \
--nomodel --call-summits --nolambda --keep-dup all -p 0.01
done

bsub -q normal_1week  -n 1 -e K9.err -o K9.log -J K9 \
macs2 callpeak -t 04.sambamba/ATACP10BGMT_TKD181000262.last.bam -c 04.sambamba/ATACP10BGWT_TKD181000260.last.bam \
-f BAM -g mm --keep-dup all -p 0.01 -n del_P10_MTWT -B --outdir 05.callpeak/

bsub -q normal_1week  -n 1 -e K9.err -o K9.log -J K9 \
macs2 callpeak -t 04.sambamba/ATACP10BGWT_TKD181000260.rmdup.bam -c 04.sambamba/ATACP7BGWT_TKD181000259.last.bam \
-f BAM -g mm --keep-dup all -p 0.01 -n del_P7_P10_WT -B --outdir 05.callpeak/

# bw for ucsc
cat name |while read id;do
bsub -q normal_1week -n 20 -e $id.err -o $id.log -J $id \
bamCoverage -p 20 -b 04.sambamba/$id.last.bam -o 06.heatmap/$id.bw 

## bw file from compare
bsub -q normal_1week -n 1 -e bamc.err -o bamc.log -J bamc \
bamCompare -b1 04.Align/D3_IPSeh1_FKDL202611416-1a.sorted.bam -b2 04.Align/D0_IPIgG_FKDL202611414-1a.sorted.bam \
-o 04.Align/log2_D3_IgG.bw

# deeptools for olig2 heatmap
bsub -q normal_1week -n 20 -e iol.err -o iol.log -J iol \
computeMatrix reference-point -p 20 --referencePoint TSS -b 5000 -a 5000 \
-R 05.callpeak/ATACP10BGCT_TKD181000261_narrowPeak.bed \
-S 06.heatmap/ATACP10BGMT_TKD181000262.bw 06.heatmap/ATACP10BGWT_TKD181000260.bw \
--skipZeros -o 06.heatmap/P10MT_WT_point_TSS.mat.gz

bsub -q normal_1week -n 1 -e iol.err -o iol.log -J iol \
plotHeatmap -m 06.heatmap/P7P10WT_point_TSS.mat.gz --heatmapHeight 10 --colorMap vlag \
--legendLocation upper-right -out 06.heatmap/P7P10WT_point_TSS.png

#cat name |while read id;do
# motif
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}'  05.callpeak/del_P10_MTWT_peaks.narrowPeak>05.callpeak/del_P10_MTWT_peaks.narrowPeak.bed

#cat name |while read id;do
bsub -q normal_1week -n 20 -e motif.err -o motif.log -J motif \
findMotifsGenome.pl 05.callpeak/del_P10_MTWT_peaks.narrowPeak.bed mm10 07.motif/P10_MTWT/ -size 200 -mask -p 20

cd 05.callpeak/
#cat ../name |while read id;do
annotatePeaks.pl del_P10_MTWT_peaks.narrowPeak.bed mm10 1>del_P10_MTWT_peaks.peakAnn.xls 2>del_P10_MTWT_peaks.annLog.txt 
