#!/bin/bash
#PBS -l nodes=1:ppn=20,walltime=100:00:00
#PBS -j oe
#PBS -N name
# mouse

cd /data/gpfs02/wmo/work/DNA-seq/Chip-seq/csx_bg

#chip peak calling progress

mkdir 03.QC/ 04.Align/ 05.MACS/ 06.heatmap/ 07.motif/ 08.diff_heatmap

cat name |while read id;do
bsub -q normal_1day_new -n 20 -e $id.err -o $id.log -J $id \
## adapter trimmer 
#skewer -f sanger -t 20 -z -m pe -x /share/apps/pepatac/tools/NexteraPE-PE.fa --quiet -o $id /$id/$id_1.clean.fq.gz /$id/$id_2.clean.fq.gz
## quality control 
fastqc -t 20 02.cleandata/$id/$id\_1.clean.fq.gz -o 03.QC/
fastqc -t 20 02.cleandata/$id/$id\_2.clean.fq.gz -o 03.QC/
# build index 
#bowtie2-build  ~/genome/mm10/mm10.fa ~/genome/reference/bowtie2_mm10/mm10
## align the read
bowtie2 -p 20 --very-sensitive -X 500 -x ~/genome/reference/bowtie2_mm10/mm10 -1 02.cleandata/$id/$id\_1.clean.fq.gz -2 02.cleandata/$id/$id\_2.clean.fq.gz \
-S 04.Align/$id.sam

cat name |while read id;do
bsub -q normal_1day_new -n 20 -e $id.err -o $id.log -J $id \

#--rg-id MT-K9 --rg SM:MT-K9 
samtools view -@ 20 -b -S 04.Align/$id.sam -o 04.Align/$id.bam 
samtools sort -@ 20 -o 04.Align/$id.sorted.bam 04.Align/$id.bam  #^N
samtools index -@ 20 04.Align/$id.sorted.bam #^N

ls *bam | while read id; do samtools markdup -r $id $(basename $id ".bam").sm.bam & done

#peakcalling narrowpeak
conda activate python3
cat name |while read id;do
bsub -q normal_1day_new -n 1 -e $id.err -o $id.log -J $id \
macs2 callpeak -t 04.Align/$id.sorted.bam -f BAM -g 1.87e9 -n $id -B -p 0.01 --outdir 05.MACS/

## input
bsub -q normal_1week -n 1 -e macs.err -o macs.log -J macs \
macs2 callpeak -t 04.Align/BG-h3k9me3-ChIP-CT_TKD180602102.sorted.bam -c 04.Align/BG-input-MT_TKD180602101.sorted.bam \
-f BAM -n CT_input -B -p 0.01 --outdir 06.heatmap/

macs2 callpeak -t 04.Align/BG-h3k9me3-ChIP-MT_TKD180602103.sorted.bam -c 04.Align/BG-input-MT_TKD180602101.sorted.bam \
-f BAM -n MT_Input -B -p 0.01 --outdir 06.heatmap/

macs2 callpeak -t 04.Align/BG-h3k9me3-ChIP-MT_TKD180602103.sorted.bam -c 04.Align/BG-h3k9me3-ChIP-MT_TKD180602103.sorted.bam \
-f BAM -n diff_MT_CT -B -p 0.01 --outdir 06.heatmap/


#macs2 callpeak -t /home/casual/ATAC_seq_out/mowei/$sampleName/$sampleName.bed -f BED -g mm \
#--outdir /home/casual/ATAC_seq_out/mowei/$sampleName/peak_calling -n $sampleName --shift -75 --extsize 150 --nomodel \
#--call-summits --nolambda --keep-dup all -p 0.01

# bw for ucsc
bsub -q normal_1week -n 1 -e bamc.err -o bamc.log -J bamc \
bamCompare -b1 04.Align/BG-h3k9me3-ChIP-CT_TKD180602102.sorted.bam -b2 04.Align/BG-input-MT_TKD180602101.sorted.bam \
--operation subtract -o 04.Align/subtract_CT_input.bw

bamCompare -b1 04.Align/BG-h3k9me3-ChIP-MT_TKD180602103.sorted.bam -b2 04.Align/BG-input-MT_TKD180602101.sorted.bam \
--operation subtract -o 04.Align/subtract_MT_Input.bw

bamCompare -b1 04.Align/BG-h3k9me3-ChIP-MT_TKD180602103.sorted.bam -b2 04.Align/BG-h3k9me3-ChIP-CT_TKD180602102.sorted.bam \
--operation subtract -o 04.Align/subtract_MT_CT.bw

# motif
cat name2 |while read id;do
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' 06.heatmap/diff_MT_CT_peaks.narrowPeak >06.heatmap/diff_MT_CT_peaks_homerPeaks.bed


conda activate py36
cat name2 |while read id;do
bsub -q normal_1day_new -n 20 -e mot.err -o mot.log -J mot \
findMotifsGenome.pl 06.heatmap/diff_MT_CT_peaks_homerPeaks.bed mm10 07.motif/diff_MT_CT/ -size 200 -mask -p 20
##homer annotation

#findMotifsGenome.pl $id\_homerPeaks.bed mm10 ${sample}_motifDir -len 8,10,12
cd 06.heatmap
cat ../name2 |while read id;do
annotatePeaks.pl diff_MT_CT_peaks_homerPeaks.bed mm10 1>diff_MT_CT_peaks_homerPeaks.peakAnn.xls 2>diff_MT_CT_peaks_homerPeaks.annLog.txt 



source activate py36
## heatmap
bsub -q normal -n 20 -e gz.err -o gz.log -J gz \
computeMatrix reference-point -p 20 --referencePoint center -b 5000 -a 5000 -R 06.heatmap/MT_Input_summits.bed \
-S 04.Align/subtract_CT_input.bw 04.Align/subtract_MT_Input.bw -o 06.heatmap/MT_CT_specialMT.mat.gz --skipZeros 

bsub -q normal -n 20 -e gz.err -o gz.log -J gz \
computeMatrix reference-point -p 20 --referencePoint center -b 5000 -a 5000 -R 06.heatmap/CT_input_summits.bed \
-S 04.Align/subtract_CT_input.bw 04.Align/subtract_MT_Input.bw -o 06.heatmap/MT_CT_specialCT.mat.gz --skipZeros 


bsub -q normal_1day_new -n 1 -e K9.err -o K9.log -J K9 \
plotHeatmap --missingDataColor white --heatmapHeight 20 \
-m  06.heatmap/MT_CT_specialMT.mat.gz --colorMap vlag --heatmapHeight 10 \
-out 06.heatmap/MT_CT_specialMT.pdf 
#--outFileSortedRegions 07.heatmap/K9_O_I_input_iolsite_center_Heatmap.bed

bsub -q normal_1day_new -n 1 -e K9.err -o K9.log -J K9 \
plotHeatmap --missingDataColor white --heatmapHeight 20 \
-m  06.heatmap/MT_CT_specialCT.mat.gz --colorMap vlag --heatmapHeight 10 \
-out 06.heatmap/MT_CT_specialCT.pdf 





## 用mt降低的位点做热图
source activate python3

#预测样本的fragment size
bsub -q normal_1day -n 1 -e CT.err -o CT.log -J CT \
macs2 predictd -i 04.Align/BG-h3k9me3-ChIP-CT_TKD180602102.sorted.bam
#predicted fragment length is 241 bps
bsub -q normal_1day -n 1 -e MT.err -o MT.log -J MT \
macs2 predictd -i 04.Align/BG-h3k9me3-ChIP-MT_TKD180602103.sorted.bam
#predicted fragment length is 233 bps

# 大多数的转录因子chip_seq数据，推荐值为200，大部分组蛋白修饰的chip_seq数据，推荐值为147，也可以取两个平均值，两个值要一样
# condition1
bsub -q normal -n 1 -e CT.err -o CT.log -J MT-K9 \
macs2 callpeak -B -t 04.Align/BG-h3k9me3-ChIP-CT_TKD180602102.sorted.bam -c 04.Align/BG-input-MT_TKD180602101.sorted.bam \
-n CT_input --outdir 08.diff_heatmap -g mm --nomodel --extsize 237 -p 0.01
# condition2
bsub -q normal -n 1 -e WT-K9.err -o WT-K9.log -J WT-K9 \
macs2 callpeak -B -t 04.Align/BG-h3k9me3-ChIP-MT_TKD180602103.sorted.bam -c 04.Align/BG-input-MT_TKD180602101.sorted.bam \
-n MT_input --outdir 08.diff_heatmap -g mm --nomodel --extsize 237 -p 0.01
# 在运行这一步的时候，会输出每个样本过滤之后的reads数目，用在下一步的d1和d2
#CT:tags after filtering in treatment: 7481381
#MT:tags after filtering in treatment: 6057126

# 获取差异位点 1倍
bsub -q normal -n 1 -e bdgdiff.err -o bdgdiff.log -J bdgdiff \
macs2 bdgdiff --t1 08.diff_heatmap/MT_input_treat_pileup.bdg --c1 08.diff_heatmap/MT_input_control_lambda.bdg \
--t2 08.diff_heatmap/CT_input_treat_pileup.bdg --c2 08.diff_heatmap/CT_input_control_lambda.bdg \
-C 0.5 --d1 6057126 --d2 7481381 -g 100 -l 200 --o-prefix diff_MT_vs_CT --outdir 08.diff_heatmap/
#其中-d1和-d2的值就是第二步运行时输出的reads数目，-o参数指定输出文件的前缀。运行成功后，会产生3个文件
#分别为diff_c1_vs_c2_c3.0_cond1.bed 、iff_c1_vs_c2_c3.0_cond2.bed、diff_c1_vs_c2_c3.0_common.bed
#其中, con1.bed保存了在condition1中上调的peak, con2.bed保存了在condition2中上调的peak, common.bed文件中保存的是没有达到阈值的，非显著差异peak。
#上述3个文件格式是完全相同的，最后一列的内容为log10 likehood ratio值，用来衡量两个条件之间的差异，默认阈值为3，大于阈值的peak为组间差异显著的peak, 这个阈值可以通过-c参数进行调整。

source activate py36
##下调的做位点，也就是con2的bed文件
bsub -q normal -n 20 -e gz.err -o gz.log -J gz \
computeMatrix reference-point -p 20 --referencePoint center -b 5000 -a 5000 -R 08.diff_heatmap/diff_MT_vs_CT_c0.5_cond2.bed \
-S 04.Align/subtract_MT_Input.bw 04.Align/subtract_CT_input.bw \
-o 08.diff_heatmap/MT_CT_down_center_2.mat.gz
#--skipZeros 

bsub -q normal_1day -n 1 -e K9.err -o K9.log -J K9 \
plotHeatmap --heatmapHeight 10 --colorMap vlag \
-m 08.diff_heatmap/MT_CT_down_center_2.mat.gz --legendLocation upper-right \
-out 08.diff_heatmap/MT_CT_down_center_2.pdf 
#--outFileSortedRegions 07.heatmap/K9_O_I_input_iolsite_center_Heatmap.bed


## motif
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' 08.diff_heatmap/diff_MT_vs_CT_c0.5_cond2.bed >08.diff_heatmap/diff_MT_vs_CT_c0.5_cond2_homerPeaks.bed

bsub -q normal_1day_new -n 20 -e mot.err -o mot.log -J mot \
findMotifsGenome.pl 08.diff_heatmap/diff_MT_vs_CT_c0.5_cond2_homerPeaks.bed mm10 07.motif/diff_MT_vs_CT_c0.5_cond2/ -size 200 -mask -p 20
##homer annotation


