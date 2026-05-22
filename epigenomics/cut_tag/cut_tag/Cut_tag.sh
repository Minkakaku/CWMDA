#!/bin/bash
#PBS -l nodes=1:ppn=20,walltime=100:00:00
#PBS -j oe
#PBS -N name

#cut-tag peak calling progress(LZX)

cd /data/gpfs01/wmo/work/DNA-seq/CUT_TAG/Lzx_seh1_nup133
cat name |while read id
do
#conda activate py36
bsub -q normal -n 15 -e $id.err -o $id.log -J $id cutadapt -j 15 -a GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -A CGTCGGCAGCGTCAGATGTGTATAAGAGACAG -q 20 \
-o 02.cleandata/$id/$id\_1\_trim.fq.gz -p 02.cleandata/$id/$id\_2\_trim.fq.gz 02.cleandata/$id/$id\_1.clean.fq.gz 02.cleandata/$id/$id\_2.clean.fq.gz
#ADAPTERS = ['CTGTCTCTTATACACATCT', 'AGATGTGTATAAGAGACAG', 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG', 'CTGTCTCTTATACACATCTGACGCTGCCGACGA', 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG', 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC',]
#5’-AATGATACGGCGACCACCGAGATCTACACIIIIIIIITCGTCGGCAGCGTCAGATGTGTATAAGAGACAG-NNNNNN-CTGTCTCTTATACACATCTCCGAGCCCACGAGACIIIIIIIIATCTCGTATGCCGTCTTCTGCTTG-3

#conda deactivate
#conda activate python3

cat name |while read id
do
bsub -q normal -n 1 -e $id.err -o $id.log -J $id fastqc 02.cleandata/$id/$id\_1\_trim.fq.gz
bsub -q normal -n 1 -e $id.err -o $id.log -J $id fastqc 02.cleandata/$id/$id\_2\_trim.fq.gz

#bowtie2-build genome/mm10/mm10.fa Index/bowtie2_mm10/mm10

cat name |while read id
do
bsub -q normal_1day -n 10 -e $id.err -o $id.log -J $id bowtie2 -p 10 --very-sensitive --phred33 -X 2000 -x ~/Index/bowtie2_mm10/mm10 \
-1 02.cleandata/$id/$id\_1_trim.fq.gz -2 02.cleandata/$id/$id\_2_trim.fq.gz -S 03.Align/$id.sam
bsub -q normal_1day -n 2 -e $id.err -o $id.log -J $id samtools view -bS -@ 2 03.Align/$id.sam -o 03.Align/$id.bam
bsub -q normal_1day -n 2 -e $id.err -o $id.log -J $id samtools sort -@ 2 -O bam -o 03.Align/$id.sorted.bam 03.Align/$id.bam
done

rm 03.Align/$id.sam
rm 03.Align/$id.bam


## remove duplicate reads
bsub -q normal -n 1 -e $id.err -o $id.log -J $id java -jar ~/script/picard.jar \
MarkDuplicates I=03.Align/$id.sorted.bam O=03.Align/$id.sorted.dedup.bam \
M=03.Align/$id\_picard_metrics_bam.txt REMOVE_DUPLICATES=true
## build index
#bsub -q normal_1day -n 1 -e $id.err -o $id.log -J $id java -jar ~/script/picard.jar BuildBamIndex I=03.Align/$id.sorted.dedup.bam


cat name |while read id
do
bsub -q normal -n 1 -e $id.err -o $id.log -J $id 

samtools view -h 03.Align/$id.sorted.dedup.bam -o 03.Align/$id.sorted.dedup.sam
python ~/script/removeChrom.py 03.Align/$id.sorted.dedup.sam 03.Align/$id.sorted.dedup.rmMT.sam chrM ##
samtools view -b 03.Align/$id.sorted.dedup.rmMT.sam -o 03.Align/$id.sorted.dedup.rmMT.bam
## build index for IGV
java -jar ~/script/picard.jar BuildBamIndex I=03.Align/$id.sorted.dedup.rmMT.bam
#samtools view -h /home/casual/ATAC_seq_out/mowei/$sampleName/$sampleName.sorted.dedup.rmMT.bam 
python ~/script/SAMtoBED.py -x -v -i 03.Align/$id.sorted.dedup.rmMT.sam -o 03.Align/$id.bed

rm 03.Align/*.sorted.dedup.sam 03.Align/*.sorted.dedup.rmMT.sam 


#macs peak calling
conda deactivate
conda activate python3 ##python版本是2.7,必须是2.7

cat name |while read id
do
bsub -q normal -n 1 -e $id.err -o $id.log -J $id 
macs2 callpeak -t 03.Align/$id.bed -n $id --outdir 04.macs/ -g mm --nomodel --nolambda --keep-dup all --call-summits -q 0.01

conda activate py36
##motif
bsub -q normal -n 1 -e $id.err -o $id.log -J $id 
findMotifsGenome.pl 04.macs/$id\_summits.bed mm10 05.motif/$id/ -size 200 -mask -p 10

##homer annotation
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' 04.macs/$id\_summits.bed >04.macs/$id\_homerPeaks.bed
#findMotifsGenome.pl $id\_homerPeaks.bed mm10 ${sample}_motifDir -len 8,10,12
cd 04.macs/
annotatePeaks.pl $id\_homerPeaks.bed mm10 1>$id.peakAnn.xls 2>$id.annLog.txt 
done 

scp -r wmo@hpc.xmu.edu.cn:/data/gpfs01/wmo/work/DNA-seq/CUT_TAG/Lzx_seh1_nup133_v2 /data/LZX/测序文件/cut_tag_第二批









