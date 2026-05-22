## download               
prefetch --option-file SRR_Acc_List.txt -O 00.rawdata/

## translate srr file to fastq file
cat SRR_Acc_List.txt | while read id;do
bsub -q normal_1day_new -n 1 -e $id.err -o $id.log -J $id \
fastq-dump 00.rawdata/$id/$id.sra --gzip --split-3 -O 01.fastq/
done

## QC
cat H3K27me3d | while read id;do
bsub -q normal_1week -n 16 -e $id.err -o $id.log -J $id \
fastqc -o 02.QC -t 16 --noextract -q 01.fastq/$id\_1.fastq 01.fastq/$id\_2.fastq
done

cat H3K27me3s | while read id;do
bsub -q normal_1week -n 16 -e $id.err -o $id.log -J $id \
fastqc -o 02.QC -t 16 --noextract -q 01.fastq/$id.fastq
done

multiqc 02.QC

cat SRR_Acc_List.txt | while read id;do
bsub -q normal_1day -n 10 -e $id.err -o $id.log -J $id \
trim_galore -j 20 -q 20 --phred33 --stringency 3 -e 0.1 \
01.fastq/$id.fastq.gz \
--gzip --fastqc -o 02.QC
done

cat SRR_Acc_List.txt| tail -n 6 | while read id;do
bsub -q normal_1day -n 20 -e $id.err -o $id.log -J $id \
trim_galore --cores 4 -q 20 --phred33 --stringency 3 --length 100 -e 0.1 \
--paired 01.fastq/$id\_1.fastq 01.fastq/$id\_2.fastq \
--gzip --fastqc -o 02.QC
done

# # build bowtie2 index
# bsub -q normal_1week -n 20 -e 1.err -o 1.log -J 1 \
# bowtie2-build --threads 20 /data/gpfs02/wmo/genome/mm39/mm39.fa /data/gpfs02/wmo/genome/bowtie2_index/mm39/mm39

# align
cat SRR_Acc_List.txt | while read id;do
bsub -q normal_1day_new -n 10 -e $id.err -o $id.log -J $id \
bowtie2 -p 10 --sensitive -X 500 -x /data/gpfs02/wmo/genome/bowtie2_index/mm39/mm39 -U 02.QC/$id\_trimmed.fq.gz \
-S 03.align/$id.sam
done

cat SRR_Acc_List.txt| tail -n 6 | while read id;do
bsub -q normal_1day_new -n 20 -e $id.err -o $id.log -J $id \
bowtie2 -p 20 --sensitive -X 500 -x /data/gpfs02/wmo/genome/bowtie2_index/GRCm39/GRCm39 \
-1 02.QC/$id\_1_val_1.fq.gz -2 02.QC/$id\_2_val_2.fq.gz \
-S 03.align/$id.sam
done

# convert to bam
cat SRR_Acc_List.txt | while read id;do
bsub -q normal_1day_new -n 1 -e $id.err -o $id.log -J $id \
samtools view -@ 1 -b -S 03.align/$id.sam -o 03.align/$id.bam
done

# sort
cat SRR_Acc_List.txt  | while read id ; do
bsub -q normal_1day_new -n 1 -e $id.err -o $id.log -J $id \
samtools sort -@ 1 -o 03.align/$id.sorted.bam 03.align/$id.bam
done

# index
cat SRR_Acc_List.txt | while read id ; do
bsub -q normal_1day_new -n 10 -e $id.err -o $id.log -J $id \
samtools index -@ 10 03.align/$id.sorted.bam
done

# compare
cat SRR_Acc_List.txt | while read id ; do
bsub -q normal_1day_new -n 10 -e $id.err -o $id.log -J $id \
bamCoverage --bam 03.align/$id.sorted.bam  -o 04.bw/$id.log2ratio.bw \
    --numberOfProcessors 10 \
    --binSize 5
done


# # remove PCR-duplicate
# cat SRR_Acc_List.txt | head -n 8|while read id;do
# bsub -q normal_1week -n 20 -e $id.err -o $id.log -J $id \
# picard MarkDuplicates \
# REMOVE_DUPLICATES=true \
# I=03.align/$id.sorted.bam \
# O=03.align/$id.sorted.markdup.bam \
# M=03.align/$id.markdup.txt
# done

# callpeak
### macs2 calls H3K27me3 broad peaks
cat H3K27me3d | while read i ;do
bsub -q normal_1day_new -n 30 -e $i.err -o $i.log -J $i \
macs2 callpeak -t 03.align/$i.sorted.bam \
-c 03.align/SRR9018522.sorted.bam \
-B \
-g mm -p 1e-3 --nomodel --shift 37 --extsize 73 -n $i\_H3K27me3 --outdir 04.results/$i
done


### macs2 calls H3K4me3 narrow peaks
### macs2 calls H3K9me3 narrow peaks
cat SRR_Acc_List.txt | while read i ;do
bsub -q normal -n 20 -e /data/gpfs02/wmo/work/HF/00.err_log/0926.err -o /data/gpfs02/wmo/work/HF/00.err_log/0926.log -J 0926 \
macs2 callpeak -t 03.align/$i.sorted.bam \
-g mm \
-q 0.05 \
-f BAM -B \
-n $i \
--outdir 04.results/$i
done

# visualize
cat SRR_Acc_List.txt | while read i ;do
bsub -q normal -n 20 -e /data/gpfs02/wmo/work/HF/00.err_log/0926.err -o /data/gpfs02/wmo/work/HF/00.err_log/0926.log -J 0926 \
sort -k1,1 -k2,2n 04.results/$i/$i\_treat_pileup.bdg -o 04.results/$i\_sort.bdg
done

cat SRR_Acc_List.txt  | while read i ;do
bsub -q normal -n 20 -e /data/gpfs02/wmo/work/HF/00.err_log/0926.err -o /data/gpfs02/wmo/work/HF/00.err_log/0926.log -J 0926 \
bedGraphToBigWig 04.results/$i\_sort.bdg /data/gpfs02/wmo/genome/mm39/mm39.chrom.sizes 04.results/$i\_sort.bw
done

## solve the problem that numbers about different size
bedClip MACS/$i/$i\_treat_pileup.bdg /data/gpfs02/wmo/genome/mm8/mm8.chrom.sizes MACS/$i/$i\_clean_pileup.bdg

## files
find -name 'WTh*' -exec basename {} .bam \; > name2