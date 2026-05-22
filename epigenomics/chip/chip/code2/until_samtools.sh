# transfer sra2fastq
fasterq-dump -e 24 00.rawdata/$i/$i.sra -O 00.rawdata
# qc with gzip, final name fq.gz
trim_galore  -q 30 --phred33 -J 24 \
            --paired 00.rawdata/$i\_1.fastq 00.rawdata/$i\_2.fastq \
            --gzip -o 01.qc
# align with hisat2
hisat2 -t -p 24 -x /data/gpfs02/wmo/genome/hisat2_index/hisat2_mm39/mm39 \
-1 01.qc/$i\_1_val_1.fq.gz -2 01.qc/$i\_2_val_2.fq.gz \
-S 02.bam/$i.sam
# post align qc
samtools view -@ 24 -b 02.bam/$i.sam -o 02.bam/$i.bam \
samtools sort -@ 24 -o 02.bam/$i.sorted.bam 02.bam/$i.bam \
samtools index -@ 24 02.bam/$i.sorted.bam
samtools stats @ 24 --reference $reference_stat_file \
02.bam/$i.sorted.bam;
samtools idxstats 02.bam/$i.sorted.bam
samtools flagstat --threads 24 02.bam/$i.sorted.bam
