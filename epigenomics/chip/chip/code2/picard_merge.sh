# merge resequenced bam files
ls *.bam | grep $treat1 | picard MergeSamFiles -O $treat1.bam
ls *.bam | grep $treat2 | picard MergeSamFiles -O $treat2.bam

# Mark duplicates & filter BAM files after merging
for i in $treat1 $treat1;
do
picard MarkDuplicates \
--INPUT $i.bam \
--OUTPUT $i.MarkDup.bam \
--METRICS_FILE $i.MarkDuplicates.metrics.txt
done
# Filter BAM file with BamTools
# SE: '-F 0x004' PE: '-F 0x004 -F 0x0008 -f 0x001'
samtools view -F 0x004 -F 0x0008 -f 0x001 -F 0x0400 \
-q 1 -b $i.MarkDup.bam \
| bamtools filter -out $i.MarkDup.filtered.bam -script ./bamtools_filter_pe.json