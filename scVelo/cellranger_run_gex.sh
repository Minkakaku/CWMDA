#!/bin/bash
#cellranger_run_gex.sh [fastq_dir] [reference_dir] [out_dir]
# --force-cells 9000 \
echo =======start cellranger at `date`==========

cd ${3}

ls ${1} | cut -d '_' -f 1 | uniq | while read i
do
    nohup /home/wzb/software/cellranger-7.1.0/bin/cellranger count \
    --id $i \
    --fastqs ${1} \
    --transcriptome ${2} \
    --include-introns true \
    --sample $i \
    --localcores 32 \
    --localmem 64 \
    &
done
wait
echo =========cellranger finish at `date`========
