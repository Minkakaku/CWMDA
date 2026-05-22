# bin/bash
# find . -name '*_2.fq.gz' -exec basename {} _2.fq.gz \; > SRR_Acc_List.txt &
# mkdir 00.rawdata 01.qc 02.bam 03.counts
# prefetch --option-file SRR_Acc_List.txt -O 00.rawdata

# conda activate py39
# job_path
cd /data/gpfs02/wmo/work/HF/

for i in `cat SRR_Acc_List.txt` ; do
path = /data/gpfs02/wmo/work/HF/
reference_stat_file = /data/gpfs02/wmo/genome/mm39/mm39.fa
export reference_stat_file
export path
export i
bsub -q normal_1day -n 24 -e 00.err_log/test.err -o 00.err_log/test.log -J 1207test \
bash until_samtools.sh
done

treat1 = 
treat2 = 
export treat1
export treat2
bash picard_merge.sh



cd 03.counts
Rscript 01.scripts/tools/hisat2_merge.r