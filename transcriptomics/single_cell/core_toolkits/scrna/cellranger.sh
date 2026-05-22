#下载原始文件
cat SRR_Acc_List |while read id;
do prefetch $i -O `pwd` && echo "** ${i}.sra done **"
done

#SRR文件转换为FQ
cat SRR_Acc_List |while read id;do
 bsub -q normal_1day -n 10 -e $id.err -o $id.log -J $id   \
fastq-dump 00.rawdata/$id/$id.sra --split-files --gzip -O fastq
done


#FQ文件改名，方便后续操作
cat SRR_Acc_List| while read i ;do
mv ${i}_1*.gz ${i}_S1_L001_R1_001.fastq.gz; mv ${i}_2*.gz ${i}_S1_L001_R2_001.fastq.gz; mv ${i}_3*.gz ${i}_S1_L001_I1_001.fastq.gz;
done

#下载小鼠cellranger参考组文件
wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-mm10-3.0.0.tar.gz
tar -xzvf refdata-cellranger-mm10-3.0.0.tar.gz 

#将BCL转化为FQ文件（一般不用）
cellranger mkfastq --id=bcl \
                    --run=/path/to/bcl \
                    --samplesheet=samplesheet-1.2.0.csv


#进行FQMAP和定量(高版本用这个)
cat SRR_Acc_List |while read id;
do
bsub -q normal_1day -n 20 -e $id.err -o $id.log -J $id \
cellranger count --id=${id}_SCcounts   \
--fastqs=/data/gpfs02/wmo/work/HF/jobs/YXX/fastq  \
--transcriptome=/data/gpfs02/wmo/genome/single_cell/refdata-gex-GRCh38-and-mm10-2020-A   \
--sample=${id} \
--nosecondary 
done

#对count结果进行整合（cellranger aggr命令）
#首先创建需要整合的样本以及其h5文件路径（csv格式）
cat >PH0_libraries.csv
sample_id,molecule_h5
PH0H1,/data/gpfs02/wmo/work/HF/jobs/YXX/fastq/SRR12755872_SCcounts/outs/molecule_info.h5
PH0H2,/data/gpfs02/wmo/work/HF/jobs/YXX/fastq/SRR12755873_SCcounts/outs/molecule_info.h5
PH0H3,/data/gpfs02/wmo/work/HF/jobs/YXX/fastq/SRR12755874_SCcounts/outs/molecule_info.h5
PH0H4,/data/gpfs02/wmo/work/HF/jobs/YXX/fastq/SRR12755875_SCcounts/outs/molecule_info.h5
PH0H5,/data/gpfs02/wmo/work/HF/jobs/YXX/fastq/SRR12755876_SCcounts/outs/molecule_info.h5
PH0H6,/data/gpfs02/wmo/work/HF/jobs/YXX/fastq/SRR12755877_SCcounts/outs/molecule_info.h5
PH0H7,/data/gpfs02/wmo/work/HF/jobs/YXX/fastq/SRR12755878_SCcounts/outs/molecule_info.h5
PH0H8,/data/gpfs02/wmo/work/HF/jobs/YXX/fastq/SRR12755879_SCcounts/outs/molecule_info.h5


cat >PH48h_libraries.csv
sample_id,molecule_h5
PH48H1,/data/gpfs02/wmo/work/HF/jobs/YXX/fastq/SRR12755880_SCcounts/outs/molecule_info.h5
PH48H2,/data/gpfs02/wmo/work/HF/jobs/YXX/fastq/SRR12755881_SCcounts/outs/molecule_info.h5
PH48H3,/data/gpfs02/wmo/work/HF/jobs/YXX/fastq/SRR12755882_SCcounts/outs/molecule_info.h5
PH48H4,/data/gpfs02/wmo/work/HF/jobs/YXX/fastq/SRR12755883_SCcounts/outs/molecule_info.h5
PH48H5,/data/gpfs02/wmo/work/HF/jobs/YXX/fastq/SRR12755884_SCcounts/outs/molecule_info.h5
PH48H6,/data/gpfs02/wmo/work/HF/jobs/YXX/fastq/SRR12755885_SCcounts/outs/molecule_info.h5
PH48H7,/data/gpfs02/wmo/work/HF/jobs/YXX/fastq/SRR12755886_SCcounts/outs/molecule_info.h5
PH48H8,/data/gpfs02/wmo/work/HF/jobs/YXX/fastq/SRR12755887_SCcounts/outs/molecule_info.h5

# 其中molecule_h5：文件molecule_info.h5 file的路径
# 下一步接入R包seurat流程，此处标准化选none
bsub -q normal_1week -n 20 -e PH.err -o PH.log -J PH \
cellranger aggr --id=PH48h_lib \
                --csv=PH48h_libraries.csv \
                --normalize=none \
                --nosecondary

cd /data/gpfs02/wmo/work/HF/scScripts
bsub -q normal_1week -n 30 -e hor.err -o hor.log -J hor \
Rscript 