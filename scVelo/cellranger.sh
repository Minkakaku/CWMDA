directory="/home/wzb/scRNA/FWSC20240705-08"
for dir in "$directory"/*/; do
    nohup /lvdata/wzb/pipline/function/scvelo/cellranger_run_gex.sh "$dir" /home/ztt/database/refdata-gex-mm10-2020-A /home/wzb/scRNA/output/FWSC20240705-08 &
done