directory="/home/wzb/scRNA/FW2023-504"
for dir in "$directory"/*/; do
    nohup velocyto run10x -m /home/wzb/velocyto_gtf/GRCh38_rmsk.gtf \
    "$dir" \
    /data/NAS01/wzb/reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf &
done