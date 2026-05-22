dir="/home/wzb/scRNA/FW2024-225_03/iMCD_12_after"

nohup velocyto run10x -m /home/wzb/velocyto_gtf/GRCh38_rmsk.gtf \
"$dir" \
/data/NAS01/wzb/reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf &
