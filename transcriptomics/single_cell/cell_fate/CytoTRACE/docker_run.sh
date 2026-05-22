docker run --name cytotrac_wzb --rm \
-v /data/NAS01/wzb/bak/home/:/data/NAS01/wzb/bak/home/ \
-v /home/wzb:/home/wzb \
-v /lvdata/wzb:/lvdata/wzb \
-v /data/NAS03/backup/wzb/wzb:/data/NAS03/backup/wzb/wzb \
-u $(id -u):$(id -g) \
wangzhenbo/seurat_pip:v1.1 Rscript /lvdata/wzb/pipline/function/CytoTRACE/cytotrace.R \
-r "/lvdata/wzb/scRNA/FWSC20241703/2024.12.02/h5_trans.rds" \
-o /lvdata/wzb/scRNA/FWSC20241703/2024.12.02/Topic5 -n 32