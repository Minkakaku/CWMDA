docker run --name monocle_plot_re --rm \
-v /data/NAS01/wzb/bak/home/:/data/NAS01/wzb/bak/home/ \
-v /home/wzb:/home/wzb \
-v /lvdata/wzb:/lvdata/wzb \
-v /data/NAS03/backup/wzb:/data/NAS03/backup/wzb \
wangzhenbo/monocle:20240909 /bin/bash -c "
su - wzb -c '
Rscript /lvdata/wzb/pipline/function/Pseudotime/PseudotimeRdsplot.R \
-r \"/lvdata/wzb/scRNA/FW2024-225_05/Myeloid/Pseudotime/Pseudotime.monocle.rds\" \
-R \"TRUE\" \
-S 5 \
-o \"/lvdata/wzb/scRNA/FW2024-225_05/Myeloid/Pseudotime/Reorder\"
'"
