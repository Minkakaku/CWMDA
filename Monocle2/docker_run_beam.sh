docker run --name monocle_beam_wzb --rm \
-v /data/NAS01/wzb/bak/home/:/data/NAS01/wzb/bak/home/ \
-v /home/wzb:/home/wzb \
-v /lvdata/wzb:/lvdata/wzb \
-v /data/NAS03/backup/wzb:/data/NAS03/backup/wzb \
wangzhenbo/monocle:20240909 /bin/bash -c "
su - wzb -c '
Rscript /lvdata/wzb/pipline/function/Pseudotime/Pseudo_BEAM.R \
-r \"/lvdata/wzb/scRNA/FWSC20240471/2024.11.12/Pseudotime/Pseudotime.monocle.rds\" \
-b 1 \
-T 250 \
-q 0.0001 \
-o \"/lvdata/wzb/scRNA/FWSC20240471/2024.11.12/BEAM\"
'"
