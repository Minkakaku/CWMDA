#!/bin/bash
docker run --name SCENIC_wzb --rm \
-v /home/wzb:/home/wzb \
-v /lvdata/wzb:/lvdata/wzb \
-v /data/NAS03/backup/wzb:/data/NAS03/backup/wzb \
-u $(id -u):$(id -g) \
aertslab/pyscenic:0.12.1 pyscenic grn --num_workers 32 \
--sparse \
--method grnboost2 \
--output $2 \
$1 \
/lvdata/wzb/pipline/function/SCENIC/mm_mgi_tfs.txt

docker run --name SCENIC_wzb --rm \
-v /home/wzb:/home/wzb \
-v /lvdata/wzb:/lvdata/wzb \
-v /data/NAS03/backup/wzb:/data/NAS03/backup/wzb \
-u $(id -u):$(id -g) \
aertslab/pyscenic:0.12.1 pyscenic ctx --num_workers 32 \
--output $3 \
--expression_mtx_fname $1 \
--mode "custom_multiprocessing" \
--annotations_fname /lvdata/wzb/pipline/function/SCENIC/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl \
$2 \
/lvdata/wzb/pipline/function/SCENIC/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

docker run --name SCENIC_wzb --rm \
-v /home/wzb:/home/wzb \
-v /lvdata/wzb:/lvdata/wzb \
-v /data/NAS03/backup/wzb:/data/NAS03/backup/wzb \
-u $(id -u):$(id -g) \
aertslab/pyscenic:0.12.1 pyscenic aucell --num_workers 32 \
--output $4	\
$1	\
$3