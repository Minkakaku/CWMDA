docker run --name monocle_wzb --rm \
  -v /data/NAS01/wzb/bak/home/:/data/NAS01/wzb/bak/home/ \
  -v /home/wzb:/home/wzb \
  -v /lvdata/wzb:/lvdata/wzb \
  -v /data/NAS03/backup/wzb/wzb/:/data/NAS03/backup/wzb/wzb/ \
  wangzhenbo/monocle:20240909 /bin/bash -c "
  su - wzb -c '
    # 运行 R 脚本
    Rscript /lvdata/wzb/pipline/function/Pseudotime/pseudotime.R \
      -r \"/lvdata/wzb/scRNA/FWSC20240705-08/2024.11.15/Interneurons/h5_trans.rds\" \
      -l 8 \
      -o \"/lvdata/wzb/scRNA/FWSC20240705-08/2024.11.15/Interneurons/\"
  '"
# -p 0.3 \
# -M 1500 \
# -R \"TRUE\" \