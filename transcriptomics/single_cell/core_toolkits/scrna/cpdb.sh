bsub -q normal -n 20 -e /data/gpfs02/wmo/work/HF/00.err_log/0923.err -o /data/gpfs02/wmo/work/HF/00.err_log/0923.log -J 0923 \
cellphonedb method statistical_analysis /data/gpfs02/wmo/work/HF/06.ST-seq/project/wk2/wk2_meta.txt /data/gpfs02/wmo/work/HF/06.ST-seq/project/wk2/wk2.h5ad --counts-data hgnc_symbol
