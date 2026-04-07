# Split Single-Cell Pipeline

按顺序运行下面 5 个脚本即可：

```bash
cd /Users/hongfan/Desktop/CWMDA/SingleCell_transcriptom
Rscript split_pipeline/01_qc_and_merge.R
Rscript split_pipeline/02_global_clustering.R
Rscript split_pipeline/03_find_markers.R
Rscript split_pipeline/04_cell_annotation.R
Rscript split_pipeline/05_cell_proportion_and_roe.R
```

每一步都会读取上一步的 `RDS`，并在 `04_rds` 下输出新的中间结果：

- `01.sc_after_QC_merged.rds`
- `02.sc_global_clustered.rds`
- `03.all_markers.rds`
- `04.sc_annotated.rds`
- `05.sc_final.rds`

为了兼容你原来依赖旧文件名的下游脚本，`04_cell_annotation.R` 还会额外保存一份：

- `04_rds/01.sc_annoted.rds`

主要输出目录：

- `01_QC`
- `02_GlobalClustering`
- `03_FindMarkers`
- `04_CellAnnotation`
- `05_CellProportion`
