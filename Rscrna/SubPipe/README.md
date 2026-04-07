# Split Generic Subcluster Pipeline

这套脚本是通用单亚群模板，不是只给 Epithelial 用的。

先在 [`00_subcluster_common.R`](/Users/hongfan/Desktop/CWMDA/SingleCell_transcriptom/split_subcluster/00_subcluster_common.R) 里改：

- `PARAM$target_celltype`
- `PARAM$input`
- 如果需要的话，`PARAM$use_harmony`、`dims`、`resolution`

然后按顺序运行：

```bash
cd /Users/hongfan/Desktop/CWMDA/SingleCell_transcriptom
Rscript split_subcluster/01_subcluster_cluster.R
Rscript split_subcluster/02_subcluster_find_markers.R
Rscript split_subcluster/03_subcluster_annotation.R
Rscript split_subcluster/04_subcluster_proportion_and_roe.R
```

推荐使用顺序：

1. 先跑 `01_subcluster_cluster.R`，看 `03_Subcluster/<目标主群>/UMAP_clusters.pdf`
2. 再跑 `02_subcluster_find_markers.R`，看 `markers_all_clusters.csv`
3. 然后编辑 `03_subcluster_annotation.R` 里的 `ANNOTATION`
4. 最后跑 `04_subcluster_proportion_and_roe.R`

核心输出：

- `03_Subcluster/<目标主群>/subcluster_clustered.rds`
- `03_Subcluster/<目标主群>/markers_all_clusters.csv`
- `03_Subcluster/<目标主群>/subcluster_annotated.rds`
- `03_Subcluster/<目标主群>/sample_by_sub_celltype_count.csv`
- `03_Subcluster/<目标主群>/ROE_matrix.csv`
