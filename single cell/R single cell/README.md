# R Single Cell Workflow

这里放 R 版本的单细胞主流程代码。主流程的参数集中在 `00_config.R`，包装函数集中在 `00_functions.R`，步骤脚本只作为入口。

## 目录

- `00_config.R`：所有可调参数，包括 QC 阈值、Harmony/CCA、marker 阈值、注释、污染 cluster 去除和比例图设置。
- `00_functions.R`：所有包装函数和步骤函数。
- `main/`：主群分析步骤。
- `subcluster/`：亚群分析步骤。
- `tools/sample_discovery.R`：R 版样本自动发现和自动分组工具。
- `run_main.sh`：主群流程入口。
- `run_subcluster.sh`：亚群流程入口。

## 主群流程

推荐先跑 QC1，检查图后修改 `00_config.R`，再继续过滤和后续分析：

```bash
bash "single cell/R single cell/run_main.sh" qc1
bash "single cell/R single cell/run_main.sh" after_qc
```

也可以按单步运行：

```bash
bash "single cell/R single cell/run_main.sh" preflight
bash "single cell/R single cell/run_main.sh" qc2
bash "single cell/R single cell/run_main.sh" cluster
bash "single cell/R single cell/run_main.sh" markers
bash "single cell/R single cell/run_main.sh" annotation
bash "single cell/R single cell/run_main.sh" proportion
```

## 亚群流程

先在 `00_config.R` 设置 `SUBCLUSTER$target_celltype`、分辨率、去批次方法、注释和污染 cluster 去除参数，然后运行：

```bash
bash "single cell/R single cell/run_subcluster.sh" all
```

也可以单步运行：

```bash
bash "single cell/R single cell/run_subcluster.sh" cluster
bash "single cell/R single cell/run_subcluster.sh" markers
bash "single cell/R single cell/run_subcluster.sh" annotation
bash "single cell/R single cell/run_subcluster.sh" proportion
```

## R 版样本自动识别

这个工具从 Python 的 `sample_discovery.py` 迁移而来，命令保持同样的 `discover / suggest / assign` 思路。

```bash
Rscript "single cell/R single cell/tools/sample_discovery.R" discover \
  --data-root data2analysis \
  --output-manifest sample_manifest.tsv

Rscript "single cell/R single cell/tools/sample_discovery.R" suggest \
  --manifest sample_manifest.tsv \
  --output-report group_suggestions.tsv

Rscript "single cell/R single cell/tools/sample_discovery.R" assign \
  --manifest sample_manifest.tsv \
  --group-count 4 \
  --output-sample-sheet sample_sheet.tsv \
  --output-report group_suggestions.tsv
```

输出的 `sample_sheet.tsv` 包含 `sample/group/sample_dir/outs_dir/matrix_dir/matrix_type`，可以作为后续流程的样本表。
