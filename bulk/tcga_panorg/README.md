# TCGA PANoRG Analysis

这是整理进 `bulk` 目录后的 TCGA bulk 分析子模块，主要完成：

- 下载 TCGA 表达和临床数据
- 进行用户风格的预处理
- 基于 PANoRG 基因集计算 `PANoScore`
- 做聚类、PCA、生存分析、Cox 和临床关联分析

## 目录结构

- `config`
  - `panorgs.txt`
  - `tcga_cohorts.csv`
  - `molecular_entity_candidates.csv`
- `R`
  - `tcga_user_style_utils.R`
  - `tcga_user_style_prepare.R`
  - `tcga_user_style_panorg.R`
- `scripts`
  - `00_install_packages.R`
  - `01_download_tcga_data.R`
  - `02_run_panorg_analysis.R`
  - `03_run_all.R`

## 运行方式

```bash
Rscript bulk/tcga_panorg/scripts/00_install_packages.R
Rscript bulk/tcga_panorg/scripts/01_download_tcga_data.R
Rscript bulk/tcga_panorg/scripts/02_run_panorg_analysis.R
```

也可以直接一步跑完：

```bash
Rscript bulk/tcga_panorg/scripts/03_run_all.R
```

## 输出

结果会写到 `bulk/tcga_panorg/results/<cohort_id>/`，包括：

- `tables/sample_annotations.csv`
- `figures/PANoRG_heatmap.pdf`
- `figures/PCA_cluster.pdf`
- `figures/KM_OS_cluster.pdf`
- `figures/KM_OS_PANoScore.pdf`
- `tables/clinical_associations.csv`
- `tables/cox_OS_univariate.csv`
- `tables/cox_OS_multivariate.csv`

## 说明

- 现在脚本入口已经改成按脚本位置自动定位项目根目录，不要求你必须先 `cd` 到模块目录。
- 原仓库里的 `TCGA_analysis` 目前还保留着，作为迁移前的历史版本。
