# TCGA PANoRG Analysis

这个目录用于复现类似论文 Fig.3 的 PANoptosis 相关分析流程，当前默认支持：

- `TCGA_CRC`：合并 `TCGA-COAD` + `TCGA-READ`
- `TCGA_PRAD`：`TCGA-PRAD`

核心步骤包括：

- 基于 PANoRGs 做一致性聚类，得到 `C1/C2`
- 基于同一基因集计算 `ssGSEA` 分数 `PANoScore`
- 绘制 `PCA`
- 做 `KM` 生存分析和 `Cox` 回归
- 评估 `PANoScore` 与临床特征的关联
- 比较不同分子实体中的 `PANoScore` 分布

## 目录结构

- `config/panorgs.txt`：PANoRG 基因集
- `config/tcga_cohorts.csv`：CRC / PRAD 队列定义
- `config/custom_molecular_entities_template.csv`：自定义分子实体模板
- `R/utils.R`：通用辅助函数
- `R/tcga_download.R`：TCGA 下载与整理
- `R/panorg_analysis.R`：聚类、打分、生存与作图
- `scripts/00_install_packages.R`：安装依赖包
- `scripts/01_download_tcga_data.R`：下载并导出表达矩阵与临床表
- `scripts/02_run_panorg_analysis.R`：运行主分析

## 推荐运行顺序

在 `E:\CWMDA\TCGA_analysis` 下运行：

```r
Rscript scripts/00_install_packages.R
Rscript scripts/01_download_tcga_data.R
Rscript scripts/02_run_panorg_analysis.R
```

## 产出文件

下载步骤会生成：

- `data_raw/TCGA_CRC_expr.tsv.gz`
- `data_raw/TCGA_CRC_clin.tsv.gz`
- `data_raw/TCGA_PRAD_expr.tsv.gz`
- `data_raw/TCGA_PRAD_clin.tsv.gz`

主分析会在 `results/TCGA_CRC` 和 `results/TCGA_PRAD` 下输出：

- `sample_annotations.tsv.gz`
- `panorg_heatmap.png`
- `pca_clusters.png`
- `km_cluster.pdf`
- `km_panoscore.pdf`
- `cox_univariate.tsv`
- `cox_multivariate.tsv`
- `clinical_associations.tsv`
- `PANoScore_by_<clinical_feature>.png`
- `PANoScore_by_molecular_entity.png`

## 分子实体说明

原文中“5类分子实体”的定义来自对应癌种的分子分型背景。对 CRC / PRAD，这里优先尝试使用 `TCGAbiolinks::PanCancerAtlas_subtypes()` 或 `TCGAquery_subtype()` 的注释，并统一映射到临床表中的 `molecular_entity` 列。

如果你希望严格按你自己的“五分类”来比较 `PANoScore`，可以把模板复制为：

```text
config/custom_molecular_entities.csv
```

然后填入三列：

- `disease`
- `patient_id`
- `molecular_entity`

脚本会优先使用你提供的 `molecular_entity` 覆盖默认注释。

## 方法假设

- 表达数据默认使用 `HTSeq - FPKM-UQ`
- 主分析默认对表达矩阵做 `log2(x + 1)` 转换
- 聚类使用 `ConsensusClusterPlus`，默认固定读取 `k = 2`
- `PANoScore` 用同一 PANoRG 基因集做 `ssGSEA`
- `PANoScore` 高低组默认按中位数切分
- `PRAD` 的 `OS` 事件较少，若你后续想改为 `PFI/DFS`，可以在 `config/tcga_cohorts.csv` 和临床整理函数里调整终点定义
