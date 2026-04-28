# Bulk Codebase

这里是仓库内 bulk 相关代码的新统一目录，当前按两类内容整理：

- `tcga_panorg`
  - 一个相对完整的 bulk 主流程子模块
  - 负责 TCGA 下载、表达预处理、PANoRG 打分、聚类、生存和临床关联分析
- `r`
  - `cibersort_core.R` / `cibersort_pipeline.R`
  - `drug_score_fgsea.R`
  - `common_utils.R`

## 目录结构

- `r`
  - `common_utils.R`：CLI 参数、目录和表格读取工具
  - `cibersort_core.R`：CIBERSORT 核心函数
  - `cibersort_pipeline.R`：CIBERSORT 标准入口
  - `drug_score_fgsea.R`：药物 FGSEA 打分和组合评分入口
- `resources`
  - `LM22.txt`：CIBERSORT 常用签名矩阵
- `tcga_panorg`
  - `config`
  - `R`
  - `scripts`
  - `README.md`

## 典型入口

### 0. 安装 bulk 依赖

```bash
Rscript bulk/scripts/00_install_bulk_packages.R
```

### 1. CIBERSORT

```bash
Rscript bulk/r/cibersort_pipeline.R \
  --sig-matrix bulk/resources/LM22.txt \
  --mixture-file your_expression.tsv \
  --out-dir result/cibersort \
  --perm 100 \
  --qn false
```

### 2. Drug Score FGSEA

先准备一个 manifest，例如 `drug_manifest.tsv`：

```tsv
drug_name	dep_file
DrugA	path/to/drugA_dep.tsv
DrugB	path/to/drugB_dep.tsv
```

然后运行：

```bash
Rscript bulk/r/drug_score_fgsea.R \
  --manifest drug_manifest.tsv \
  --gmt-file pathways.gmt \
  --pathways ICD,Antigen,Glycolysis,Lipogenesis \
  --out-dir result/drug_score
```

### 3. TCGA PANoRG

```bash
Rscript bulk/tcga_panorg/scripts/00_install_packages.R
Rscript bulk/tcga_panorg/scripts/01_download_tcga_data.R
Rscript bulk/tcga_panorg/scripts/02_run_panorg_analysis.R
```

或者直接：

```bash
Rscript bulk/tcga_panorg/scripts/03_run_all.R
```

## 说明

- 这轮整理以 bulk 主流程和常用工具统一为主，还没有删除旧 bulk 目录。
- 旧目录目前仍保留，主要是为了避免打断你现有分析；后面如果你要，我可以再做一次“bulk 旧代码删除和清扫”。
