# TME Deconvolution (IOBR)

本目录提供基于 IOBR 的 TME 去卷积流程，覆盖 8 个方法：

- CIBERSORT
- CIBERSORT-ABS
- quanTIseq
- xCell
- MCP-counter
- EPIC
- TIMER
- ESTIMATE

## 快速开始

1. 修改同目录的 `config.yaml`：
   - `expression_file`: 你的表达矩阵（基因在行，样本在列）
   - `separator`: 表达矩阵的分隔符（默认 `\t`）
   - `gene_column`: 基因列名称或列序号（默认第 1 列）
   - `output_dir`: 结果输出目录

2. 运行脚本：

```bash
Rscript "TME_deconvolution_iobr.R"
```

## 在项目内调用（R 会话）

```r
source("transcriptomics/bulk_rnaseq/deconvolution/TME_deconvolution_iobr.R")
```

脚本会自动运行 8 种方法，并将结果输出到 `output_dir` 目录下：

- `tme_deconvolution_<method>.tsv`
- `tme_deconvolution_all_methods.rds`

## 依赖

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("IOBR")
install.packages("yaml")
```
