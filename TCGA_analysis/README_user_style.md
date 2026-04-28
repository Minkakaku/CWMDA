# TCGA PANoRG Analysis: User-Style Pipeline

这份说明对应当前实际运行的新版脚本。

## 处理思路

完全按你的 `04tcga.R` 方式处理 TCGA：

- `GDCquery` + `GDCdownload` + `GDCprepare`
- `STAR - Counts`
- `rowData(se)$gene_id` 去版本号
- duplicated `Ensembl` 用 `rowsum()` 合并
- `colData(se)` 直接取临床和 subtype 信息
- `shortLetterCode == "TP"` 保留原发肿瘤
- `edgeR + voom` 生成表达矩阵

## 当前癌种

- `TCGA_COAD`
- `TCGA_CRC` = `TCGA-COAD + TCGA-READ`
- `TCGA_PRAD`

## 实际运行脚本

- `scripts/00_install_packages.R`
- `scripts/01_download_tcga_data.R`
- `scripts/02_run_panorg_analysis.R`

## 实际运行模块

- `R/tcga_user_style_utils.R`
- `R/tcga_user_style_prepare.R`
- `R/tcga_user_style_panorg.R`

## 运行

```r
Rscript scripts/00_install_packages.R
Rscript scripts/01_download_tcga_data.R
Rscript scripts/02_run_panorg_analysis.R
```

## 说明

旧版 `R/tcga_download.R`、`R/panorg_analysis.R` 和旧 README 没有被新版脚本调用。
