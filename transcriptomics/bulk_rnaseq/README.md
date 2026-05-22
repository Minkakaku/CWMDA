# Bulk RNA-seq workflow

本目录包含 Bulk RNA-seq 的上游处理与下游分析脚本，并对流程进行统一整理。

## 上游处理（SRA → FASTQ → QC → 比对 → 计数 → 差异分析）

入口脚本：`omics/transcriptomics/scripts/run_bulk_rnaseq_pipeline.py`。

### 需要的输入

- `omics/transcriptomics/config/transcriptome_config.yaml` 中 `bulk_rnaseq` 配置段
- 样本元信息表：`bulk_rnaseq.sample_sheet`
- SRA accession 列表：`bulk_rnaseq.accession_list`（可由样本表生成）
- STAR 索引目录：`bulk_rnaseq.alignment.star_index`
- GTF 注释：`bulk_rnaseq.counting.reference_gtf`
- 差异分析样本信息表：`bulk_rnaseq.differential.sample_sheet`

### 输出目录结构

所有输出写入 `bulk_rnaseq.output_dir` 下的分层目录：

- `00.rawdata`：SRA 下载数据
- `01.fastq`：转换后的 FASTQ
- `02.qc`：清洗后的 reads 及 QC 报告
- `03.align`：STAR 比对产物与 BAM 索引
- `04.counts`：featureCounts 结果
- `05.diff`：差异分析输出

### 运行方式

```bash
python -m omics.transcriptomics.scripts.run_bulk_rnaseq_pipeline
```

更完整的参数说明见：`omics/transcriptomics/docs/bulk_rnaseq_pipeline.md`。

## 一站式下游分析入口

使用统一脚本串联 QC、差异分析以及下游分析（可选）：

```bash
Rscript transcriptomics/bulk_rnaseq/run_bulk_rnaseq_analysis.R \
  --counts counts_matrix.tsv \
  --sample-sheet sample_sheet.tsv \
  --contrast case,control \
  --method deseq2 \
  --run-enrichment \
  --run-mfuzz \
  --run-wgcna \
  --outdir results/bulk_rnaseq
```

如果已有 QC 输出，可使用 `--skip-qc` 并提供 `--counts`（差异分析）或 `--expression`（下游分析）作为输入。

## 下游分析（表达矩阵质控与扩展分析）

1. `01_counts/main_function.r`：表达矩阵 QC 与差异分析前准备。
2. `02_differential/differential_expression.r`：基于 DESeq2、edgeR 或 limma 的差异分析。
3. `pseudotime/` 与 `wgcna/`：拟时序（Mfuzz 等）与 WGCNA 等扩展分析。

### Main QC and matrix preparation

```bash
Rscript transcriptomics/bulk_rnaseq/01_counts/main_function.r \
  --counts counts_matrix.tsv \
  --sample-sheet sample_sheet.tsv \
  --sample-id-col sample_id \
  --group-col group \
  --min-counts 1 \
  --min-samples 2 \
  --outdir results/expression_qc
```

### Differential expression

```bash
Rscript transcriptomics/bulk_rnaseq/02_differential/differential_expression.r \
  --counts counts_matrix.tsv \
  --sample-sheet sample_sheet.tsv \
  --contrast case,control \
  --method deseq2 \
  --outdir results/differential_expression
```

### Mfuzz clustering

```bash
Rscript transcriptomics/bulk_rnaseq/pseudotime/mfuzz/mfuzz.r \
  --expression results/expression_qc/log2cpm.filtered.tsv \
  --sample-sheet sample_sheet.tsv \
  --time-col time \
  --clusters 6 \
  --outdir results/mfuzz
```

### WGCNA

```bash
Rscript transcriptomics/bulk_rnaseq/wgcna/wgcna.r \
  --expression results/expression_qc/log2cpm.filtered.tsv \
  --outdir results/wgcna
```

### GO/KEGG enrichment

```bash
Rscript transcriptomics/bulk_rnaseq/enrichment/GO&KEGG.R \
  --diff results/differential_expression/deseq2/all_results.tsv \
  --logfc-threshold 1 \
  --padj-threshold 0.05 \
  --keytype SYMBOL \
  --orgdb org.Mm.eg.db \
  --organism mmu \
  --outdir results/enrichment
```
