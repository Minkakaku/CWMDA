# 基因组分析流程说明

## 流程概览

1. 输入检查（BAM + 参考基因组）
2. 质量控制（输出统一的 TSV + HTML）
3. 变异检测（占位说明，后续接入真实工具）

## 核心脚本

- `omics/genomics/scripts/run_genome_pipeline.py`：基因组流程入口示例。

## 输入输出

- 输入：
  - `omics/genomics/config/genome_config.yaml` 中配置 `input_bam`、`reference_fasta`。
- 输出：
  - `data/genome/output/qc_stats.tsv`
  - `data/genome/output/qc_report.html`

## 运行方式

```bash
python -m omics.genomics.scripts.run_genome_pipeline
```

## 常见问题

- 若报错 `Missing 输入 BAM`：请检查 `omics/genomics/config/genome_config.yaml` 中的 `input_bam` 路径。
- 若报错 `Missing 参考基因组`：请检查 `reference_fasta`。
