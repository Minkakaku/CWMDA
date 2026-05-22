# RNA-seq Upstream Pipeline

This folder now provides a single unified upstream pipeline script that supports:

- **Aligner**: STAR or HISAT2
- **Quantification**: featureCounts or RSEM

## Usage

```bash
./upstream_pipeline.sh --aligner star --quant featurecounts \
  --genome-dir /path/to/star/index \
  --gtf /path/to/genes.gtf
```

```bash
./upstream_pipeline.sh --aligner hisat2 --quant rsem \
  --hisat2-index /path/to/hisat2/index \
  --rsem-ref /path/to/rsem/reference
```

### Input metadata

The script uses a sample sheet by default:

```
data/transcriptome/metadata/sample_sheet.tsv
```

The file must include columns `sample_id`, `layout` (SE/PE), and `sra`.
If you do not have a sample sheet, provide `--acc-list` and (optionally)
`--default-layout` to process a simple list of SRA accessions.

### Output

The script writes to the working directory:

- `00.rawdata/` – downloaded SRA files
- `01.fastq/` – FASTQ files
- `02.cleandata/` – cleaned reads from fastp
- `03.align/` – aligned BAM files
- `04.quant/` – quantification output
- `logs/` – logs per step

## Notes

- Use `--skip-download` or `--skip-qc` to reuse existing data.
- STAR requires `--genome-dir` and HISAT2 requires `--hisat2-index`.
- featureCounts requires `--gtf` and RSEM requires `--rsem-ref`.
