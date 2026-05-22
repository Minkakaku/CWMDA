# Mfuzz clustering

## Overview
This directory hosts an Mfuzz pipeline script for time-course clustering.

## Inputs
- **Expression matrix**: genes x samples with normalized expression values.
- **Time or condition metadata**: vector or table describing sample ordering.

## Outputs
- **Cluster assignments**: gene-to-cluster membership.
- **Cluster profiles**: centroid expression patterns per cluster.

## Dependencies
- R (>= 3.6 recommended)
- Mfuzz
- Biobase

## Example command
```bash
Rscript Mfuzz_pip.R \
  --expression data/transcriptome/timecourse/expression_matrix.tsv \
  --outdir results/clustering/mfuzz
```

> The current Mfuzz script is a template; customize parameters (cluster number, fuzzifier) to match your experiment.
