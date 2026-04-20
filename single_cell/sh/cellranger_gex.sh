#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <fastq_dir> <reference_dir> <out_dir> [cellranger_bin] [localcores] [localmem]"
  exit 1
fi

FASTQ_DIR="$1"
REFERENCE_DIR="$2"
OUT_DIR="$3"
CELLRANGER_BIN="${4:-cellranger}"
LOCAL_CORES="${5:-32}"
LOCAL_MEM="${6:-64}"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

mapfile -t SAMPLE_IDS < <(
  find "${FASTQ_DIR}" -maxdepth 1 -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) -printf "%f\n" |
    cut -d "_" -f 1 |
    sort -u
)

if [[ ${#SAMPLE_IDS[@]} -eq 0 ]]; then
  echo "No FASTQ files found under ${FASTQ_DIR}"
  exit 1
fi

echo "======= start cellranger at $(date) ======="

for SAMPLE_ID in "${SAMPLE_IDS[@]}"; do
  "${CELLRANGER_BIN}" count \
    --id "${SAMPLE_ID}" \
    --fastqs "${FASTQ_DIR}" \
    --transcriptome "${REFERENCE_DIR}" \
    --include-introns true \
    --sample "${SAMPLE_ID}" \
    --localcores "${LOCAL_CORES}" \
    --localmem "${LOCAL_MEM}"
done

echo "======= cellranger finished at $(date) ======="
