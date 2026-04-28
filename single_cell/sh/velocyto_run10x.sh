#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "Usage: $0 <sample_root_dir> <reference_dir> <out_root_dir> <cellranger_script>"
  exit 1
fi

SAMPLE_ROOT_DIR="$1"
REFERENCE_DIR="$2"
OUT_ROOT_DIR="$3"
CELLRANGER_SCRIPT="$4"

mkdir -p "${OUT_ROOT_DIR}"

for DIR in "${SAMPLE_ROOT_DIR}"/*/; do
  [[ -d "${DIR}" ]] || continue
  "${CELLRANGER_SCRIPT}" "${DIR}" "${REFERENCE_DIR}" "${OUT_ROOT_DIR}"
done
