#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <output_root_dir>"
  exit 1
fi

OUTPUT_ROOT_DIR="$1"

for DIR in "${OUTPUT_ROOT_DIR}"/*/; do
  [[ -d "${DIR}" ]] || continue
  rm -rf "${DIR}/velocyto"
done
