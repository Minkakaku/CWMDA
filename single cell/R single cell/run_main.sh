#!/usr/bin/env bash
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODE="${1:-all}"

run_preflight() {
  Rscript "$PIPELINE_DIR/main/00_preflight_check.R"
}

run_qc1() {
  Rscript "$PIPELINE_DIR/main/01_qc1.R"
}

run_after_qc() {
  Rscript "$PIPELINE_DIR/main/01_qc2.R"
  Rscript "$PIPELINE_DIR/main/02_globalcluster.R"
  Rscript "$PIPELINE_DIR/main/03_findmarkers.R"
  Rscript "$PIPELINE_DIR/main/04_cell_annotation.R"
  Rscript "$PIPELINE_DIR/main/05_cell_proportion.R"
}

case "$MODE" in
  all)
    run_preflight
    run_qc1
    run_after_qc
    ;;
  preflight)
    run_preflight
    ;;
  qc1)
    run_qc1
    ;;
  after_qc)
    run_after_qc
    ;;
  qc2)
    Rscript "$PIPELINE_DIR/main/01_qc2.R"
    ;;
  cluster)
    Rscript "$PIPELINE_DIR/main/02_globalcluster.R"
    ;;
  markers)
    Rscript "$PIPELINE_DIR/main/03_findmarkers.R"
    ;;
  annotation)
    Rscript "$PIPELINE_DIR/main/04_cell_annotation.R"
    ;;
  proportion)
    Rscript "$PIPELINE_DIR/main/05_cell_proportion.R"
    ;;
  *)
    echo "Usage: $0 {all|preflight|qc1|after_qc|qc2|cluster|markers|annotation|proportion}" >&2
    exit 2
    ;;
esac
