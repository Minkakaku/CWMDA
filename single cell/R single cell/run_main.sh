#!/usr/bin/env bash
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
case "$(uname -s 2>/dev/null || echo unknown)" in
  MINGW*|MSYS*|CYGWIN*)
    PROJECT_DIR="${PROJECT_DIR:-$(pwd -P)}"
    ;;
  *)
    PROJECT_DIR="${PROJECT_DIR:-$(cd "$PIPELINE_DIR/.." && pwd)}"
    ;;
esac
SAMPLE_SHEET="${SAMPLE_SHEET:-$PROJECT_DIR/01_QC/sample_groups.csv}"
QC_FILTER_SHEET="${QC_FILTER_SHEET:-$PROJECT_DIR/01_QC/sample_qc_thresholds.csv}"
MODE="${1:-all}"

run_preflight() {
  Rscript "$PIPELINE_DIR/main/00_preflight_check.R"
}

run_qc1() {
  Rscript "$PIPELINE_DIR/main/01_qc1.R"
}

run_after_qc() {
  if [[ ! -f "$QC_FILTER_SHEET" ]]; then
    echo "Sample QC threshold CSV not found: $QC_FILTER_SHEET" >&2
    echo "Run qc1 first, edit threshold columns, then rerun: $0 after_qc" >&2
    exit 2
  fi

  local had_sample_sheet=0
  if [[ -f "$SAMPLE_SHEET" ]]; then
    had_sample_sheet=1
  fi

  Rscript "$PIPELINE_DIR/main/01_qc2.R"
  if [[ "$had_sample_sheet" -eq 0 ]]; then
    echo "Sample group CSV created: $SAMPLE_SHEET" >&2
    echo "Edit the group column, then rerun: $0 after_qc" >&2
    exit 0
  fi

  Rscript "$PIPELINE_DIR/main/02_globalcluster.R"
  Rscript "$PIPELINE_DIR/main/03_findmarkers.R"
  Rscript "$PIPELINE_DIR/main/04_cell_annotation.R"
  Rscript "$PIPELINE_DIR/main/05_cell_proportion.R"
}

case "$MODE" in
  all)
    had_qc_filter_sheet=0
    if [[ -f "$QC_FILTER_SHEET" ]]; then
      had_qc_filter_sheet=1
    fi

    run_preflight
    run_qc1
    if [[ "$had_qc_filter_sheet" -eq 0 ]]; then
      echo "Sample QC threshold CSV created: $QC_FILTER_SHEET" >&2
      echo "Edit threshold columns, then rerun: $0 after_qc" >&2
      exit 0
    fi

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
    if [[ ! -f "$QC_FILTER_SHEET" ]]; then
      echo "Sample QC threshold CSV not found: $QC_FILTER_SHEET" >&2
      echo "Run qc1 first, edit threshold columns, then rerun: $0 qc2" >&2
      exit 2
    fi
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
