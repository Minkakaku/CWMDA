#!/usr/bin/env bash
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODE="${1:-all}"

case "$MODE" in
  all)
    Rscript "$PIPELINE_DIR/subcluster/00_subcluster_clustering.R"
    Rscript "$PIPELINE_DIR/subcluster/01_subcluster_findmarkers.R"
    Rscript "$PIPELINE_DIR/subcluster/02_subcluster_annotation.R"
    Rscript "$PIPELINE_DIR/subcluster/03_subcluster_proportion.R"
    ;;
  cluster)
    Rscript "$PIPELINE_DIR/subcluster/00_subcluster_clustering.R"
    ;;
  markers)
    Rscript "$PIPELINE_DIR/subcluster/01_subcluster_findmarkers.R"
    ;;
  annotation)
    Rscript "$PIPELINE_DIR/subcluster/02_subcluster_annotation.R"
    ;;
  proportion)
    Rscript "$PIPELINE_DIR/subcluster/03_subcluster_proportion.R"
    ;;
  *)
    echo "Usage: $0 {all|cluster|markers|annotation|proportion}" >&2
    exit 2
    ;;
esac
