#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 7 ]]; then
  echo "Usage: $0 <expression_mtx> <grn_out> <ctx_out> <auc_out> <tf_list> <motif_annotations> <ranking_db> [pyscenic_bin] [workers]"
  exit 1
fi

EXPR_MTX="$1"
GRN_OUT="$2"
CTX_OUT="$3"
AUC_OUT="$4"
TF_LIST="$5"
MOTIF_ANNOTATIONS="$6"
RANKING_DB="$7"
PYSCENIC_BIN="${8:-pyscenic}"
WORKERS="${9:-32}"

"${PYSCENIC_BIN}" grn \
  --num_workers "${WORKERS}" \
  --sparse \
  --method grnboost2 \
  --output "${GRN_OUT}" \
  "${EXPR_MTX}" \
  "${TF_LIST}"

"${PYSCENIC_BIN}" ctx \
  --num_workers "${WORKERS}" \
  --output "${CTX_OUT}" \
  --expression_mtx_fname "${EXPR_MTX}" \
  --mode custom_multiprocessing \
  --annotations_fname "${MOTIF_ANNOTATIONS}" \
  "${GRN_OUT}" \
  "${RANKING_DB}"

"${PYSCENIC_BIN}" aucell \
  --num_workers "${WORKERS}" \
  --output "${AUC_OUT}" \
  "${EXPR_MTX}" \
  "${CTX_OUT}"
