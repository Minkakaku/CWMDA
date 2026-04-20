#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${1:-sc_pipeline_r45}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
ENV_FILE="${2:-${ROOT_DIR}/environment.yml}"

if ! command -v conda >/dev/null 2>&1; then
  echo "Error: conda is not available in PATH."
  exit 1
fi

if [[ ! -f "${ENV_FILE}" ]]; then
  echo "Error: environment file not found: ${ENV_FILE}"
  exit 1
fi

TMP_YML="$(mktemp "${TMPDIR:-/tmp}/environment.XXXXXX.yml")"
sed "s/^name: .*/name: ${ENV_NAME}/" "${ENV_FILE}" > "${TMP_YML}"

TMP_CONDARC="$(mktemp "${TMPDIR:-/tmp}/condarc.XXXXXX.yml")"
cat > "${TMP_CONDARC}" <<'EOF'
show_channel_urls: true
channel_priority: flexible
channels:
  - conda-forge
  - bioconda
  - defaults
default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/msys2
custom_channels:
  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  pytorch: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  nvidia: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
EOF

export CONDARC="${TMP_CONDARC}"

if conda env list | awk 'NR>2 {print $1}' | grep -Fxq "${ENV_NAME}"; then
  conda env update -n "${ENV_NAME}" -f "${TMP_YML}" --prune
else
  conda env create -f "${TMP_YML}"
fi

rm -f "${TMP_YML}"
rm -f "${TMP_CONDARC}"

R_REPO='options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))'

conda run -n "${ENV_NAME}" Rscript -e "${R_REPO}; if (!requireNamespace('remotes', quietly = TRUE)) stop('r-remotes is missing'); remotes::install_github('satijalab/seurat-wrappers', upgrade = 'never', dependencies = FALSE)"
conda run -n "${ENV_NAME}" Rscript -e "${R_REPO}; if (!requireNamespace('remotes', quietly = TRUE)) stop('r-remotes is missing'); remotes::install_github('Japrin/STARTRAC', upgrade = 'never', dependencies = FALSE)"

conda run -n "${ENV_NAME}" Rscript -e "${R_REPO}; pkgs <- c('Seurat', 'SeuratWrappers', 'SingleCellExperiment', 'scDblFinder', 'celda', 'dplyr', 'ggplot2', 'patchwork', 'scales', 'harmony', 'Startrac', 'pheatmap'); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; if (length(missing) > 0) stop('Missing packages: ', paste(missing, collapse = ', ')); cat('All required packages are available.\n')"

echo "Conda environment ready: ${ENV_NAME}"
echo "Activate with: conda activate ${ENV_NAME}"
