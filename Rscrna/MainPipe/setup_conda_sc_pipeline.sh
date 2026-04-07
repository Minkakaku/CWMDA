#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${1:-sc_pipeline_r45}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v conda >/dev/null 2>&1; then
  echo "Error: conda 不在 PATH 里，请先加载 conda。"
  exit 1
fi

# 临时 environment.yml
TMP_YML="$(mktemp "${TMPDIR:-/tmp}/environment.XXXXXX.yml")"
sed "s/^name: .*/name: ${ENV_NAME}/" "${SCRIPT_DIR}/environment.yml" > "${TMP_YML}"

# 临时 condarc：优先使用清华源
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
  conda env update --environment-spec environment.yml -n "${ENV_NAME}" -f "${TMP_YML}" --prune
else
  conda env create --environment-spec environment.yml -f "${TMP_YML}"
fi

rm -f "${TMP_YML}"
rm -f "${TMP_CONDARC}"

# R 使用清华 CRAN 镜像
R_REPO='options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))'

conda run -n "${ENV_NAME}" Rscript -e "${R_REPO}; if (!requireNamespace('remotes', quietly = TRUE)) stop('r-remotes 未安装'); remotes::install_github('satijalab/seurat-wrappers', upgrade = 'never', dependencies = FALSE)"

conda run -n "${ENV_NAME}" Rscript -e "${R_REPO}; if (!requireNamespace('remotes', quietly = TRUE)) stop('r-remotes 未安装'); remotes::install_github('Japrin/STARTRAC', upgrade = 'never', dependencies = FALSE)"

conda run -n "${ENV_NAME}" Rscript -e "${R_REPO}; pkgs <- c('Seurat', 'SeuratWrappers', 'SingleCellExperiment', 'scDblFinder', 'celda', 'dplyr', 'ggplot2', 'patchwork', 'scales', 'harmony', 'Startrac', 'pheatmap'); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; if (length(missing) > 0) stop('以下包未安装成功: ', paste(missing, collapse = ', ')); cat('All required packages are available.\n')"

echo "Conda 环境已准备好：${ENV_NAME}"
echo "激活命令：conda activate ${ENV_NAME}"