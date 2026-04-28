#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${1:-miptd}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

if ! command -v conda >/dev/null 2>&1; then
  echo "conda was not found in PATH." >&2
  exit 1
fi

echo "[MIPTD install] Step 1: create minimal conda environment: ${ENV_NAME}"
conda create -n "${ENV_NAME}" -y python=3.11 pip

echo "[MIPTD install] Step 2: install mamba into ${ENV_NAME}"
conda install -n "${ENV_NAME}" -y -c conda-forge mamba

echo "[MIPTD install] Step 3: install conda dependencies with mamba"
conda run -n "${ENV_NAME}" conda install -y -c pytorch -c bioconda -c conda-forge \
  nodejs \
  pytorch \
  rdkit \
  r-base \
  r-tidyverse \
  r-readxl \
  r-writexl \
  r-ggvenn \
  r-jsonlite \
  bioconductor-clusterprofiler \
  bioconductor-org.hs.eg.db \
  r-ggnewscale \
  r-circlize \
  bioconductor-complexheatmap \
  matplotlib-base \
  pillow

echo "[MIPTD install] Step 4: install R package legendry"
conda run -n "${ENV_NAME}" Rscript -e "install.packages('legendry', repos='https://cloud.r-project.org')"

echo "[MIPTD install] Step 5: install Python package in editable mode"
conda run -n "${ENV_NAME}" pip install -e "${PROJECT_ROOT}" chemprop==2.2.2

echo "[MIPTD install] Completed."
echo "Activate with: conda activate ${ENV_NAME}"
echo "Verify with: conda run -n ${ENV_NAME} MIPTD --help"
