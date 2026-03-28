#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV_NAME="vdjdb-redcea"
PYTHON_VERSION="3.11"
TCREMP_GIT_URL="git+https://github.com/antigenomics/tcremp.git"

conda create -n "$CONDA_ENV_NAME" "python=$PYTHON_VERSION" -y

conda run -n "$CONDA_ENV_NAME" python -m pip install --upgrade pip
conda run -n "$CONDA_ENV_NAME" python -m pip install "$TCREMP_GIT_URL"
conda run -n "$CONDA_ENV_NAME" python -m pip install -e redcea

echo "REDCEA conda environment is ready: $CONDA_ENV_NAME"
echo "Activate it with: conda activate $CONDA_ENV_NAME"
