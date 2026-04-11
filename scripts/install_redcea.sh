#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

CONDA_ENV_NAME="vdjdb-redcea"
PYTHON_VERSION="3.11"

if ! conda env list | awk '{print $1}' | grep -Fxq "$CONDA_ENV_NAME"; then
  conda create -n "$CONDA_ENV_NAME" "python=$PYTHON_VERSION" -y
fi

conda run -n "$CONDA_ENV_NAME" python -m pip install --upgrade pip
conda run -n "$CONDA_ENV_NAME" python -m pip install -e redcea

echo "REDCEA conda environment is ready: $CONDA_ENV_NAME"
echo "Activate it with: conda activate $CONDA_ENV_NAME"
