#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TCRNET_NPROC="${TCRNET_NPROC:-8}"

cd "$ROOT_DIR"
export TCRNET_NPROC
Rscript -e "rmarkdown::render('tcrnet/compute_vdjdb_motifs.Rmd')"
python scripts/generate_tcrnet_reports.py
