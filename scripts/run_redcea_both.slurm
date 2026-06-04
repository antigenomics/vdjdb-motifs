#!/bin/bash
#SBATCH --job-name=redcea
#SBATCH --array=0-1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/redcea_%A_%a.out
#SBATCH --error=logs/redcea_%A_%a.err

set -euo pipefail

# Edit this path for your cluster checkout location if needed.
ROOT_DIR="${ROOT_DIR:-/path/to/vdjdb-motifs}"

case "${SLURM_ARRAY_TASK_ID}" in
  0) CHAIN="TRA" ;;
  1) CHAIN="TRB" ;;
  *)
    echo "Unsupported SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
    ;;
esac

cd "${ROOT_DIR}"
mkdir -p logs

export REDCEA_CHAIN="${CHAIN}"
export REDCEA_NPROC="${SLURM_CPUS_PER_TASK}"

echo "Running RedCEA for chain=${REDCEA_CHAIN} with nproc=${REDCEA_NPROC}"
bash ./scripts/run_redcea.sh
