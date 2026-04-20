#!/bin/bash
#SBATCH --job-name=redcea-grid
#SBATCH --array=0-35
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --partition=short
#SBATCH --output=logs/redcea_grid_%A_%a.out
#SBATCH --error=logs/redcea_grid_%A_%a.err

set -euo pipefail

ROOT_DIR="/home/evlasova/vdjdb-motifs"
REDCEA_CHAINS="${REDCEA_CHAINS:-TRA TRB}"
REDCEA_K_GRID="${REDCEA_K_GRID:-3 5 10 15 30 50}"
REDCEA_RESOLUTION_GRID="${REDCEA_RESOLUTION_GRID:-0.5 1 2}"
REDCEA_EPITOPES="${REDCEA_EPITOPES:-GLCTLVAML YLQPRTFLL}"
REDCEA_OUTPUT="${REDCEA_OUTPUT:-results/redcea_param_grid}"
REDCEA_TCREMP_CACHE_DIR="${REDCEA_TCREMP_CACHE_DIR:-/home/evlasova/vdjdb-motifs/results/redcea/tcremp}"
REDCEA_SKIP_UMAP="${REDCEA_SKIP_UMAP:-1}"

read -r -a chains <<< "$REDCEA_CHAINS"
read -r -a k_grid <<< "$REDCEA_K_GRID"
read -r -a resolution_grid <<< "$REDCEA_RESOLUTION_GRID"

n_chains="${#chains[@]}"
n_k="${#k_grid[@]}"
n_resolution="${#resolution_grid[@]}"
total_tasks=$((n_chains * n_k * n_resolution))
task_id="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

if (( task_id < 0 || task_id >= total_tasks )); then
  echo "Task ${task_id} is outside grid size ${total_tasks}" >&2
  exit 1
fi

resolution_index=$((task_id % n_resolution))
k_index=$(((task_id / n_resolution) % n_k))
chain_index=$((task_id / (n_resolution * n_k)))

chain="${chains[$chain_index]}"
k_neighbors="${k_grid[$k_index]}"
leiden_resolution="${resolution_grid[$resolution_index]}"
resolution_tag="${leiden_resolution//./p}"
output_tag="k${k_neighbors}_res${resolution_tag}"

cd "$ROOT_DIR"
mkdir -p logs "$REDCEA_OUTPUT"

export REDCEA_CHAIN="$chain"
export REDCEA_K_NEIGHBORS="$k_neighbors"
export REDCEA_LEIDEN_RESOLUTION="$leiden_resolution"
export REDCEA_OUTPUT_TAG="$output_tag"
export REDCEA_OUTPUT
export REDCEA_TCREMP_CACHE_DIR
export REDCEA_EPITOPES
export REDCEA_SKIP_UMAP
export REDCEA_NPROC="${SLURM_CPUS_PER_TASK}"

echo "Running task ${task_id}/${total_tasks}: chain=${REDCEA_CHAIN}, k=${REDCEA_K_NEIGHBORS}, resolution=${REDCEA_LEIDEN_RESOLUTION}"
echo "Epitopes: ${REDCEA_EPITOPES}"
echo "TCRemP cache: ${REDCEA_TCREMP_CACHE_DIR}"
echo "Output: ${REDCEA_OUTPUT}/cluster_members_${REDCEA_CHAIN}_${REDCEA_OUTPUT_TAG}.txt"

bash ./scripts/run_redcea_ylq.sh
