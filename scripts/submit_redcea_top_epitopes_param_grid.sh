#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

REDCEA_VDJDB="${REDCEA_VDJDB:-vdjdb_release/vdjdb.slim.txt}"
REDCEA_SPECIES="${REDCEA_SPECIES:-HomoSapiens}"
REDCEA_TOP_EPITOPES_N="${REDCEA_TOP_EPITOPES_N:-10}"
REDCEA_TOP_EPITOPES_MIN_COUNT="${REDCEA_TOP_EPITOPES_MIN_COUNT:-1}"
REDCEA_CHAINS="${REDCEA_CHAINS:-TRB}"
REDCEA_INCLUDE_NONCANONICAL="${REDCEA_INCLUDE_NONCANONICAL:-0}"
REDCEA_SBATCH_EXTRA_ARGS="${REDCEA_SBATCH_EXTRA_ARGS:-}"
REDCEA_VDBSCAN_LEIDEN_K_GRID="${REDCEA_VDBSCAN_LEIDEN_K_GRID:-5 8 12}"
REDCEA_VDBSCAN_LEIDEN_EPS_K_GRID="${REDCEA_VDBSCAN_LEIDEN_EPS_K_GRID:-5 8 12}"
REDCEA_VDBSCAN_LEIDEN_CORE_MIN_SAMPLES_GRID="${REDCEA_VDBSCAN_LEIDEN_CORE_MIN_SAMPLES_GRID:-2 3 5}"
REDCEA_VDBSCAN_LEIDEN_RESOLUTION_GRID="${REDCEA_VDBSCAN_LEIDEN_RESOLUTION_GRID:-1.0 1.5 2.0}"
REDCEA_VDBSCAN_LEIDEN_SUB_RESOLUTION_GRID="${REDCEA_VDBSCAN_LEIDEN_SUB_RESOLUTION_GRID:-1.0 1.5}"
REDCEA_VDBSCAN_LEIDEN_EPS_MODE_GRID="${REDCEA_VDBSCAN_LEIDEN_EPS_MODE_GRID:-sample background all}"
REDCEA_VDBSCAN_LEIDEN_SYM_RULE_GRID="${REDCEA_VDBSCAN_LEIDEN_SYM_RULE_GRID:-asymmetric min max}"

read -r -a chains <<< "$REDCEA_CHAINS"
read -r -a k_grid <<< "$REDCEA_VDBSCAN_LEIDEN_K_GRID"
read -r -a eps_k_grid <<< "$REDCEA_VDBSCAN_LEIDEN_EPS_K_GRID"
read -r -a min_samples_grid <<< "$REDCEA_VDBSCAN_LEIDEN_CORE_MIN_SAMPLES_GRID"
read -r -a resolution_grid <<< "$REDCEA_VDBSCAN_LEIDEN_RESOLUTION_GRID"
read -r -a sub_resolution_grid <<< "$REDCEA_VDBSCAN_LEIDEN_SUB_RESOLUTION_GRID"
read -r -a eps_mode_grid <<< "$REDCEA_VDBSCAN_LEIDEN_EPS_MODE_GRID"
read -r -a sym_rule_grid <<< "$REDCEA_VDBSCAN_LEIDEN_SYM_RULE_GRID"

helper_args=(
  --vdjdb "$REDCEA_VDJDB"
  --chains "${chains[@]}"
  --species "$REDCEA_SPECIES"
  --top-n "$REDCEA_TOP_EPITOPES_N"
  --min-count "$REDCEA_TOP_EPITOPES_MIN_COUNT"
  --format tsv
)
if [[ "$REDCEA_INCLUDE_NONCANONICAL" == "1" ]]; then
  helper_args+=(--include-noncanonical)
fi

mapfile -t epitope_rows < <(
  python scripts/select_top_vdjdb_epitopes.py "${helper_args[@]}" | tail -n +2
)

if [[ "${#epitope_rows[@]}" -eq 0 ]]; then
  echo "No epitopes selected; nothing to submit." >&2
  exit 1
fi

n_epitopes="${#epitope_rows[@]}"
n_k="${#k_grid[@]}"
n_eps_k="${#eps_k_grid[@]}"
n_min_samples="${#min_samples_grid[@]}"
n_resolution="${#resolution_grid[@]}"
n_sub_resolution="${#sub_resolution_grid[@]}"
n_eps_mode="${#eps_mode_grid[@]}"
n_sym_rule="${#sym_rule_grid[@]}"
total_tasks=$((n_epitopes * n_k * n_eps_k * n_min_samples * n_resolution * n_sub_resolution * n_eps_mode * n_sym_rule))

echo "Selected ${n_epitopes} chain/epitope targets:"
printf '%s\n' "${epitope_rows[@]}"
echo
echo "Grid sizes:"
echo "  k=${n_k}"
echo "  eps_k=${n_eps_k}"
echo "  min_samples=${n_min_samples}"
echo "  resolution=${n_resolution}"
echo "  sub_resolution=${n_sub_resolution}"
echo "  eps_mode=${n_eps_mode}"
echo "  sym_rule=${n_sym_rule}"
echo "Total tasks: ${total_tasks}"
echo

sbatch_args=(--array "0-$((total_tasks - 1))")
if [[ -n "$REDCEA_SBATCH_EXTRA_ARGS" ]]; then
  # shellcheck disable=SC2206
  extra_args=( $REDCEA_SBATCH_EXTRA_ARGS )
  sbatch_args+=("${extra_args[@]}")
fi

echo "Submitting: sbatch ${sbatch_args[*]} scripts/run_redcea_top_epitopes_param_grid.slurm"
sbatch "${sbatch_args[@]}" scripts/run_redcea_top_epitopes_param_grid.slurm
