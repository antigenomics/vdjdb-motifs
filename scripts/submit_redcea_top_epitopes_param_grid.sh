#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

REDCEA_VDJDB="${REDCEA_VDJDB:-vdjdb_release/vdjdb.slim.txt}"
REDCEA_SPECIES="${REDCEA_SPECIES:-HomoSapiens}"
REDCEA_EPITOPE_ROWS_FILE="${REDCEA_EPITOPE_ROWS_FILE:-}"
REDCEA_INCLUDE_NONCANONICAL="${REDCEA_INCLUDE_NONCANONICAL:-0}"
REDCEA_SBATCH_EXTRA_ARGS="${REDCEA_SBATCH_EXTRA_ARGS:-}"
REDCEA_VDBSCAN_LEIDEN_K_GRID="${REDCEA_VDBSCAN_LEIDEN_K_GRID:-4 8 12}"
REDCEA_VDBSCAN_LEIDEN_EPS_K_GRID="${REDCEA_VDBSCAN_LEIDEN_EPS_K_GRID:-4 8 12}"
REDCEA_VDBSCAN_LEIDEN_CORE_MIN_SAMPLES_GRID="${REDCEA_VDBSCAN_LEIDEN_CORE_MIN_SAMPLES_GRID:-3}"
REDCEA_VDBSCAN_LEIDEN_RESOLUTION_GRID="${REDCEA_VDBSCAN_LEIDEN_RESOLUTION_GRID:-0.5 1.0 1.5}"
REDCEA_VDBSCAN_LEIDEN_EPS_MODE_GRID="${REDCEA_VDBSCAN_LEIDEN_EPS_MODE_GRID:-sample all}"
REDCEA_VDBSCAN_LEIDEN_SYM_RULE_GRID="${REDCEA_VDBSCAN_LEIDEN_SYM_RULE_GRID:-asymmetric min}"

read -r -a k_grid <<< "$REDCEA_VDBSCAN_LEIDEN_K_GRID"
read -r -a eps_k_grid <<< "$REDCEA_VDBSCAN_LEIDEN_EPS_K_GRID"
read -r -a min_samples_grid <<< "$REDCEA_VDBSCAN_LEIDEN_CORE_MIN_SAMPLES_GRID"
read -r -a resolution_grid <<< "$REDCEA_VDBSCAN_LEIDEN_RESOLUTION_GRID"
read -r -a eps_mode_grid <<< "$REDCEA_VDBSCAN_LEIDEN_EPS_MODE_GRID"
read -r -a sym_rule_grid <<< "$REDCEA_VDBSCAN_LEIDEN_SYM_RULE_GRID"

load_epitope_rows() {
  if [[ -n "$REDCEA_EPITOPE_ROWS_FILE" ]]; then
    if [[ ! -f "$REDCEA_EPITOPE_ROWS_FILE" ]]; then
      echo "Epitope rows file not found: $REDCEA_EPITOPE_ROWS_FILE" >&2
      exit 1
    fi
    mapfile -t epitope_rows < "$REDCEA_EPITOPE_ROWS_FILE"
    return
  fi

  mapfile -t epitope_rows <<'EOF'
TRB	GILGFVFTL	5189
TRB	NLVPMVATV	4613
TRB	AVFDRKSDAK	1647
TRB	ELAGIGILTV	1426
TRB	RAKFKQLL	1378
EOF
}

load_epitope_rows

if [[ "${#epitope_rows[@]}" -eq 0 ]]; then
  echo "No epitopes selected; nothing to submit." >&2
  exit 1
fi

n_epitopes="${#epitope_rows[@]}"
n_k="${#k_grid[@]}"
n_eps_k="${#eps_k_grid[@]}"
n_min_samples="${#min_samples_grid[@]}"
n_resolution="${#resolution_grid[@]}"
n_eps_mode="${#eps_mode_grid[@]}"
n_sym_rule="${#sym_rule_grid[@]}"
total_tasks=$((n_epitopes * n_k * n_eps_k * n_min_samples * n_resolution * n_eps_mode * n_sym_rule))

echo "Selected ${n_epitopes} chain/epitope targets:"
printf '%s\n' "${epitope_rows[@]}"
echo
echo "Grid sizes:"
echo "  k=${n_k}"
echo "  eps_k=${n_eps_k}"
echo "  min_samples=${n_min_samples}"
echo "  resolution=${n_resolution}"
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
