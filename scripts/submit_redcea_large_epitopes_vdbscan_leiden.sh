#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

REDCEA_VDJDB="${REDCEA_VDJDB:-vdjdb_release/vdjdb.slim.txt}"
REDCEA_SPECIES="${REDCEA_SPECIES:-HomoSapiens}"
REDCEA_CHAINS="${REDCEA_CHAINS:-TRB}"
REDCEA_TOP_N="${REDCEA_TOP_N:-20}"
REDCEA_MIN_EPITOPE_COUNT="${REDCEA_MIN_EPITOPE_COUNT:-1000}"
REDCEA_EPITOPE_ROWS_FILE="${REDCEA_EPITOPE_ROWS_FILE:-}"
REDCEA_GENERATED_ROWS_DIR="${REDCEA_GENERATED_ROWS_DIR:-.tmp/redcea_large_epitopes}"
REDCEA_OUTPUT_ROOT="${REDCEA_OUTPUT_ROOT:-results/redcea_vdbscan_leiden_large_epitopes}"
REDCEA_TCREMP_CACHE_ROOT="${REDCEA_TCREMP_CACHE_ROOT:-results/redcea_vdbscan_leiden_large_epitopes_cache}"
REDCEA_PREP_ROOT="${REDCEA_PREP_ROOT:-results/redcea_vdbscan_leiden_large_epitopes_prep}"
REDCEA_INCLUDE_NONCANONICAL="${REDCEA_INCLUDE_NONCANONICAL:-0}"
REDCEA_SBATCH_EXTRA_ARGS="${REDCEA_SBATCH_EXTRA_ARGS:-}"

REDCEA_VDBSCAN_LEIDEN_K_GRID="${REDCEA_VDBSCAN_LEIDEN_K_GRID:-8 12}"
REDCEA_VDBSCAN_LEIDEN_EPS_K_GRID="${REDCEA_VDBSCAN_LEIDEN_EPS_K_GRID:-4 8}"
REDCEA_VDBSCAN_LEIDEN_CORE_MIN_SAMPLES_GRID="${REDCEA_VDBSCAN_LEIDEN_CORE_MIN_SAMPLES_GRID:-3}"
REDCEA_VDBSCAN_LEIDEN_RESOLUTION_GRID="${REDCEA_VDBSCAN_LEIDEN_RESOLUTION_GRID:-0.5 1.0}"
REDCEA_VDBSCAN_LEIDEN_EPS_MODE_GRID="${REDCEA_VDBSCAN_LEIDEN_EPS_MODE_GRID:-sample background}"
REDCEA_VDBSCAN_LEIDEN_SYM_RULE_GRID="${REDCEA_VDBSCAN_LEIDEN_SYM_RULE_GRID:-asymmetric}"

read -r -a chains <<< "$REDCEA_CHAINS"
read -r -a k_grid <<< "$REDCEA_VDBSCAN_LEIDEN_K_GRID"
read -r -a eps_k_grid <<< "$REDCEA_VDBSCAN_LEIDEN_EPS_K_GRID"
read -r -a min_samples_grid <<< "$REDCEA_VDBSCAN_LEIDEN_CORE_MIN_SAMPLES_GRID"
read -r -a resolution_grid <<< "$REDCEA_VDBSCAN_LEIDEN_RESOLUTION_GRID"
read -r -a eps_mode_grid <<< "$REDCEA_VDBSCAN_LEIDEN_EPS_MODE_GRID"
read -r -a sym_rule_grid <<< "$REDCEA_VDBSCAN_LEIDEN_SYM_RULE_GRID"

if [[ -z "$REDCEA_EPITOPE_ROWS_FILE" ]]; then
  mkdir -p "$REDCEA_GENERATED_ROWS_DIR"
  REDCEA_EPITOPE_ROWS_FILE="${REDCEA_GENERATED_ROWS_DIR}/large_epitopes.tsv"
  tmp_rows="${REDCEA_EPITOPE_ROWS_FILE}.tmp"
  python scripts/select_top_vdjdb_epitopes.py \
    --vdjdb "$REDCEA_VDJDB" \
    --chains "${chains[@]}" \
    --species "$REDCEA_SPECIES" \
    --top-n "$REDCEA_TOP_N" \
    --min-count "$REDCEA_MIN_EPITOPE_COUNT" \
    --format tsv > "$tmp_rows"
  tail -n +2 "$tmp_rows" > "$REDCEA_EPITOPE_ROWS_FILE"
  rm -f "$tmp_rows"
fi

if [[ ! -f "$REDCEA_EPITOPE_ROWS_FILE" ]]; then
  echo "Epitope rows file not found: $REDCEA_EPITOPE_ROWS_FILE" >&2
  exit 1
fi

mapfile -t epitope_rows < <(grep -v $'^chain\t' "$REDCEA_EPITOPE_ROWS_FILE" | sed '/^[[:space:]]*$/d')
if [[ "${#epitope_rows[@]}" -eq 0 ]]; then
  echo "No large epitopes selected; check REDCEA_MIN_EPITOPE_COUNT or REDCEA_EPITOPE_ROWS_FILE." >&2
  exit 1
fi

echo "Selected ${#epitope_rows[@]} large epitopes:"
printf '%s\n' "${epitope_rows[@]}"
echo

submission_count=0
for k in "${k_grid[@]}"; do
  for eps_k in "${eps_k_grid[@]}"; do
    for min_samples in "${min_samples_grid[@]}"; do
      for resolution in "${resolution_grid[@]}"; do
        for eps_mode in "${eps_mode_grid[@]}"; do
          for sym_rule in "${sym_rule_grid[@]}"; do
            resolution_tag="${resolution//./p}"
            job_name="rvdbl_k${k}_ek${eps_k}_r${resolution_tag}_m${min_samples}_${eps_mode}_${sym_rule}"
            sbatch_args=(--job-name "$job_name")
            if [[ -n "$REDCEA_SBATCH_EXTRA_ARGS" ]]; then
              # shellcheck disable=SC2206
              extra_args=( $REDCEA_SBATCH_EXTRA_ARGS )
              sbatch_args+=("${extra_args[@]}")
            fi

            echo "Submitting ${job_name}"
            sbatch \
              "${sbatch_args[@]}" \
              --export="ALL,ROOT_DIR=${ROOT_DIR},REDCEA_VDJDB=${REDCEA_VDJDB},REDCEA_SPECIES=${REDCEA_SPECIES},REDCEA_EPITOPE_ROWS_FILE=${REDCEA_EPITOPE_ROWS_FILE},REDCEA_OUTPUT_ROOT=${REDCEA_OUTPUT_ROOT},REDCEA_TCREMP_CACHE_ROOT=${REDCEA_TCREMP_CACHE_ROOT},REDCEA_PREP_ROOT=${REDCEA_PREP_ROOT},REDCEA_INCLUDE_NONCANONICAL=${REDCEA_INCLUDE_NONCANONICAL},REDCEA_K_NEIGHBORS=${k},REDCEA_EPS_K_NEIGHBORS=${eps_k},REDCEA_CORE_MIN_SAMPLES=${min_samples},REDCEA_LEIDEN_RESOLUTION=${resolution},REDCEA_EPS_MODE=${eps_mode},REDCEA_SYM_RULE=${sym_rule}" \
              scripts/run_redcea_large_epitopes_vdbscan_leiden.slurm
            submission_count=$((submission_count + 1))
          done
        done
      done
    done
  done
done

echo
echo "Submitted ${submission_count} jobs."
echo "Epitope rows file: ${REDCEA_EPITOPE_ROWS_FILE}"
