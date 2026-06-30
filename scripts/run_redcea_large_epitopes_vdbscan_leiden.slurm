#!/bin/bash
#SBATCH --job-name=redcea-large-vdbl
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=short
#SBATCH --output=logs/redcea_large_vdbl_%x_%j.out
#SBATCH --error=logs/redcea_large_vdbl_%x_%j.err

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-/home/evlasova/vdjdb-motifs}"
REDCEA_VDJDB="${REDCEA_VDJDB:-vdjdb_release/vdjdb.slim.txt}"
REDCEA_SPECIES="${REDCEA_SPECIES:-HomoSapiens}"
REDCEA_EPITOPE_ROWS_FILE="${REDCEA_EPITOPE_ROWS_FILE:?REDCEA_EPITOPE_ROWS_FILE is required}"
REDCEA_OUTPUT_ROOT="${REDCEA_OUTPUT_ROOT:-results/redcea_vdbscan_leiden_large_epitopes}"
REDCEA_TCREMP_CACHE_ROOT="${REDCEA_TCREMP_CACHE_ROOT:-results/redcea_vdbscan_leiden_large_epitopes_cache}"
REDCEA_PREP_ROOT="${REDCEA_PREP_ROOT:-results/redcea_vdbscan_leiden_large_epitopes_prep}"
REDCEA_INCLUDE_NONCANONICAL="${REDCEA_INCLUDE_NONCANONICAL:-0}"
REDCEA_METRICS="${REDCEA_METRICS:-dissimilarity}"
REDCEA_CLUSTER_PC_COMPONENTS="${REDCEA_CLUSTER_PC_COMPONENTS:-50}"
REDCEA_N_BG_POINTS="${REDCEA_N_BG_POINTS:-100000}"
REDCEA_ENRICHMENT_TEST="${REDCEA_ENRICHMENT_TEST:-zbinom}"
REDCEA_RANDOM_SEED="${REDCEA_RANDOM_SEED:-42}"
REDCEA_DEBUG_SAVE_INTERMEDIATE="${REDCEA_DEBUG_SAVE_INTERMEDIATE:-0}"
REDCEA_ADD_AUXILIARY_CLUSTER_METRICS="${REDCEA_ADD_AUXILIARY_CLUSTER_METRICS:-0}"

REDCEA_K_NEIGHBORS="${REDCEA_K_NEIGHBORS:?REDCEA_K_NEIGHBORS is required}"
REDCEA_EPS_K_NEIGHBORS="${REDCEA_EPS_K_NEIGHBORS:?REDCEA_EPS_K_NEIGHBORS is required}"
REDCEA_CORE_MIN_SAMPLES="${REDCEA_CORE_MIN_SAMPLES:?REDCEA_CORE_MIN_SAMPLES is required}"
REDCEA_LEIDEN_RESOLUTION="${REDCEA_LEIDEN_RESOLUTION:?REDCEA_LEIDEN_RESOLUTION is required}"
REDCEA_EPS_MODE="${REDCEA_EPS_MODE:?REDCEA_EPS_MODE is required}"
REDCEA_SYM_RULE="${REDCEA_SYM_RULE:?REDCEA_SYM_RULE is required}"

TRA_BG_AIRR="${TRA_BG_AIRR:-redcea/data/backgrounds/tra_background_100k.tsv}"
TRA_BG_EMBEDDING="${TRA_BG_EMBEDDING:-redcea/data/backgrounds/tra_background_embeddings.parquet}"
REDCEA_PREPARE_VJ_BACKGROUND="${REDCEA_PREPARE_VJ_BACKGROUND:-1}"
REDCEA_VJ_BACKGROUND_CHAIN="${REDCEA_VJ_BACKGROUND_CHAIN:-TRB}"
REDCEA_VJ_BACKGROUND_SEED="${REDCEA_VJ_BACKGROUND_SEED:-42}"
REDCEA_VJ_BACKGROUND_KEEP_ALLELES="${REDCEA_VJ_BACKGROUND_KEEP_ALLELES:-0}"
REDCEA_VJ_BACKGROUND_ALLOW_SHORTFALL="${REDCEA_VJ_BACKGROUND_ALLOW_SHORTFALL:-0}"
TRB_BG_SOURCE_AIRR="${TRB_BG_SOURCE_AIRR:-redcea/data/backgrounds/trb_background_100k.tsv}"
TRB_BG_SOURCE_EMBEDDING="${TRB_BG_SOURCE_EMBEDDING:-redcea/data/backgrounds/trb_background_embeddings.parquet}"
TRA_BG_SOURCE_AIRR="${TRA_BG_SOURCE_AIRR:-$TRA_BG_AIRR}"
TRA_BG_SOURCE_EMBEDDING="${TRA_BG_SOURCE_EMBEDDING:-$TRA_BG_EMBEDDING}"
TRB_BG_AIRR="${TRB_BG_AIRR:-$REDCEA_PREP_ROOT/backgrounds/trb_background_vj.tsv}"
TRB_BG_EMBEDDING="${TRB_BG_EMBEDDING:-$REDCEA_PREP_ROOT/backgrounds/trb_background_vj_embeddings.parquet}"

cd "$ROOT_DIR"
mkdir -p logs

require_path() {
  local required_path="$1"
  if [[ ! -e "$required_path" ]]; then
    echo "Required input not found: $required_path" >&2
    exit 1
  fi
}

prepare_vj_background_if_needed() {
  local requested_chain="$1"
  local source_airr=""
  local source_embedding=""
  local output_airr=""
  local output_embedding=""
  local lock_dir=""
  local -a prep_args=()

  if [[ "$REDCEA_PREPARE_VJ_BACKGROUND" != "1" ]]; then
    return
  fi

  if [[ "$requested_chain" != "$REDCEA_VJ_BACKGROUND_CHAIN" ]]; then
    return
  fi

  case "$REDCEA_VJ_BACKGROUND_CHAIN" in
    TRA)
      source_airr="$TRA_BG_SOURCE_AIRR"
      source_embedding="$TRA_BG_SOURCE_EMBEDDING"
      output_airr="$TRA_BG_AIRR"
      output_embedding="$TRA_BG_EMBEDDING"
      ;;
    TRB)
      source_airr="$TRB_BG_SOURCE_AIRR"
      source_embedding="$TRB_BG_SOURCE_EMBEDDING"
      output_airr="$TRB_BG_AIRR"
      output_embedding="$TRB_BG_EMBEDDING"
      ;;
    *)
      echo "Unsupported REDCEA_VJ_BACKGROUND_CHAIN value: $REDCEA_VJ_BACKGROUND_CHAIN" >&2
      exit 1
      ;;
  esac

  if [[ -f "$output_airr" && -f "$output_embedding" ]]; then
    return
  fi

  require_path "$source_airr"
  require_path "$source_embedding"

  mkdir -p "$(dirname "$output_airr")" "$(dirname "$output_embedding")"
  lock_dir="${output_airr}.lockdir"
  while ! mkdir "$lock_dir" 2>/dev/null; do
    if [[ -f "$output_airr" && -f "$output_embedding" ]]; then
      return
    fi
    sleep 5
  done

  if [[ -f "$output_airr" && -f "$output_embedding" ]]; then
    rmdir "$lock_dir"
    return
  fi

  prep_args=(
    -u
    -m
    vdjdb_redcea.sample_beta_background_by_vj
    --chain "$REDCEA_VJ_BACKGROUND_CHAIN"
    --source "$source_airr"
    --source-embedding "$source_embedding"
    --target "$REDCEA_VDJDB"
    --output "$output_airr"
    --output-embedding "$output_embedding"
    --seed "$REDCEA_VJ_BACKGROUND_SEED"
  )
  if [[ "$REDCEA_VJ_BACKGROUND_KEEP_ALLELES" == "1" ]]; then
    prep_args+=(--keep-alleles)
  fi
  if [[ "$REDCEA_VJ_BACKGROUND_ALLOW_SHORTFALL" == "1" ]]; then
    prep_args+=(--allow-shortfall)
  fi

  python "${prep_args[@]}"
  rmdir "$lock_dir"
}

resolve_background_paths() {
  local chain="$1"
  case "$chain" in
    TRA)
      printf '%s\t%s\n' "$TRA_BG_AIRR" "$TRA_BG_EMBEDDING"
      ;;
    TRB)
      printf '%s\t%s\n' "$TRB_BG_AIRR" "$TRB_BG_EMBEDDING"
      ;;
    *)
      echo "Unsupported chain: $chain" >&2
      exit 1
      ;;
  esac
}

sanitize_token() {
  local raw="$1"
  local clean
  clean="$(printf '%s' "$raw" | tr -cd 'A-Za-z0-9_.-')"
  if [[ -z "$clean" ]]; then
    echo "unknown"
  else
    echo "$clean"
  fi
}

require_path "$REDCEA_EPITOPE_ROWS_FILE"
require_path "$REDCEA_VDJDB"

mapfile -t epitope_rows < <(grep -v $'^chain\t' "$REDCEA_EPITOPE_ROWS_FILE" | sed '/^[[:space:]]*$/d')

if [[ "${#epitope_rows[@]}" -eq 0 ]]; then
  echo "No epitope rows found in $REDCEA_EPITOPE_ROWS_FILE" >&2
  exit 1
fi

resolution_tag="${REDCEA_LEIDEN_RESOLUTION//./p}"
combo_tag="vdbscan_leiden_k${REDCEA_K_NEIGHBORS}_ek${REDCEA_EPS_K_NEIGHBORS}_res${resolution_tag}_min${REDCEA_CORE_MIN_SAMPLES}_eps${REDCEA_EPS_MODE}_sym${REDCEA_SYM_RULE}"
mkdir -p "$REDCEA_OUTPUT_ROOT" "$REDCEA_TCREMP_CACHE_ROOT" "$REDCEA_PREP_ROOT/backgrounds"

success_count=0
declare -a failures=()

for row in "${epitope_rows[@]}"; do
  IFS=$'\t' read -r chain epitope epitope_count <<< "$row"
  epitope_tag="$(sanitize_token "$epitope")"
  run_output_dir="${REDCEA_OUTPUT_ROOT}/${chain}/${epitope_tag}/${combo_tag}"
  shared_cache_dir="${REDCEA_TCREMP_CACHE_ROOT}/${chain}/${epitope_tag}"
  prefix="${chain,,}_vdjdb_${epitope_tag}"
  sample_airr="${shared_cache_dir}/${prefix}.tsv"
  sample_embedding="${shared_cache_dir}/${prefix}_sample_embeddings.parquet"

  mkdir -p "$run_output_dir" "$shared_cache_dir"
  prepare_vj_background_if_needed "$chain"

  IFS=$'\t' read -r background_airr background_embedding <<< "$(resolve_background_paths "$chain")"
  require_path "$background_airr"
  require_path "$background_embedding"

  export_args=(
    --vdjdb "$REDCEA_VDJDB"
    --epitope "$epitope"
    --chain "$chain"
    --species "$REDCEA_SPECIES"
    --output "$sample_airr"
  )
  if [[ "$REDCEA_INCLUDE_NONCANONICAL" == "1" ]]; then
    export_args+=(--include-noncanonical)
  fi

  if ! python scripts/export_vdjdb_epitope_airr.py "${export_args[@]}" >/dev/null; then
    failures+=("${chain}"$'\t'"${epitope}"$'\t'"${epitope_count}"$'\t'"export_airr"$'\t'"export_failed"$'\t'"${run_output_dir}")
    continue
  fi

  cmd=(
    python
    -u
    -m
    redcea.redcea
    --sample "$sample_airr"
    --background "$background_airr"
    --output "$run_output_dir"
    --prefix "$prefix"
    --chain "$chain"
    --species "$REDCEA_SPECIES"
    --metrics "$REDCEA_METRICS"
    --sample-embedding "$sample_embedding"
    --background-embedding "$background_embedding"
    --n-bg-points "$REDCEA_N_BG_POINTS"
    --cluster-pc-components "$REDCEA_CLUSTER_PC_COMPONENTS"
    --cluster-algo vdbscan_leiden
    --core-min-samples "$REDCEA_CORE_MIN_SAMPLES"
    --k-neighbors "$REDCEA_K_NEIGHBORS"
    --eps-k-neighbors "$REDCEA_EPS_K_NEIGHBORS"
    --leiden-resolution "$REDCEA_LEIDEN_RESOLUTION"
    --leiden-sub-resolution 1.0
    --eps-estimation-based-on "$REDCEA_EPS_MODE"
    --vdbscan-sym-rule "$REDCEA_SYM_RULE"
    --enrichment-test "$REDCEA_ENRICHMENT_TEST"
    --random-seed "$REDCEA_RANDOM_SEED"
    --nproc "${SLURM_CPUS_PER_TASK:-8}"
  )
  if [[ "$REDCEA_DEBUG_SAVE_INTERMEDIATE" == "1" ]]; then
    cmd+=(--debug-save-intermediate)
  fi
  if [[ "$REDCEA_ADD_AUXILIARY_CLUSTER_METRICS" == "1" ]]; then
    cmd+=(--add-auxiliary-cluster-metrics)
  fi

  echo "Running combo=${combo_tag}: chain=${chain}, epitope=${epitope}, count=${epitope_count}"
  if ! "${cmd[@]}"; then
    failures+=("${chain}"$'\t'"${epitope}"$'\t'"${epitope_count}"$'\t'"redcea"$'\t'"run_failed"$'\t'"${run_output_dir}")
    continue
  fi

  success_count=$((success_count + 1))
  cat > "${run_output_dir}/run_metadata.tsv" <<EOF
key	value
chain	${chain}
epitope	${epitope}
epitope_count	${epitope_count}
k_neighbors	${REDCEA_K_NEIGHBORS}
eps_k_neighbors	${REDCEA_EPS_K_NEIGHBORS}
core_min_samples	${REDCEA_CORE_MIN_SAMPLES}
leiden_resolution	${REDCEA_LEIDEN_RESOLUTION}
leiden_sub_resolution	1.0
eps_estimation_based_on	${REDCEA_EPS_MODE}
vdbscan_sym_rule	${REDCEA_SYM_RULE}
cluster_algo	vdbscan_leiden
output_tag	${combo_tag}
output_dir	${run_output_dir}
tcremp_cache_dir	${shared_cache_dir}
sample_airr	${sample_airr}
sample_embedding	${sample_embedding}
background_airr	${background_airr}
background_embedding	${background_embedding}
slurm_job_id	${SLURM_JOB_ID:-}
EOF
done

if [[ "${#failures[@]}" -gt 0 ]]; then
  failure_report="${REDCEA_OUTPUT_ROOT}/${combo_tag}_failed_epitopes.tsv"
  {
    printf 'chain\tepitope\tepitope_count\tstage\tstatus\toutput_dir\n'
    printf '%s\n' "${failures[@]}"
  } > "$failure_report"
  echo "Wrote failure report to ${failure_report}"
fi

if (( success_count == 0 )); then
  echo "All epitopes failed for combo ${combo_tag}" >&2
  exit 1
fi

echo "Completed combo ${combo_tag}: successful_epitopes=${success_count}, failed_epitopes=${#failures[@]}"
