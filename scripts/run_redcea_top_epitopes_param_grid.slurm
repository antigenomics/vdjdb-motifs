#!/bin/bash
#SBATCH --job-name=redcea-top-epi-grid
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --partition=short
#SBATCH --output=logs/redcea_top_epi_grid_%A_%a.out
#SBATCH --error=logs/redcea_top_epi_grid_%A_%a.err

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-/home/evlasova/vdjdb-motifs}"
REDCEA_VDJDB="${REDCEA_VDJDB:-vdjdb_release/vdjdb.slim.txt}"
REDCEA_SPECIES="${REDCEA_SPECIES:-HomoSapiens}"
REDCEA_EPITOPE_ROWS_FILE="${REDCEA_EPITOPE_ROWS_FILE:-}"
REDCEA_OUTPUT_ROOT="${REDCEA_OUTPUT_ROOT:-results/redcea_vdbscan_leiden_top_epitope_grid}"
REDCEA_TCREMP_CACHE_ROOT="${REDCEA_TCREMP_CACHE_ROOT:-results/redcea_vdbscan_leiden_top_epitope_grid_cache}"
REDCEA_PREP_ROOT="${REDCEA_PREP_ROOT:-results/redcea_vdbscan_leiden_top_epitope_grid_prep}"
REDCEA_INCLUDE_NONCANONICAL="${REDCEA_INCLUDE_NONCANONICAL:-0}"
REDCEA_METRICS="${REDCEA_METRICS:-dissimilarity}"
REDCEA_CLUSTER_PC_COMPONENTS="${REDCEA_CLUSTER_PC_COMPONENTS:-50}"
REDCEA_N_BG_POINTS="${REDCEA_N_BG_POINTS:-100000}"
REDCEA_ENRICHMENT_TEST="${REDCEA_ENRICHMENT_TEST:-zbinom}"
REDCEA_RANDOM_SEED="${REDCEA_RANDOM_SEED:-42}"
REDCEA_DEBUG_SAVE_INTERMEDIATE="${REDCEA_DEBUG_SAVE_INTERMEDIATE:-0}"
REDCEA_ADD_AUXILIARY_CLUSTER_METRICS="${REDCEA_ADD_AUXILIARY_CLUSTER_METRICS:-0}"
REDCEA_VDBSCAN_LEIDEN_K_GRID="${REDCEA_VDBSCAN_LEIDEN_K_GRID:-4 8 12}"
REDCEA_VDBSCAN_LEIDEN_EPS_K_GRID="${REDCEA_VDBSCAN_LEIDEN_EPS_K_GRID:-4 8 12}"
REDCEA_VDBSCAN_LEIDEN_CORE_MIN_SAMPLES_GRID="${REDCEA_VDBSCAN_LEIDEN_CORE_MIN_SAMPLES_GRID:-3}"
REDCEA_VDBSCAN_LEIDEN_RESOLUTION_GRID="${REDCEA_VDBSCAN_LEIDEN_RESOLUTION_GRID:-0.5 1.0 1.5}"
REDCEA_VDBSCAN_LEIDEN_EPS_MODE_GRID="${REDCEA_VDBSCAN_LEIDEN_EPS_MODE_GRID:-sample all}"
REDCEA_VDBSCAN_LEIDEN_SYM_RULE_GRID="${REDCEA_VDBSCAN_LEIDEN_SYM_RULE_GRID:-asymmetric min}"

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

read -r -a k_grid <<< "$REDCEA_VDBSCAN_LEIDEN_K_GRID"
read -r -a eps_k_grid <<< "$REDCEA_VDBSCAN_LEIDEN_EPS_K_GRID"
read -r -a min_samples_grid <<< "$REDCEA_VDBSCAN_LEIDEN_CORE_MIN_SAMPLES_GRID"
read -r -a resolution_grid <<< "$REDCEA_VDBSCAN_LEIDEN_RESOLUTION_GRID"
read -r -a eps_mode_grid <<< "$REDCEA_VDBSCAN_LEIDEN_EPS_MODE_GRID"
read -r -a sym_rule_grid <<< "$REDCEA_VDBSCAN_LEIDEN_SYM_RULE_GRID"

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

load_epitope_rows() {
  if [[ -n "$REDCEA_EPITOPE_ROWS_FILE" ]]; then
    require_path "$REDCEA_EPITOPE_ROWS_FILE"
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
  echo "No epitopes were configured for the grid" >&2
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
task_id="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

if (( task_id < 0 || task_id >= total_tasks )); then
  echo "Task ${task_id} is outside grid size ${total_tasks}; exiting without work."
  exit 0
fi

epitope_index=$((task_id % n_epitopes))
k_index=$(((task_id / n_epitopes) % n_k))
eps_k_index=$(((task_id / (n_epitopes * n_k)) % n_eps_k))
min_samples_index=$(((task_id / (n_epitopes * n_k * n_eps_k)) % n_min_samples))
resolution_index=$(((task_id / (n_epitopes * n_k * n_eps_k * n_min_samples)) % n_resolution))
eps_mode_index=$(((task_id / (n_epitopes * n_k * n_eps_k * n_min_samples * n_resolution)) % n_eps_mode))
sym_rule_index=$((task_id / (n_epitopes * n_k * n_eps_k * n_min_samples * n_resolution * n_eps_mode)))

row="${epitope_rows[$epitope_index]}"
IFS=$'\t' read -r chain epitope epitope_count <<< "$row"
k_neighbors="${k_grid[$k_index]}"
eps_k_neighbors="${eps_k_grid[$eps_k_index]}"
core_min_samples="${min_samples_grid[$min_samples_index]}"
leiden_resolution="${resolution_grid[$resolution_index]}"
eps_mode="${eps_mode_grid[$eps_mode_index]}"
sym_rule="${sym_rule_grid[$sym_rule_index]}"

epitope_tag="$(sanitize_token "$epitope")"
resolution_tag="${leiden_resolution//./p}"
output_tag="vdbscan_leiden_k${k_neighbors}_ek${eps_k_neighbors}_res${resolution_tag}_min${core_min_samples}_eps${eps_mode}_sym${sym_rule}"
run_output_dir="${REDCEA_OUTPUT_ROOT}/${chain}/${epitope_tag}/${output_tag}"
shared_cache_dir="${REDCEA_TCREMP_CACHE_ROOT}/${chain}/${epitope_tag}"
prefix="${chain,,}_vdjdb_${epitope_tag}"
sample_airr="${shared_cache_dir}/${prefix}.tsv"
sample_embedding="${shared_cache_dir}/${prefix}_sample_embeddings.parquet"
mkdir -p "$run_output_dir" "$shared_cache_dir" "$REDCEA_PREP_ROOT/backgrounds"

prepare_vj_background_if_needed "$chain"

IFS=$'\t' read -r background_airr background_embedding <<< "$(resolve_background_paths "$chain")"
require_path "$background_airr"
require_path "$background_embedding"
require_path "$REDCEA_VDJDB"

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
python scripts/export_vdjdb_epitope_airr.py "${export_args[@]}" >/dev/null

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
  --core-min-samples "$core_min_samples"
  --k-neighbors "$k_neighbors"
  --eps-k-neighbors "$eps_k_neighbors"
  --leiden-resolution "$leiden_resolution"
  --leiden-sub-resolution 1.0
  --eps-estimation-based-on "$eps_mode"
  --vdbscan-sym-rule "$sym_rule"
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

echo "Running task ${task_id}/${total_tasks}: chain=${chain}, epitope=${epitope}, count=${epitope_count}, k=${k_neighbors}, eps_k=${eps_k_neighbors}, min_samples=${core_min_samples}, res=${leiden_resolution}, eps_mode=${eps_mode}, sym_rule=${sym_rule}"
echo "Output dir: ${run_output_dir}"
echo "Sample AIRR: ${sample_airr}"
echo "Sample embedding cache: ${sample_embedding}"
echo "Background AIRR: ${background_airr}"
echo "Background embedding: ${background_embedding}"

"${cmd[@]}"

cat > "${run_output_dir}/run_metadata.tsv" <<EOF
key	value
chain	${chain}
epitope	${epitope}
epitope_count	${epitope_count}
k_neighbors	${k_neighbors}
eps_k_neighbors	${eps_k_neighbors}
core_min_samples	${core_min_samples}
leiden_resolution	${leiden_resolution}
leiden_sub_resolution	1.0
eps_estimation_based_on	${eps_mode}
vdbscan_sym_rule	${sym_rule}
cluster_algo	vdbscan_leiden
output_tag	${output_tag}
output_dir	${run_output_dir}
tcremp_cache_dir	${shared_cache_dir}
sample_airr	${sample_airr}
sample_embedding	${sample_embedding}
background_airr	${background_airr}
background_embedding	${background_embedding}
slurm_array_task_id	${task_id}
EOF
