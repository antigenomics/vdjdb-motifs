#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

REDCEA_VDJDB="${REDCEA_VDJDB:-vdjdb_release/vdjdb.slim.txt}"
REDCEA_SPECIES="${REDCEA_SPECIES:-HomoSapiens}"
REDCEA_NPROC="${REDCEA_NPROC:-8}"
REDCEA_MIN_EPITOPE_CLONOTYPES="${REDCEA_MIN_EPITOPE_CLONOTYPES:-100}"
REDCEA_CHAIN="${REDCEA_CHAIN:-both}"
REDCEA_K_NEIGHBORS="${REDCEA_K_NEIGHBORS:-15}"
REDCEA_LEIDEN_RESOLUTION="${REDCEA_LEIDEN_RESOLUTION:-1.0}"
REDCEA_CLUSTER_MIN_SAMPLES="${REDCEA_CLUSTER_MIN_SAMPLES:-5}"

# Plotting UMAP defaults can differ per chain. If a global REDCEA_UMAP_* override
# is set in the environment, it still wins over the chain-specific defaults.
REDCEA_UMAP_N_NEIGHBORS_TRA="${REDCEA_UMAP_N_NEIGHBORS_TRA:-25}"
REDCEA_UMAP_MIN_DIST_TRA="${REDCEA_UMAP_MIN_DIST_TRA:-0.4}"
REDCEA_UMAP_N_NEIGHBORS_TRB="${REDCEA_UMAP_N_NEIGHBORS_TRB:-50}"
REDCEA_UMAP_MIN_DIST_TRB="${REDCEA_UMAP_MIN_DIST_TRB:-0.4}"

TRA_BG_AIRR="${TRA_BG_AIRR:-redcea/data/backgrounds/tra_background_100k.tsv}"
TRA_BG_EMBEDDING="${TRA_BG_EMBEDDING:-redcea/data/backgrounds/tra_background_embeddings.parquet}"
TRB_BG_AIRR="${TRB_BG_AIRR:-redcea/data/backgrounds/trb_background_100k.tsv}"
TRB_BG_EMBEDDING="${TRB_BG_EMBEDDING:-redcea/data/backgrounds/trb_background_embeddings.parquet}"

REDCEA_OUTPUT="${REDCEA_OUTPUT:-results/redcea}"

mkdir -p "$REDCEA_OUTPUT"

case "$REDCEA_CHAIN" in
  TRA|TRB|both)
    ;;
  *)
    echo "Unsupported REDCEA_CHAIN value: $REDCEA_CHAIN (expected TRA, TRB, or both)" >&2
    exit 1
    ;;
esac

require_path() {
  local required_path="$1"
  if [[ ! -e "$required_path" ]]; then
    echo "Required input not found: $required_path" >&2
    exit 1
  fi
}

require_path "$REDCEA_VDJDB"

case "$REDCEA_CHAIN" in
  TRA)
    require_path "$TRA_BG_AIRR"
    require_path "$TRA_BG_EMBEDDING"
    ;;
  TRB)
    require_path "$TRB_BG_AIRR"
    require_path "$TRB_BG_EMBEDDING"
    ;;
  both)
    require_path "$TRA_BG_AIRR"
    require_path "$TRA_BG_EMBEDDING"
    require_path "$TRB_BG_AIRR"
    require_path "$TRB_BG_EMBEDDING"
    ;;
esac

activate_redcea_env() {
  if ! command -v conda >/dev/null 2>&1; then
    echo "conda is required but was not found in PATH" >&2
    exit 1
  fi

  # Initialize conda for this non-interactive bash process, then activate once.
  eval "$(conda shell.bash hook)"
  conda activate vdjdb-redcea
}

run_chain() {
  local chain="$1"
  local bg_airr="$2"
  local bg_embedding="$3"
  local output_dir="$4"
  local umap_n_neighbors
  local umap_min_dist

  case "$chain" in
    TRA)
      umap_n_neighbors="${REDCEA_UMAP_N_NEIGHBORS:-$REDCEA_UMAP_N_NEIGHBORS_TRA}"
      umap_min_dist="${REDCEA_UMAP_MIN_DIST:-$REDCEA_UMAP_MIN_DIST_TRA}"
      ;;
    TRB)
      umap_n_neighbors="${REDCEA_UMAP_N_NEIGHBORS:-$REDCEA_UMAP_N_NEIGHBORS_TRB}"
      umap_min_dist="${REDCEA_UMAP_MIN_DIST:-$REDCEA_UMAP_MIN_DIST_TRB}"
      ;;
    *)
      echo "Unsupported chain value for UMAP config: $chain" >&2
      exit 1
      ;;
  esac

  python -u -m vdjdb_redcea.vdjdb_clusters_launch_with_transform \
    --vdjdb "$REDCEA_VDJDB" \
    --background-airr "$bg_airr" \
    --background-embedding "$bg_embedding" \
    --output "$output_dir" \
    --chain "$chain" \
    --species "$REDCEA_SPECIES" \
    --min-epitope-clonotypes "$REDCEA_MIN_EPITOPE_CLONOTYPES" \
    --k-neighbors "$REDCEA_K_NEIGHBORS" \
    --leiden-resolution "$REDCEA_LEIDEN_RESOLUTION" \
    --cluster-min-samples "$REDCEA_CLUSTER_MIN_SAMPLES" \
    --umap-n-neighbors "$umap_n_neighbors" \
    --umap-min-dist "$umap_min_dist" \
    --nproc "$REDCEA_NPROC"
}

activate_redcea_env

if [[ "$REDCEA_CHAIN" == "TRA" ]]; then
  run_chain TRA "$TRA_BG_AIRR" "$TRA_BG_EMBEDDING" "$REDCEA_OUTPUT"
elif [[ "$REDCEA_CHAIN" == "TRB" ]]; then
  run_chain TRB "$TRB_BG_AIRR" "$TRB_BG_EMBEDDING" "$REDCEA_OUTPUT"
else
  run_chain TRA "$TRA_BG_AIRR" "$TRA_BG_EMBEDDING" "$REDCEA_OUTPUT" &
  pid_tra=$!

  run_chain TRB "$TRB_BG_AIRR" "$TRB_BG_EMBEDDING" "$REDCEA_OUTPUT" &
  pid_trb=$!

  wait "$pid_tra"
  wait "$pid_trb"
fi
