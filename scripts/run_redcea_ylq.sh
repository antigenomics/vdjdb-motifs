#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

REDCEA_VDJDB="${REDCEA_VDJDB:-vdjdb_release/vdjdb.slim.txt}"
REDCEA_SPECIES="${REDCEA_SPECIES:-HomoSapiens}"
REDCEA_NPROC="${REDCEA_NPROC:-8}"
REDCEA_CHAIN="${REDCEA_CHAIN:-TRB}"
REDCEA_OUTPUT="${REDCEA_OUTPUT:-results/redcea_ylq_glc}"
REDCEA_EPITOPES="${REDCEA_EPITOPES:-${REDCEA_EPITOPE:-GLCTLVAML YLQPRTFLL}}"
REDCEA_MIN_EPITOPE_CLONOTYPES="${REDCEA_MIN_EPITOPE_CLONOTYPES:-1}"
REDCEA_SKIP_UMAP="${REDCEA_SKIP_UMAP:-1}"
REDCEA_K_NEIGHBORS="${REDCEA_K_NEIGHBORS:-15}"
REDCEA_LEIDEN_RESOLUTION="${REDCEA_LEIDEN_RESOLUTION:-1.0}"
REDCEA_OUTPUT_TAG="${REDCEA_OUTPUT_TAG:-}"
REDCEA_TCREMP_CACHE_DIR="${REDCEA_TCREMP_CACHE_DIR:-}"

# Plotting UMAP defaults are kept for compatibility, but can be skipped entirely
# with REDCEA_SKIP_UMAP=1 for fast tuning runs.
REDCEA_UMAP_N_NEIGHBORS_TRA="${REDCEA_UMAP_N_NEIGHBORS_TRA:-25}"
REDCEA_UMAP_MIN_DIST_TRA="${REDCEA_UMAP_MIN_DIST_TRA:-0.4}"
REDCEA_UMAP_N_NEIGHBORS_TRB="${REDCEA_UMAP_N_NEIGHBORS_TRB:-50}"
REDCEA_UMAP_MIN_DIST_TRB="${REDCEA_UMAP_MIN_DIST_TRB:-0.4}"

TRA_BG_AIRR="${TRA_BG_AIRR:-redcea/data/backgrounds/tra_background_100k.tsv}"
TRA_BG_EMBEDDING="${TRA_BG_EMBEDDING:-redcea/data/backgrounds/tra_background_embeddings.parquet}"
TRB_BG_AIRR="${TRB_BG_AIRR:-redcea/data/backgrounds/trb_background_100k.tsv}"
TRB_BG_EMBEDDING="${TRB_BG_EMBEDDING:-redcea/data/backgrounds/trb_background_embeddings.parquet}"

mkdir -p "$REDCEA_OUTPUT"

case "$REDCEA_CHAIN" in
  TRA|TRB)
    ;;
  *)
    echo "Unsupported REDCEA_CHAIN value: $REDCEA_CHAIN (expected TRA or TRB)" >&2
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

activate_redcea_env() {
  if ! command -v conda >/dev/null 2>&1; then
    echo "conda is required but was not found in PATH" >&2
    exit 1
  fi

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
  local -a extra_args=()

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

  if [[ "$REDCEA_SKIP_UMAP" == "1" ]]; then
    extra_args+=(--skip-umap)
  fi

  read -r -a epitope_args <<< "$REDCEA_EPITOPES"
  local output_tag="$REDCEA_OUTPUT_TAG"
  if [[ -z "$output_tag" ]]; then
    output_tag="k${REDCEA_K_NEIGHBORS}_res${REDCEA_LEIDEN_RESOLUTION//./p}"
  fi
  if [[ -n "$REDCEA_TCREMP_CACHE_DIR" ]]; then
    extra_args+=(--tcremp-cache-dir "$REDCEA_TCREMP_CACHE_DIR")
  fi

  python -u -m vdjdb_redcea.vdjdb_clusters_launch_with_transform \
    --vdjdb "$REDCEA_VDJDB" \
    --background-airr "$bg_airr" \
    --background-embedding "$bg_embedding" \
    --output "$output_dir" \
    --chain "$chain" \
    --species "$REDCEA_SPECIES" \
    --epitopes "${epitope_args[@]}" \
    --min-epitope-clonotypes "$REDCEA_MIN_EPITOPE_CLONOTYPES" \
    --k-neighbors "$REDCEA_K_NEIGHBORS" \
    --leiden-resolution "$REDCEA_LEIDEN_RESOLUTION" \
    --output-tag "$output_tag" \
    --umap-n-neighbors "$umap_n_neighbors" \
    --umap-min-dist "$umap_min_dist" \
    --nproc "$REDCEA_NPROC" \
    "${extra_args[@]}"
}

activate_redcea_env

if [[ "$REDCEA_CHAIN" == "TRA" ]]; then
  require_path "$TRA_BG_AIRR"
  require_path "$TRA_BG_EMBEDDING"
  run_chain TRA "$TRA_BG_AIRR" "$TRA_BG_EMBEDDING" "$REDCEA_OUTPUT"
else
  require_path "$TRB_BG_AIRR"
  require_path "$TRB_BG_EMBEDDING"
  run_chain TRB "$TRB_BG_AIRR" "$TRB_BG_EMBEDDING" "$REDCEA_OUTPUT"
fi
