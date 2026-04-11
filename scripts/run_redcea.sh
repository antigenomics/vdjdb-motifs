#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

REDCEA_VDJDB="${REDCEA_VDJDB:-vdjdb_release/vdjdb.slim.txt}"
REDCEA_SPECIES="${REDCEA_SPECIES:-HomoSapiens}"
REDCEA_NPROC="${REDCEA_NPROC:-8}"
REDCEA_CONDA_ENV="${REDCEA_CONDA_ENV:-vdjdb-redcea}"
REDCEA_MIN_EPITOPE_CLONOTYPES="${REDCEA_MIN_EPITOPE_CLONOTYPES:-100}"
REDCEA_CHAIN="${REDCEA_CHAIN:-both}"

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

run_chain() {
  local chain="$1"
  local bg_airr="$2"
  local bg_embedding="$3"
  local output_dir="$4"
  local log_path="$output_dir/${chain,,}.log"

  conda run -n "$REDCEA_CONDA_ENV" python -u -m vdjdb_redcea.vdjdb_clusters_launch_with_transform \
    --vdjdb "$REDCEA_VDJDB" \
    --background-airr "$bg_airr" \
    --background-embedding "$bg_embedding" \
    --output "$output_dir" \
    --chain "$chain" \
    --species "$REDCEA_SPECIES" \
    --min-epitope-clonotypes "$REDCEA_MIN_EPITOPE_CLONOTYPES" \
    --nproc "$REDCEA_NPROC" \
    2>&1 | tee "$log_path"
}

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
