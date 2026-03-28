#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

REDCEA_VDJDB="redcea/data/vdjdb_full.txt"
REDCEA_SPECIES="HomoSapiens"
REDCEA_NPROC="8"
REDCEA_CONDA_ENV="vdjdb-redcea"

TRA_BG_AIRR="redcea/data/backgrounds/tra_background_100k.tsv"
TRA_BG_EMBEDDING="redcea/data/backgrounds/tra_background_embeddings.parquet"
TRB_BG_AIRR="redcea/data/backgrounds/trb_background_100k.tsv"
TRB_BG_EMBEDDING="redcea/data/backgrounds/trb_background_embeddings.parquet"

REDCEA_OUTPUT="results/redcea"

mkdir -p "$REDCEA_OUTPUT"

run_chain() {
  local chain="$1"
  local bg_airr="$2"
  local bg_embedding="$3"
  local output_dir="$4"

  conda run -n "$REDCEA_CONDA_ENV" python -m vdjdb_redcea.vdjdb_clusters_launch_with_transform \
    --vdjdb "$REDCEA_VDJDB" \
    --background-airr "$bg_airr" \
    --background-embedding "$bg_embedding" \
    --output "$output_dir" \
    --chain "$chain" \
    --species "$REDCEA_SPECIES" \
    --nproc "$REDCEA_NPROC"
}

run_chain TRA "$TRA_BG_AIRR" "$TRA_BG_EMBEDDING" "$REDCEA_OUTPUT" &
pid_tra=$!

run_chain TRB "$TRB_BG_AIRR" "$TRB_BG_EMBEDDING" "$REDCEA_OUTPUT" &
pid_trb=$!

wait "$pid_tra"
wait "$pid_trb"
