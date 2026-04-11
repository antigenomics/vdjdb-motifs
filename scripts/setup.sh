#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

TARGET="${1:-all}"
REDCEA_BG_URL="${2:-}"

usage() {
  cat <<'EOF'
Usage:
  ./scripts/setup.sh [tcrnet|redcea|all] [redcea_background_url]

Examples:
  ./scripts/setup.sh tcrnet
  ./scripts/setup.sh redcea
  ./scripts/setup.sh all
  ./scripts/setup.sh redcea https://zenodo.org/record/19262060/files/redcea_bg.gz
EOF
}

bootstrap_tcrnet() {
  Rscript scripts/install_tcrnet_deps.R
  ./scripts/fetch_tcrnet_backgrounds.sh
}

bootstrap_redcea() {
  ./scripts/install_redcea.sh
  if [[ -n "$REDCEA_BG_URL" ]]; then
    ./scripts/fetch_redcea_backgrounds.sh "$REDCEA_BG_URL"
  else
    ./scripts/fetch_redcea_backgrounds.sh
  fi
}

case "$TARGET" in
  tcrnet)
    bootstrap_tcrnet
    ;;
  redcea)
    bootstrap_redcea
    ;;
  all)
    bootstrap_tcrnet
    bootstrap_redcea
    ;;
  -h|--help|help)
    usage
    ;;
  *)
    echo "Unknown bootstrap target: $TARGET" >&2
    usage >&2
    exit 1
    ;;
esac
