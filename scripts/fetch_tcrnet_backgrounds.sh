#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

URL="https://zenodo.org/record/6339774/files/pools.zip"
OUTPUT_DIR="tcrnet"

mkdir -p "$OUTPUT_DIR"

ARCHIVE_PATH="$OUTPUT_DIR/pools.zip"

echo "Downloading TCRNET background pools from $URL"
curl -L "$URL" -o "$ARCHIVE_PATH"

echo "Extracting $ARCHIVE_PATH"
unzip -o "$ARCHIVE_PATH" -d "$OUTPUT_DIR"
