#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

URL="${1:-https://zenodo.org/record/19262060/files/redcea_bg.gz}"
OUTPUT_DIR="redcea/data/backgrounds"
ARCHIVE_PATH="$OUTPUT_DIR/$(basename "$URL")"
UNPACK_DIR="$OUTPUT_DIR/redcea_bundle"

mkdir -p "$OUTPUT_DIR"

required_backgrounds=(
  "$OUTPUT_DIR/tra_background_100k.tsv"
  "$OUTPUT_DIR/tra_background_embeddings.parquet"
  "$OUTPUT_DIR/trb_background_100k.tsv"
  "$OUTPUT_DIR/trb_background_embeddings.parquet"
)

all_backgrounds_present=true
for required_path in "${required_backgrounds[@]}"; do
  if [[ ! -f "$required_path" ]]; then
    all_backgrounds_present=false
    break
  fi
done

if [[ "$all_backgrounds_present" == true ]]; then
  echo "REDCEA background files already exist in $OUTPUT_DIR, skipping download"
  exit 0
fi

rm -rf "$UNPACK_DIR"
mkdir -p "$UNPACK_DIR"

if [[ -f "$ARCHIVE_PATH" ]]; then
  echo "Using existing REDCEA archive $ARCHIVE_PATH"
else
  echo "Downloading REDCEA bundle from $URL"
  curl -L "$URL" -o "$ARCHIVE_PATH"
fi

echo "Extracting $ARCHIVE_PATH"
tar -xzf "$ARCHIVE_PATH" -C "$UNPACK_DIR"

SOURCE_DIR="$UNPACK_DIR"
if [[ -d "$UNPACK_DIR/download" ]]; then
  SOURCE_DIR="$UNPACK_DIR/download"
fi

echo "Installing REDCEA background files"

cp "$SOURCE_DIR/tra_background_100k.tsv" "$OUTPUT_DIR/tra_background_100k.tsv"
cp "$SOURCE_DIR/tra_background_100k_embeddings.parquet" "$OUTPUT_DIR/tra_background_embeddings.parquet"
cp "$SOURCE_DIR/trb_background_100k.tsv" "$OUTPUT_DIR/trb_background_100k.tsv"
cp "$SOURCE_DIR/trb_background_100k_embeddings.parquet" "$OUTPUT_DIR/trb_background_embeddings.parquet"

echo "REDCEA background files are ready in $OUTPUT_DIR"
