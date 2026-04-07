#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

URL="${1:-https://zenodo.org/record/6339774/files/redcea_bg.gz}"
OUTPUT_DIR="redcea/data/backgrounds"
ARCHIVE_PATH="$OUTPUT_DIR/$(basename "$URL")"
UNPACK_DIR="$OUTPUT_DIR/redcea_bundle"
TCREMP_DIR="results/redcea/tcremp"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$TCREMP_DIR"

rm -rf "$UNPACK_DIR"

echo "Downloading REDCEA bundle from $URL"
curl -L "$URL" -o "$ARCHIVE_PATH"

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

cp "$SOURCE_DIR/tra_background_transform.joblib" "$TCREMP_DIR/tra_background_transform.joblib"
cp "$SOURCE_DIR/trb_background_transform.joblib" "$TCREMP_DIR/trb_background_transform.joblib"

cp "$SOURCE_DIR"/tra_background_transform_bg_umap_*.npy "$TCREMP_DIR/"
cp "$SOURCE_DIR"/trb_background_transform_bg_umap_*.npy "$TCREMP_DIR/"

echo "REDCEA background files are ready in $OUTPUT_DIR"
