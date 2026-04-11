# VDJdb Motif Clustering

This repository now keeps two clustering pipelines side by side:

- `tcrnet/`: the R Markdown workflow based on VDJtools/TCRNET
- `redcea/`: the new Python package based on TCRemP + REDCEA-style shared-background clustering

The goal is to keep the TCRNET pipeline runnable as before while letting the new implementation evolve independently.

## Repository Layout

```text
vdjdb-motifs/
  tcrnet/   TCRNET clustering workflow
  redcea/   installable Python package
  results/  generated tabular outputs
  figures/  generated PDF figures
  scripts/  run and setup helpers
```

## Simple Commands

If you just want one command per pipeline:

```bash
./scripts/setup.sh tcrnet
./scripts/run_tcrnet.sh

./scripts/setup.sh redcea
./scripts/run_redcea.sh
```

`run_redcea.sh` launches `TRA` and `TRB` in parallel.

Thread count is configurable for both pipelines:

```bash
TCRNET_NPROC=16 ./scripts/run_tcrnet.sh
REDCEA_NPROC=16 ./scripts/run_redcea.sh
```

If you want to prepare everything in one go:

```bash
./scripts/setup.sh all
```

Assumptions:

- `vdjdb_release/vdjdb.slim.txt` already exists before running either pipeline.
- `VDJtools` is already installed and available to TCRNET.
- the repository is run in a Unix-like environment.

## Background Data

- `tcrnet` backgrounds are public and downloaded automatically from Zenodo by `scripts/fetch_tcrnet_backgrounds.sh`.
- `redcea` needs two prepared background files:
  - AIRR repertoire table
  - embedding parquet file
- The repository now includes a downloader script for `redcea`.
- By default it downloads `https://zenodo.org/record/19262060/files/redcea_bg.gz`.
- You can still pass a different Zenodo archive URL as the first argument if needed.
- `scripts/setup.sh redcea <url>` passes that URL through to the REDCEA background fetch step.
- `scripts/install_redcea.sh` creates a separate conda environment and installs this repository's `redcea` package together with its declared git-based dependencies.

## REDCEA Zenodo Files

If you want `redcea` to start quickly from downloaded assets, the default Zenodo bundle is:

- `https://zenodo.org/record/19262060/files/redcea_bg.gz`

Required per chain:
- `tra_background_100k.tsv`
- `tra_background_100k_embeddings.parquet`
- `trb_background_100k.tsv`
- `trb_background_100k_embeddings.parquet`

What these do:
- `*_background_100k.tsv` is the actual background repertoire input passed to `--background-airr`.
- `*_background_100k_embeddings.parquet` is the precomputed background embedding input passed to `--background-embedding` after renaming during install.

What you do not need to upload:
- background transform caches
- per-epitope sample embeddings
- per-run cluster tables
- HTML visualizations
- `cluster_members_*.txt`

Practical recommendation:
- if you run both `TRA` and `TRB`, include both required chain-specific AIRR/embedding pairs in the bundle
- `./scripts/run_redcea.sh` uses the default local paths above directly

UMAP tuning note:
- the REDCEA plotting layout now accepts `--umap-n-neighbors` and `--umap-min-dist`
- `scripts/run_redcea.sh` keeps separate per-chain defaults via `REDCEA_UMAP_N_NEIGHBORS_TRA` / `REDCEA_UMAP_MIN_DIST_TRA` and `REDCEA_UMAP_N_NEIGHBORS_TRB` / `REDCEA_UMAP_MIN_DIST_TRB`
- if needed, a global `REDCEA_UMAP_N_NEIGHBORS` or `REDCEA_UMAP_MIN_DIST` environment override still takes precedence for ad hoc runs
- if sample clusters form detached islands far from the grey background cloud, try increasing these values, for example `--umap-n-neighbors 50 --umap-min-dist 0.4`
- changing either value invalidates the previous plotting transform on purpose, so the cached background transform is recomputed automatically

## Notes

- TCRNET outputs now go to `results/tcrnet/` and keep the standard names `cluster_members.txt` and `motif_pwms.txt`.
- TCRNET PDF figures go to the repository-level `figures/` directory.
- TCRNET respects `TCRNET_NPROC` for parallel R steps.
- REDCEA writes both chains into `results/redcea/`.
- REDCEA HTML visualizations are collected in `results/redcea/viz/`.
- REDCEA `cluster_members_TRA.txt` and `cluster_members_TRB.txt` are written directly into `results/redcea/`.
- REDCEA uses `vdjdb_release/vdjdb.slim.txt` as the default VDJdb input table.
- REDCEA respects `REDCEA_NPROC`; when `REDCEA_CHAIN=both`, `TRA` and `TRB` are still launched in parallel as separate jobs.
