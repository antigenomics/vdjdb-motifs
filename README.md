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
`setup.sh redcea` requires a working `conda` installation because it creates and updates the `vdjdb-redcea` environment automatically.

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
- the repository is run in a Unix-like environment.

Required command-line tools:

- shared: `bash`
- TCRNET: `Rscript`, `python`, `VDJtools`, `wget`, `unzip`
- REDCEA: `conda`, `curl`, `tar`

Notes on these requirements:

- `scripts/setup.sh tcrnet` installs the needed R packages automatically, but it does not install `VDJtools` itself.
- `scripts/run_tcrnet.sh` now also runs a small Python post-processing step that generates per-epitope HTML visualizations.
- `scripts/setup.sh redcea` creates and updates the `vdjdb-redcea` conda environment automatically.

## Background Data

- `tcrnet` backgrounds are public and downloaded automatically from Zenodo by `scripts/fetch_tcrnet_backgrounds.sh`.
- `redcea` needs two prepared background files:
  - AIRR repertoire table
  - embedding parquet file
- The repository now includes a downloader script for `redcea`.
- By default it downloads `https://zenodo.org/records/19520535/files/redcea_bg.gz`.
- You can still pass a different Zenodo archive URL as the first argument if needed.
- `scripts/setup.sh redcea <url>` passes that URL through to the REDCEA background fetch step.
- `scripts/install_redcea.sh` creates a separate conda environment and installs this repository's `redcea` package together with its declared git-based dependencies.

## REDCEA Zenodo Files

If you want `redcea` to start quickly from downloaded assets, the default Zenodo bundle is:

- `https://zenodo.org/records/19520535/files/redcea_bg.gz`

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

Per-epitope clustering note:
- `python -m vdjdb_redcea.vdjdb_epitope_clustering` now accepts `--epitope-config <json>`
- the JSON can override `cluster_algo`, `k_neighbors`, `eps_k_neighbors`, `leiden_resolution`, `cluster_min_samples`, `eps_estimation_based_on`, and `vdbscan_sym_rule` per epitope while still fitting one shared background transform and one joint plotting UMAP across the whole selected epitope set
- this is useful when you want one common background view for multiple epitopes but still keep each epitope's best clustering hyperparameters
- accepted JSON shapes are either a top-level epitope mapping or `{ "epitopes": { ... } }`

Example:

```json
{
  "epitopes": {
    "YLQPRTFLL": {
      "cluster_algo": "vdbscan_leiden",
      "k_neighbors": 12,
      "eps_k_neighbors": 8,
      "leiden_resolution": 1.0,
      "cluster_min_samples": 3,
      "eps_estimation_based_on": "background",
      "vdbscan_sym_rule": "asymmetric"
    },
    "GLCTLVAML": {
      "cluster_algo": "vdbscan_leiden",
      "k_neighbors": 8,
      "eps_k_neighbors": 8,
      "leiden_resolution": 0.5,
      "cluster_min_samples": 3,
      "eps_estimation_based_on": "sample",
      "vdbscan_sym_rule": "asymmetric"
    }
  }
}
```

## Notes

- TCRNET outputs now go to `results/tcrnet/` and keep the standard names `cluster_members.txt` and `motif_pwms.txt`.
- TCRNET HTML visualizations are collected in `results/tcrnet/viz/`.
- TCRNET PDF figures go to the repository-level `figures/` directory.
- TCRNET respects `TCRNET_NPROC` for parallel R steps.
- REDCEA writes both chains into `results/redcea/`.
- REDCEA HTML visualizations are collected in `results/redcea/viz/`.
- REDCEA `cluster_members_TRA.txt` and `cluster_members_TRB.txt` are written directly into `results/redcea/`.
- REDCEA also writes `<chain>_vdjdb_clonotype_coords_2d.tsv` with exported 2D sample clonotype coordinates and `<chain>_background_coords_2d.tsv` with the shared background layout for downstream custom visualization.
- REDCEA uses `vdjdb_release/vdjdb.slim.txt` as the default VDJdb input table.
- REDCEA respects `REDCEA_NPROC`; when `REDCEA_CHAIN=both`, `TRA` and `TRB` are still launched in parallel as separate jobs.

## Possig Heatmaps

For large-epitope `vdbscan_leiden` reruns that already contain per-run
`run_metadata.tsv`, `*_summary_tcrempnet.tsv`, and
`knn_sample_sample__*.distances.npy`, use:

```bash
python scripts/build_redcea_possig_heatmaps.py \
  --runs-root <path-to-redcea_vdbscan_leiden_large_epitopes_rerun> \
  --output-dir results/redcea_possig_heatmaps
```

This script computes the repo-native `redcea_possig_density_score` components
and by default derives `d_ref` from the current rerun itself via
`median(d_epi)` across all discovered runs. If you need the original analytic
mode from the reporting scripts, pass `--d-ref-mode glc_ylq --glc-knn ... --ylq-knn ...`.

It writes:

- `run_level_metric_breakdown.tsv`
- `redcea_possig_parameter_runs.tsv`
- `run_level_distances.tsv`
- `epitope_distance_summary.tsv`
- `epitope_tightness_summary.tsv`
- `sample_clonotype_nn.tsv`
- `possig_cluster_lfc.tsv`
- `epitope_object_metric_summary.tsv`
- `parameter_robustness_summary.tsv`
- `best_parameter_by_epitope.tsv`
- `parameter_tightness_dependence.tsv`
- `selected_parameter_contrasts.tsv`
- `metric_collection_coverage.tsv`
- publication-style heatmaps in `png` / `pdf` / `svg`
- extra heatmaps with epitopes sorted by nearest-neighbor tightness and a reduced top-config panel chosen to contrast tight vs loose epitopes
- per-epitope diagnostic panels under `epitope_metric_panels/`, where:
  - the left subplot is the sample-clonotype Euclidean nearest-neighbor distance distribution
  - the right subplot is the `log_fold_change` distribution over all positive-significant enriched clusters
