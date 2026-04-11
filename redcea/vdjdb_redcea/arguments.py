from __future__ import annotations

import argparse

from .config import CHAIN_COLS


def get_arguments_vdjdb_clusters(argv: list[str] | None = None):
    """Parse CLI arguments for the VDJdb REDCEA pipeline."""
    parser = argparse.ArgumentParser(
        prog="vdjdb-redcea",
        description="Cluster VDJdb clonotypes per epitope with REDCEA/TCRemP-style background correction.",
    )
    parser.add_argument("--vdjdb", required=True, help="Path to the VDJdb slim table.")
    parser.add_argument(
        "--background-airr",
        dest="background_airr",
        required=True,
        help="Path to the AIRR-format background repertoire table.",
    )
    parser.add_argument(
        "--background-embedding",
        dest="background_embedding",
        required=True,
        help="Path to precomputed background embeddings parquet.",
    )
    parser.add_argument("--output", required=True, help="Output directory root.")
    parser.add_argument("--chain", required=True, choices=sorted(CHAIN_COLS), help="TCR chain to analyze.")
    parser.add_argument("--species", default="HomoSapiens", help="Species name understood by SegmentLibrary.")
    parser.add_argument(
        "--epitopes",
        nargs="+",
        default=None,
        help="Optional list of epitopes to restrict analysis to.",
    )
    parser.add_argument(
        "--min-epitope-clonotypes",
        type=int,
        default=None,
        help="Skip epitopes with fewer clonotypes than this threshold.",
    )
    parser.add_argument(
        "--prototypes-path",
        default=None,
        help="Optional path to a prototype repertoire file. Defaults to TCRemP built-ins.",
    )
    parser.add_argument(
        "--index-col",
        default="id",
        help="Prototype repertoire column used as clonotype identifier.",
    )
    parser.add_argument(
        "--n-prototypes",
        type=int,
        default=None,
        help="Optional number of prototype clonotypes to retain.",
    )
    parser.add_argument(
        "--sample-random-clonotypes",
        action="store_true",
        help="Randomly sample prototypes when --n-prototypes is set.",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=42,
        help="Random seed for background transform and optional prototype sampling.",
    )
    parser.add_argument(
        "--umap-n-neighbors",
        type=int,
        default=15,
        help="UMAP neighborhood size for the joint plotting layout. Larger values usually make the map more global and compact.",
    )
    parser.add_argument(
        "--umap-min-dist",
        type=float,
        default=0.1,
        help="UMAP minimum distance for the joint plotting layout. Larger values usually reduce tight isolated islands.",
    )
    parser.add_argument("--nproc", type=int, default=1, help="Number of worker threads/processes.")
    parser.set_defaults(
        prefix=None,
        metrics="dissimilarity",
        sample=None,
        sample_embedding=None,
        background_embedding=None,
        n_clonotypes=None,
        lower_len_cdr3=None,
        higher_len_cdr3=None,
        cluster_algo="leiden",
        k_neighbors=15,
        cluster_min_samples=5,
        eps_k_neighbors=4,
        leiden_resolution=1.0,
        leiden_sub_resolution=1.0,
        eps_estimation_based_on="sample",
        vdbscan_sym_rule="asymmetric",
        cluster_pc_components=50,
        n_bg_points=None,
    )
    return parser.parse_args(argv)
