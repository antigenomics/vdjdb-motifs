from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


CHAIN_COLS = {
    "TRA": {"cdr3": "cdr3", "v": "v.segm", "j": "j.segm", "locus": "alpha", "gene": "alpha"},
    "TRB": {"cdr3": "cdr3", "v": "v.segm", "j": "j.segm", "locus": "beta", "gene": "beta"},
}

DEFAULT_PLOT_BG_POINTS = 100_000


@dataclass(frozen=True)
class EpitopeClusteringParams:
    cluster_algo: str
    k_neighbors: int
    eps_k_neighbors: int
    leiden_resolution: float
    cluster_min_samples: int
    eps_estimation_based_on: str
    vdbscan_sym_rule: str
    leiden_sub_resolution: float

    @classmethod
    def from_args(cls, args) -> "EpitopeClusteringParams":
        return cls(
            cluster_algo=str(args.cluster_algo),
            k_neighbors=int(args.k_neighbors),
            eps_k_neighbors=int(args.eps_k_neighbors),
            leiden_resolution=float(args.leiden_resolution),
            cluster_min_samples=int(args.cluster_min_samples),
            eps_estimation_based_on=str(args.eps_estimation_based_on),
            vdbscan_sym_rule=str(args.vdbscan_sym_rule),
            leiden_sub_resolution=float(args.leiden_sub_resolution),
        )

    def with_overrides(
        self,
        *,
        cluster_algo: str | None = None,
        k_neighbors: int | None = None,
        eps_k_neighbors: int | None = None,
        leiden_resolution: float | None = None,
        cluster_min_samples: int | None = None,
        eps_estimation_based_on: str | None = None,
        vdbscan_sym_rule: str | None = None,
        leiden_sub_resolution: float | None = None,
    ) -> "EpitopeClusteringParams":
        return EpitopeClusteringParams(
            cluster_algo=self.cluster_algo if cluster_algo is None else str(cluster_algo),
            k_neighbors=self.k_neighbors if k_neighbors is None else int(k_neighbors),
            eps_k_neighbors=self.eps_k_neighbors if eps_k_neighbors is None else int(eps_k_neighbors),
            leiden_resolution=self.leiden_resolution if leiden_resolution is None else float(leiden_resolution),
            cluster_min_samples=self.cluster_min_samples if cluster_min_samples is None else int(cluster_min_samples),
            eps_estimation_based_on=(
                self.eps_estimation_based_on
                if eps_estimation_based_on is None
                else str(eps_estimation_based_on)
            ),
            vdbscan_sym_rule=self.vdbscan_sym_rule if vdbscan_sym_rule is None else str(vdbscan_sym_rule),
            leiden_sub_resolution=(
                self.leiden_sub_resolution if leiden_sub_resolution is None else float(leiden_sub_resolution)
            ),
        )


@dataclass(frozen=True)
class OutputPaths:
    output_root: Path
    viz_dir: Path
    airr_dir: Path
    airr_processed_dir: Path
    tcremp_dir: Path
    tcrempnet_dir: Path
