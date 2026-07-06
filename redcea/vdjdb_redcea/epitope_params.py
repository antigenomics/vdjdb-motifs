from __future__ import annotations

import json
from pathlib import Path

from .config import EpitopeClusteringParams


ALLOWED_EPITOPE_PARAM_KEYS = {
    "cluster_algo",
    "k_neighbors",
    "eps_k_neighbors",
    "leiden_resolution",
    "cluster_min_samples",
    "eps_estimation_based_on",
    "vdbscan_sym_rule",
    "leiden_sub_resolution",
}

REQUIRED_EPITOPE_PARAM_KEYS = ALLOWED_EPITOPE_PARAM_KEYS - {"leiden_sub_resolution"}


def load_epitope_clustering_overrides(path: str | Path) -> dict[str, EpitopeClusteringParams]:
    config_path = Path(path).resolve()
    payload = json.loads(config_path.read_text(encoding="utf-8"))

    if not isinstance(payload, dict):
        raise ValueError(f"Expected top-level object in epitope config: {config_path}")

    if "epitopes" in payload:
        payload = payload["epitopes"]

    if not isinstance(payload, dict):
        raise ValueError(f"Expected 'epitopes' mapping in epitope config: {config_path}")

    overrides: dict[str, EpitopeClusteringParams] = {}
    for epitope, raw_params in payload.items():
        if not isinstance(epitope, str) or not epitope.strip():
            raise ValueError(f"Epitope names must be non-empty strings in {config_path}")
        if not isinstance(raw_params, dict):
            raise ValueError(f"Expected parameter object for epitope {epitope!r} in {config_path}")

        unknown_keys = sorted(set(raw_params) - ALLOWED_EPITOPE_PARAM_KEYS)
        if unknown_keys:
            raise ValueError(
                f"Unsupported per-epitope parameter(s) for {epitope!r} in {config_path}: {', '.join(unknown_keys)}"
            )

        missing_keys = sorted(REQUIRED_EPITOPE_PARAM_KEYS - set(raw_params))
        if missing_keys:
            raise ValueError(
                f"Missing per-epitope parameter(s) for {epitope!r} in {config_path}: {', '.join(missing_keys)}"
            )

        overrides[epitope] = EpitopeClusteringParams(
            cluster_algo=str(raw_params["cluster_algo"]),
            k_neighbors=int(raw_params["k_neighbors"]),
            eps_k_neighbors=int(raw_params["eps_k_neighbors"]),
            leiden_resolution=float(raw_params["leiden_resolution"]),
            cluster_min_samples=int(raw_params["cluster_min_samples"]),
            eps_estimation_based_on=str(raw_params["eps_estimation_based_on"]),
            vdbscan_sym_rule=str(raw_params["vdbscan_sym_rule"]),
            leiden_sub_resolution=float(raw_params.get("leiden_sub_resolution", 1.0)),
        )

    return overrides
