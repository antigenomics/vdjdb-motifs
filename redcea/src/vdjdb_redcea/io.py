from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from .compat import BackgroundTransform
from .config import CHAIN_COLS, OutputPaths


def sanitize_filename_token(value: str) -> str:
    """Sanitize a string to be safe for use in filenames.

    Args:
        value: The string to sanitize.

    Returns:
        A sanitized string safe for filenames.
    """
    safe = "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in str(value).strip())
    return safe.strip("._") or "unknown"


def build_airr_from_epitope(ep_df: pd.DataFrame, chain: str) -> pd.DataFrame:
    """Build AIRR format table from epitope DataFrame for a given chain.

    Args:
        ep_df: DataFrame containing epitope data.
        chain: The chain type ('TRA' or 'TRB').

    Returns:
        DataFrame in AIRR format.
    """
    cfg = CHAIN_COLS[chain]
    out = ep_df[[cfg["cdr3"], cfg["v"], cfg["j"]]].copy()
    out.columns = ["junction_aa", "v_call", "j_call"]
    out["locus"] = cfg["locus"]
    logging.info("Built AIRR table with %d rows for chain %s", len(out), chain)
    out.dropna(subset=["junction_aa"], inplace=True)
    logging.info("After dropping rows with missing CDR3, %d rows remain", len(out))
    return out.reset_index(drop=True)


def prepare_output_dirs(output_root: Path) -> OutputPaths:
    """Prepare output directories and return OutputPaths object.

    Args:
        output_root: Root directory for outputs.

    Returns:
        OutputPaths object with all necessary directories.
    """
    paths = OutputPaths(
        output_root=output_root,
        viz_dir=output_root / "viz",
        airr_dir=output_root / "airr_format",
        tcremp_dir=output_root / "tcremp",
        tcrempnet_dir=output_root / "tcrempnet",
    )
    for path in (paths.viz_dir, paths.airr_dir, paths.tcremp_dir, paths.tcrempnet_dir):
        path.mkdir(parents=True, exist_ok=True)
    return paths


def get_background_transform_path(*, args, output_root: Path) -> Path:
    """Get the path for the background transform file.

    Args:
        args: Arguments object containing background_transform.
        output_root: Root output directory.

    Returns:
        Path to the background transform file.
    """
    default_path = output_root / "tcremp" / f"{args.chain.lower()}_background_transform.joblib"
    return Path(getattr(args, "background_transform", None) or default_path)


def get_background_umap_cache_path(*, transform_path: Path, n_bg_points: int) -> Path:
    """Get the cache path for background UMAP.

    Args:
        transform_path: Path to the transform file.
        n_bg_points: Number of background points.

    Returns:
        Path to the cached UMAP file.
    """
    return transform_path.with_name(f"{transform_path.stem}_bg_umap_{int(n_bg_points)}.npy")


def load_or_fit_background_transform(*, args, output_root: Path, bg_emb: pd.DataFrame) -> tuple[BackgroundTransform, Path]:
    """Load existing background transform or fit a new one.

    Args:
        args: Arguments object.
        output_root: Root output directory.
        bg_emb: Background embeddings DataFrame.

    Returns:
        Tuple of BackgroundTransform and its path.
    """
    transform_path = get_background_transform_path(args=args, output_root=output_root)

    def _transform_matches_requested_umap(transform: BackgroundTransform) -> bool:
        return (
            int(getattr(transform, "umap_n_neighbors", 15)) == int(args.umap_n_neighbors)
            and float(getattr(transform, "umap_min_dist", 0.1)) == float(args.umap_min_dist)
        )

    def _fit_and_save_transform() -> BackgroundTransform:
        logging.info(
            "Preparing to fit background transform: chain=%s, n_bg=%d, n_features=%s, path=%s",
            args.chain,
            len(bg_emb),
            getattr(bg_emb, "shape", None),
            transform_path,
        )
        logging.info("Instantiating BackgroundTransform")
        transform = BackgroundTransform(
            chain=args.chain,
            n_pca_components=args.cluster_pc_components,
            umap_n_neighbors=args.umap_n_neighbors,
            umap_min_dist=args.umap_min_dist,
            random_state=args.random_seed,
        )
        logging.info("Calling BackgroundTransform.fit(...)")
        transform.fit(bg_emb)
        logging.info(
            "BackgroundTransform.fit(...) finished: background_pca_shape=%s",
            getattr(getattr(transform, "background_pca_", None), "shape", None),
        )
        logging.info("Saving background transform to %s", transform_path)
        transform.save(transform_path)
        logging.info("Saved background transform to %s", transform_path)
        return transform

    if transform_path.exists():
        logging.info("Loading background transform from %s", transform_path)
        try:
            transform = BackgroundTransform.load(transform_path)
            logging.info(
                "Loaded cached background transform successfully: background_pca_shape=%s",
                getattr(getattr(transform, "background_pca_", None), "shape", None),
            )
            if not _transform_matches_requested_umap(transform):
                logging.info(
                    "Cached background transform UMAP parameters do not match requested values "
                    "(cached n_neighbors=%s, min_dist=%s; requested n_neighbors=%s, min_dist=%s). Recomputing.",
                    getattr(transform, "umap_n_neighbors", None),
                    getattr(transform, "umap_min_dist", None),
                    args.umap_n_neighbors,
                    args.umap_min_dist,
                )
                transform = _fit_and_save_transform()
        except ModuleNotFoundError as exc:
            logging.warning(
                "Failed to load cached background transform from %s due to legacy module path (%s). Recomputing it.",
                transform_path,
                exc,
            )
            transform = _fit_and_save_transform()
    else:
        transform = _fit_and_save_transform()

    return transform, transform_path


def fit_joint_umap(
    transform: BackgroundTransform,
    bg_pca: np.ndarray,
    sample_pca_blocks: list[np.ndarray],
    n_bg_points: int,
    *,
    transform_path: Path,
) -> tuple[np.ndarray, list[np.ndarray]]:
    """Fit UMAP on background PCA plus all epitope PCA blocks for plotting.

    PCA is still assumed to be fitted on background only. This function fits the
    2D plotting layout on the joint PCA representation and returns the
    corresponding coordinates split back into background and per-epitope blocks.
    """
    bg_pca_subset = np.asarray(bg_pca[:n_bg_points], dtype=np.float32)
    if bg_pca_subset.size == 0:
        raise ValueError("Cannot fit joint UMAP for empty background PCA array")

    sample_blocks = [np.asarray(block, dtype=np.float32) for block in sample_pca_blocks]
    joint_pca = np.vstack([bg_pca_subset] + sample_blocks).astype(np.float32, copy=False)

    logging.info(
        "Fitting plotting UMAP on %d background + %d sample clonotypes",
        len(bg_pca_subset),
        int(sum(len(block) for block in sample_blocks)),
    )
    joint_umap = np.asarray(transform.fit_umap(joint_pca), dtype=np.float32)
    transform.save(transform_path)
    logging.info("Saved background transform with fitted joint UMAP to %s", transform_path)

    bg_end = len(bg_pca_subset)
    bg_umap = joint_umap[:bg_end]

    sample_umap_blocks: list[np.ndarray] = []
    cursor = bg_end
    for block in sample_blocks:
        block_len = len(block)
        sample_umap_blocks.append(joint_umap[cursor : cursor + block_len])
        cursor += block_len

    cache_path = get_background_umap_cache_path(transform_path=transform_path, n_bg_points=len(bg_umap))
    np.save(cache_path, np.asarray(bg_umap, dtype=np.float32))
    logging.info("Saved cached background UMAP to %s", cache_path)

    return bg_umap, sample_umap_blocks


def build_sample_members_table(
    sample_cluster_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    chain: str,
    epitope_df: pd.DataFrame,
    epitope: str,
    sample_ids: pd.Series,
    sample_umap: np.ndarray,
) -> pd.DataFrame:
    """Build table of sample cluster members with metadata.

    Args:
        sample_cluster_df: DataFrame of sample clusters.
        summary_df: Summary DataFrame.
        chain: Chain type.
        epitope_df: Epitope DataFrame.
        epitope: Epitope name.
        sample_ids: Sample IDs.
        sample_umap: Sample UMAP coordinates.

    Returns:
        DataFrame of cluster members.
    """
    cfg = CHAIN_COLS[chain]
    cdr3_col = f"cdr3aa_{cfg['gene']}"
    v_col = f"v_{cfg['gene']}"
    j_col = f"j_{cfg['gene']}"

    meta = epitope_df.iloc[0]
    cluster_sizes = summary_df.set_index("cluster_id")["cluster_size"].to_dict()
    chain_tag = "A" if chain == "TRA" else "B"
    sample_cluster_df = sample_cluster_df.copy()
    umap_df = pd.DataFrame(
        {
            "clone_id": sample_ids.to_numpy(),
            "x": sample_umap[:, 0],
            "y": sample_umap[:, 1],
        }
    )
    sample_cluster_df = sample_cluster_df.merge(umap_df, on="clone_id", how="left")

    return pd.DataFrame(
        {
            "species": meta["species"],
            "antigen.epitope": epitope,
            "antigen.gene": meta["antigen.gene"],
            "antigen.species": meta["antigen.species"],
            "mhc.a": meta["mhc.a"],
            "mhc.b": meta["mhc.b"],
            "mhc.class": meta["mhc.class"],
            "gene": chain,
            "cdr3aa": sample_cluster_df[cdr3_col].values,
            "x": sample_cluster_df["x"].values,
            "y": sample_cluster_df["y"].values,
            "cid": [f"H.{chain_tag}.{epitope}.{int(cid)}" for cid in sample_cluster_df["cluster_id"].values],
            "csz": [int(cluster_sizes[int(cid)]) for cid in sample_cluster_df["cluster_id"].values],
            "v.segm": sample_cluster_df[v_col].values,
            "j.segm": sample_cluster_df[j_col].values,
            "v.end": pd.NA,
            "j.start": pd.NA,
            "v.segm.repr": sample_cluster_df[v_col].values,
            "j.segm.repr": sample_cluster_df[j_col].values,
        }
    )
