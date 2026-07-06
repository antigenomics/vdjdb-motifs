#!/usr/bin/env python
from __future__ import annotations

import gc
import logging
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .arguments import get_arguments_vdjdb_clusters
from .config import CHAIN_COLS, DEFAULT_PLOT_BG_POINTS, EpitopeClusteringParams


@dataclass(frozen=True)
class EpitopeClusteringArtifacts:
    cluster_df: pd.DataFrame
    summary_df: pd.DataFrame
    significant_cluster_ids: set[int]
    sample_cluster_df: pd.DataFrame
    enriched_sample_cluster_df: pd.DataFrame


def _normalize_vdjdb_columns(vdjdb_df: pd.DataFrame, chain: str) -> pd.DataFrame:
    """Normalize VDJdb inputs to generic per-chain columns.

    After selecting the requested chain, the rest of the pipeline works only
    with ``cdr3``/``v.segm``/``j.segm`` regardless of whether the source table
    was a slim export or a full VDJdb dump.
    """
    normalized = vdjdb_df.copy()

    if "gene" in normalized.columns:
        before_gene_filter = len(normalized)
        normalized = normalized[normalized["gene"].astype(str).str.upper() == chain].copy()
        logging.info(
            "Filtered generic VDJdb table by gene=%s: %d -> %d",
            chain,
            before_gene_filter,
            len(normalized),
        )

    generic_targets = {
        "cdr3": {"TRA": "cdr3.alpha", "TRB": "cdr3.beta"}[chain],
        "v.segm": {"TRA": "v.alpha", "TRB": "v.beta"}[chain],
        "j.segm": {"TRA": "j.alpha", "TRB": "j.beta"}[chain],
    }
    renamed_cols = {}
    for target_col, source_col in generic_targets.items():
        if target_col not in normalized.columns and source_col in normalized.columns:
            renamed_cols[source_col] = target_col

    if renamed_cols:
        normalized = normalized.rename(columns=renamed_cols)
        logging.info("Renamed VDJdb columns to generic names for %s chain: %s", chain, renamed_cols)

    return normalized


def _cluster_members_filename(chain: str, output_tag: str | None, *, prefix: str = "cluster_members") -> str:
    if output_tag is None:
        return f"{prefix}_{chain}.txt"
    clean_tag = "".join(char if char.isalnum() or char in {"-", "_"} else "_" for char in output_tag).strip("_")
    if not clean_tag:
        return f"{prefix}_{chain}.txt"
    return f"{prefix}_{chain}_{clean_tag}.txt"


def resolve_joint_knn(
    sample_pca: np.ndarray,
    bg_pca: np.ndarray,
    k_neighbors: int,
    nproc: int,
    sample_index_path: Path,
    bg_index_path: Path,
):
    """Resolve joint kNN for sample and background PCA.

    Args:
        sample_pca: Sample PCA array.
        bg_pca: Background PCA array.
        k_neighbors: Number of neighbors.
        nproc: Number of processes.
        sample_index_path: Path to sample index.
        bg_index_path: Path to background index.

    Returns:
        KNN result.
    """
    from .compat import build_joint_knn_from_split, compute_blockwise_knn_merged

    if sample_index_path is not None and bg_index_path is not None:
        try:
            result = compute_blockwise_knn_merged(
                bg=bg_pca,
                sample=sample_pca,
                k_neighbors=k_neighbors,
                bg_index_path=bg_index_path,
                sample_index_path=sample_index_path,
                rebuild_bg=False,
                rebuild_sample=False,
                save_blocks=False,
                output_dir=None,
                nproc=nproc,
            )
            if isinstance(result, tuple) and len(result) == 2:
                return result
            if isinstance(result, tuple) and len(result) == 8:
                dist_ss, ind_ss, dist_bb, ind_bb, dist_sb, ind_sb, dist_bs, ind_bs = result
                return build_joint_knn_from_split(
                    dist_ss=dist_ss,
                    ind_ss=ind_ss,
                    dist_bb=dist_bb,
                    ind_bb=ind_bb,
                    dist_sb=dist_sb,
                    ind_sb=ind_sb,
                    dist_bs=dist_bs,
                    ind_bs=ind_bs,
                    k_out=k_neighbors,
                )
            logging.warning(
                "Unexpected blockwise kNN result for sample index %s and background index %s; "
                "falling back to direct joint FAISS search.",
                sample_index_path,
                bg_index_path,
            )
        except Exception as exc:
            logging.warning(
                "Failed to load cached FAISS indexes (%s, %s): %s. "
                "Falling back to direct joint FAISS search.",
                sample_index_path,
                bg_index_path,
                exc,
            )
    else:
        logging.warning("Missing cached FAISS indexes; falling back to direct joint FAISS search.")

    import faiss

    joint = np.vstack([sample_pca, bg_pca]).astype("float32", copy=False)
    index = faiss.IndexFlatL2(joint.shape[1])
    faiss.omp_set_num_threads(int(nproc))
    index.add(joint)
    return index.search(joint, k_neighbors)


def _compute_sample_embeddings(
    args, genes: list[str], locus: str, lib, proto, paths, chain: str, prefix: str, airr_path: Path
):
    """Compute or load sample embeddings."""
    from .compat import compute_embeddings_if_needed, load_embedding_artifacts, normalize_config

    args.sample = str(airr_path.resolve())
    args.sample_embedding = str((paths.tcremp_dir / f"{prefix}_sample_embeddings.parquet").resolve())
    args.prefix = prefix
    config = normalize_config(args)

    compute_embeddings_if_needed(
        path=config.sample,
        config=config,
        is_sample=True,
        proto=proto,
        chain=genes,
        lib=lib,
        locus=locus,
        prefix=prefix,
        output_path=paths.tcremp_dir,
    )

    return load_embedding_artifacts(
        path=config.sample,
        args=config,
        is_sample=True,
        lib=lib,
        locus=locus,
        prefix=prefix,
        output_path=paths.tcremp_dir,
    )


def _load_precomputed_sample_embeddings(args, lib, locus: str, paths, prefix: str, airr_path: Path):
    """Load previously computed per-epitope embeddings from saved files."""
    from .compat import load_embedding_artifacts, normalize_config

    args.sample = str(airr_path.resolve())
    args.sample_embedding = str((paths.tcremp_dir / f"{prefix}_sample_embeddings.parquet").resolve())
    args.prefix = prefix
    config = normalize_config(args)

    return load_embedding_artifacts(
        path=config.sample,
        args=config,
        is_sample=True,
        lib=lib,
        locus=locus,
        prefix=prefix,
        output_path=paths.tcremp_dir,
    )


def _perform_clustering(
    sample_pca: np.ndarray,
    bg_pca: np.ndarray,
    clustering_params: EpitopeClusteringParams,
    nproc: int,
    sample_index_path: Path,
    bg_index_path: Path,
    sample_ids: pd.Series,
    bg_ids: pd.Series,
):
    """Perform clustering on sample and background."""
    from .compat import run_leiden_clustering

    distances, indices = resolve_joint_knn(
        sample_pca,
        bg_pca,
        clustering_params.k_neighbors,
        nproc,
        sample_index_path=sample_index_path,
        bg_index_path=bg_index_path,
    )
    labels = run_leiden_clustering(
        knn_indices=indices,
        knn_distances=distances,
        resolution=clustering_params.leiden_resolution,
        n_threads=nproc,
        min_cluster_size=clustering_params.cluster_min_samples,
        min_cluster_size_mask=np.arange(len(sample_pca) + len(bg_pca)) < len(sample_pca),
    )
    return sample_pca, labels


def _resolve_epitope_clustering_params(args, epitope: str) -> EpitopeClusteringParams:
    params = EpitopeClusteringParams.from_args(args)
    overrides = getattr(args, "epitope_clustering_overrides", {}) or {}
    override = overrides.get(epitope)
    if override is None:
        return params
    return params.with_overrides(
        cluster_algo=override.cluster_algo,
        k_neighbors=override.k_neighbors,
        eps_k_neighbors=override.eps_k_neighbors,
        leiden_resolution=override.leiden_resolution,
        cluster_min_samples=override.cluster_min_samples,
        eps_estimation_based_on=override.eps_estimation_based_on,
        vdbscan_sym_rule=override.vdbscan_sym_rule,
        leiden_sub_resolution=override.leiden_sub_resolution,
    )


def _build_clonotype_coords_table(
    *,
    epitope: str,
    chain: str,
    sample_reps: pd.DataFrame,
    sample_ids: pd.Series,
    sample_labels: np.ndarray,
    summary_df: pd.DataFrame,
    sample_umap,
    clustering_params: EpitopeClusteringParams,
) -> pd.DataFrame:
    cfg = CHAIN_COLS[chain]
    cdr3_col = f"cdr3aa_{cfg['gene']}"
    v_col = f"v_{cfg['gene']}"
    j_col = f"j_{cfg['gene']}"

    df = sample_reps.copy().reset_index(drop=True)
    df["clone_id"] = sample_ids.to_numpy()
    df["cluster_id"] = np.asarray(sample_labels[: len(sample_ids)], dtype=np.int32)
    if sample_umap is None:
        df["x"] = pd.NA
        df["y"] = pd.NA
    else:
        df["x"] = sample_umap[:, 0]
        df["y"] = sample_umap[:, 1]

    summary_by_cluster = summary_df.set_index("cluster_id")
    df["significant"] = df["cluster_id"].map(summary_by_cluster["significant"]).fillna(False)
    df["cluster_size_sample"] = df["cluster_id"].map(summary_by_cluster["sample"])
    df["log_fold_change"] = df["cluster_id"].map(summary_by_cluster["log_fold_change"])

    return pd.DataFrame(
        {
            "chain": chain,
            "epitope": epitope,
            "clone_id": df["clone_id"],
            "cdr3aa": df[cdr3_col] if cdr3_col in df.columns else pd.NA,
            "v.segm": df[v_col] if v_col in df.columns else pd.NA,
            "j.segm": df[j_col] if j_col in df.columns else pd.NA,
            "x": df["x"],
            "y": df["y"],
            "cluster_id": df["cluster_id"],
            "significant": df["significant"],
            "cluster_size_sample": df["cluster_size_sample"],
            "log_fold_change": df["log_fold_change"],
            "cluster_algo": clustering_params.cluster_algo,
            "k_neighbors": clustering_params.k_neighbors,
            "eps_k_neighbors": clustering_params.eps_k_neighbors,
            "leiden_resolution": clustering_params.leiden_resolution,
            "cluster_min_samples": clustering_params.cluster_min_samples,
            "eps_estimation_based_on": clustering_params.eps_estimation_based_on,
            "vdbscan_sym_rule": clustering_params.vdbscan_sym_rule,
            "leiden_sub_resolution": clustering_params.leiden_sub_resolution,
        }
    )


def _build_background_coords_table(*, chain: str, bg_umap: np.ndarray) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chain": chain,
            "background_index": np.arange(len(bg_umap), dtype=np.int32),
            "x": bg_umap[:, 0],
            "y": bg_umap[:, 1],
        }
    )


def _format_epitope_param_tag(params: EpitopeClusteringParams) -> str:
    resolution_tag = str(params.leiden_resolution).replace(".", "p")
    eps_mode = "".join(ch for ch in params.eps_estimation_based_on if ch.isalnum() or ch in {"-", "_"})
    sym_rule = "".join(ch for ch in params.vdbscan_sym_rule if ch.isalnum() or ch in {"-", "_"})
    algo = "".join(ch for ch in params.cluster_algo if ch.isalnum() or ch in {"-", "_"})
    return (
        f"{algo}_k{params.k_neighbors}_ek{params.eps_k_neighbors}_"
        f"res{resolution_tag}_min{params.cluster_min_samples}_eps{eps_mode}_sym{sym_rule}"
    )


def _build_sample_labels(sample_ids: pd.Series, sample_cluster_df: pd.DataFrame) -> np.ndarray:
    cluster_lookup = (
        sample_cluster_df.loc[:, ["clone_id", "cluster_id"]]
        .drop_duplicates(subset=["clone_id"])
        .set_index("clone_id")["cluster_id"]
        .to_dict()
    )
    return np.asarray([int(cluster_lookup.get(clone_id, -1)) for clone_id in sample_ids], dtype=np.int32)


def _run_external_vdbscan_clustering(
    *,
    epitope_data: dict[str, object],
    args,
    paths,
    clustering_params: EpitopeClusteringParams,
    sample_ids: pd.Series,
) -> EpitopeClusteringArtifacts:
    prefix = str(epitope_data["prefix"])
    epitope = str(epitope_data["epitope"])
    run_dir = paths.output_root / "per_epitope_cluster_runs" / prefix / _format_epitope_param_tag(clustering_params)
    run_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable,
        "-m",
        "redcea.redcea",
        "--sample",
        str(Path(epitope_data["airr_path"]).resolve()),
        "--background",
        str(Path(args.background).resolve()),
        "--output",
        str(run_dir),
        "--prefix",
        prefix,
        "--chain",
        str(args.chain),
        "--species",
        str(args.species),
        "--metrics",
        str(args.metrics),
        "--sample-embedding",
        str(Path(epitope_data["sample_embedding_path"]).resolve()),
        "--background-embedding",
        str(Path(args.background_embedding).resolve()),
        "--cluster-pc-components",
        str(args.cluster_pc_components),
        "--cluster-algo",
        str(clustering_params.cluster_algo),
        "--core-min-samples",
        str(clustering_params.cluster_min_samples),
        "--k-neighbors",
        str(clustering_params.k_neighbors),
        "--eps-k-neighbors",
        str(clustering_params.eps_k_neighbors),
        "--leiden-resolution",
        str(clustering_params.leiden_resolution),
        "--leiden-sub-resolution",
        str(clustering_params.leiden_sub_resolution),
        "--eps-estimation-based-on",
        str(clustering_params.eps_estimation_based_on),
        "--vdbscan-sym-rule",
        str(clustering_params.vdbscan_sym_rule),
        "--random-seed",
        str(args.random_seed),
        "--nproc",
        str(args.nproc),
    ]
    if args.n_bg_points is not None:
        cmd.extend(["--n-bg-points", str(args.n_bg_points)])

    logging.info("Running external clustering backend for epitope %s: %s", epitope, " ".join(cmd))
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            "Could not run external vdbscan backend because `redcea.redcea` is unavailable in this environment."
        ) from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            "External vdbscan backend failed for "
            f"{epitope} with exit code {exc.returncode}.\nSTDOUT:\n{exc.stdout}\nSTDERR:\n{exc.stderr}"
        ) from exc

    cluster_path = run_dir / f"{prefix}_tcremp_clusters.tsv"
    summary_path = run_dir / f"{prefix}_summary_tcrempnet.tsv"
    if not cluster_path.exists():
        raise FileNotFoundError(f"Missing external cluster output: {cluster_path}")
    if not summary_path.exists():
        raise FileNotFoundError(f"Missing external summary output: {summary_path}")

    cluster_df = pd.read_csv(cluster_path, sep="\t")
    summary_df = pd.read_csv(summary_path, sep="\t")
    if "significant" not in summary_df.columns:
        summary_df["significant"] = (
            (pd.to_numeric(summary_df["enrichment_fdr_zbinom"], errors="coerce") < 0.05)
            & (pd.to_numeric(summary_df["log_fold_change"], errors="coerce") > 0)
        )

    sample_id_set = set(sample_ids.astype(str))
    cluster_df["clone_id"] = cluster_df["clone_id"].astype(str)
    summary_df["cluster_id"] = pd.to_numeric(summary_df["cluster_id"], errors="coerce").astype(int)
    significant_cluster_ids = set(summary_df.loc[summary_df["significant"], "cluster_id"].astype(int))
    sample_cluster_df = cluster_df[
        cluster_df["clone_id"].isin(sample_id_set) & pd.to_numeric(cluster_df["cluster_id"], errors="coerce").ne(-1)
    ].copy()
    sample_cluster_df["cluster_id"] = pd.to_numeric(sample_cluster_df["cluster_id"], errors="coerce").astype(int)
    enriched_sample_cluster_df = sample_cluster_df[
        sample_cluster_df["cluster_id"].isin(significant_cluster_ids)
    ].copy()

    return EpitopeClusteringArtifacts(
        cluster_df=cluster_df,
        summary_df=summary_df,
        significant_cluster_ids=significant_cluster_ids,
        sample_cluster_df=sample_cluster_df,
        enriched_sample_cluster_df=enriched_sample_cluster_df,
    )


def _precompute_epitope_embeddings(
    epitope: str,
    ep_df: pd.DataFrame,
    *,
    args,
    genes: list[str],
    locus: str,
    lib: SegmentLibrary,
    proto,
    paths,
) -> dict[str, object]:
    """Compute and persist per-epitope embeddings before analysis."""
    from .io import build_airr_from_epitope, build_processed_airr_from_tcremp_representations

    chain = args.chain
    prefix = f"{chain.lower()}_vdjdb_{epitope}"
    airr_path = paths.airr_dir / f"{prefix}.tsv"
    processed_airr_path = paths.airr_processed_dir / f"{prefix}.tsv"

    logging.info("Preparing epitope %s with %d clonotypes", epitope, len(ep_df))
    raw_airr_df = build_airr_from_epitope(ep_df, chain)
    raw_airr_df.to_csv(airr_path, sep="\t", index=False)
    sample_artifacts = _compute_sample_embeddings(
        args, genes, locus, lib, proto, paths, chain, prefix, airr_path
    )
    processed_airr_df = build_processed_airr_from_tcremp_representations(sample_artifacts.representations, chain)
    processed_airr_df.to_csv(processed_airr_path, sep="\t", index=False)
    logging.info(
        "Saved processed AIRR table for epitope %s: raw_rows=%d, tcremp_rows=%d, path=%s",
        epitope,
        len(raw_airr_df),
        len(processed_airr_df),
        processed_airr_path,
    )
    return {
        "epitope": epitope,
        "ep_df": ep_df,
        "prefix": prefix,
        "airr_path": airr_path,
        "processed_airr_path": processed_airr_path,
        "sample_embedding_path": paths.tcremp_dir / f"{prefix}_sample_embeddings.parquet",
        "n_sample_rows": len(sample_artifacts.ids),
    }


def _load_epitope_input(
    epitope_info: dict[str, object],
    *,
    args,
    lib: SegmentLibrary,
    locus: str,
    paths,
) -> dict[str, object]:
    """Load precomputed per-epitope embeddings and metadata for analysis."""
    epitope = epitope_info["epitope"]
    prefix = epitope_info["prefix"]
    airr_path = epitope_info["airr_path"]

    logging.info("Loading precomputed embeddings for epitope %s", epitope)
    sample_artifacts = _load_precomputed_sample_embeddings(
        args, lib, locus, paths, prefix, airr_path
    )
    loaded = dict(epitope_info)
    loaded.update(
        {
            "sample_emb": sample_artifacts.embeddings,
            "sample_reps": sample_artifacts.representations,
            "sample_ids": sample_artifacts.ids,
            "sample_index_path": sample_artifacts.cache_path,
        }
    )
    return loaded


def _stage_precompute_epitope_embeddings(
    vdjdb_df: pd.DataFrame,
    *,
    args,
    genes: list[str],
    locus: str,
    lib: SegmentLibrary,
    proto,
    paths,
) -> list[dict[str, object]]:
    """Stage 1: compute and persist per-epitope embeddings."""
    logging.info("Stage 1/4: precomputing TCRemP embeddings for selected epitopes")
    epitope_infos: list[dict[str, object]] = []
    failed_epitopes: list[dict[str, object]] = []
    chain_lower = args.chain.lower()
    for epitope, ep_df in vdjdb_df.groupby("antigen.epitope", sort=True):
        prefix = f"{chain_lower}_vdjdb_{epitope}"
        airr_path = paths.airr_dir / f"{prefix}.tsv"
        processed_airr_path = paths.airr_processed_dir / f"{prefix}.tsv"
        try:
            epitope_info = _precompute_epitope_embeddings(
                epitope,
                ep_df,
                args=args,
                genes=genes,
                locus=locus,
                lib=lib,
                proto=proto,
                paths=paths,
            )
            epitope_infos.append(epitope_info)
        except Exception as exc:
            logging.exception(
                "Skipping epitope %s for chain %s after a preprocessing/embedding failure: %s",
                epitope,
                args.chain,
                exc,
            )
            failed_epitopes.append(
                {
                    "epitope": epitope,
                    "chain": args.chain,
                    "n_input_rows": len(ep_df),
                    "stage": "precompute_embeddings",
                    "error_type": type(exc).__name__,
                    "error_message": str(exc),
                    "airr_path": str(airr_path),
                    "processed_airr_path": str(processed_airr_path),
                }
            )
            continue

    if failed_epitopes:
        failed_path = paths.output_root / f"{chain_lower}_skipped_epitopes.tsv"
        pd.DataFrame(failed_epitopes).to_csv(failed_path, sep="\t", index=False)
        logging.warning(
            "Skipped %d epitope(s) during Stage 1/4 for chain %s; wrote failure report to %s",
            len(failed_epitopes),
            args.chain,
            failed_path,
        )

    if not epitope_infos:
        raise RuntimeError(
            f"No epitopes could be prepared successfully for chain {args.chain}. "
            f"See {paths.output_root / f'{chain_lower}_skipped_epitopes.tsv'} for failure details."
            if failed_epitopes
            else f"No epitopes could be prepared successfully for chain {args.chain}."
        )

    logging.info("Stage 1/4 done: prepared %d epitope embedding files", len(epitope_infos))
    return epitope_infos


def _stage_fit_background_transform(*, args, output_root: Path, bg_emb: pd.DataFrame):
    """Stage 2: fit or load background-only PCA transform."""
    from .io import load_or_fit_background_transform

    logging.info(
        "Stage 2/4: fitting or loading background transform for chain=%s with bg_emb_shape=%s",
        args.chain,
        getattr(bg_emb, "shape", None),
    )
    transform, transform_path = load_or_fit_background_transform(args=args, output_root=output_root, bg_emb=bg_emb)
    bg_pca = transform.background_pca_
    if bg_pca is None:
        logging.info("Cached transform has no background_pca_; calling transform.transform_pca(bg_emb)")
        bg_pca = transform.transform_pca(bg_emb)
        logging.info("transform.transform_pca(bg_emb) finished: bg_pca_shape=%s", getattr(bg_pca, "shape", None))

    # The full background embedding matrix is no longer needed after PCA/transform fit.
    del bg_emb
    gc.collect()

    logging.info("Stage 2/4 done: background PCA shape=%s", getattr(bg_pca, "shape", None))
    return transform, transform_path, bg_pca


def _stage_prepare_joint_umap(
    epitope_infos: list[dict[str, object]],
    *,
    args,
    lib: SegmentLibrary,
    locus: str,
    paths,
    transform,
    bg_pca: np.ndarray,
    transform_path: Path,
) -> tuple[list[dict[str, object]], np.ndarray]:
    """Stage 3: load precomputed sample embeddings, map to background PCA, fit joint UMAP."""
    from .io import fit_joint_umap

    logging.info("Stage 3/4: loading precomputed epitope embeddings and fitting joint plotting UMAP")
    epitope_inputs: list[dict[str, object]] = []
    sample_pca_blocks: list[np.ndarray] = []
    for epitope_info in epitope_infos:
        epitope_data = _load_epitope_input(
            epitope_info,
            args=args,
            lib=lib,
            locus=locus,
            paths=paths,
        )
        sample_pca = transform.transform_pca(epitope_data["sample_emb"])
        epitope_data["sample_pca"] = sample_pca
        epitope_inputs.append(epitope_data)
        sample_pca_blocks.append(sample_pca)
        del epitope_data["sample_emb"]

    plot_bg_points = min(len(bg_pca), args.n_bg_points or DEFAULT_PLOT_BG_POINTS)
    bg_umap, sample_umap_blocks = fit_joint_umap(
        transform,
        bg_pca,
        sample_pca_blocks,
        plot_bg_points,
        transform_path=transform_path,
    )

    for epitope_data, sample_umap in zip(epitope_inputs, sample_umap_blocks):
        epitope_data["sample_umap"] = sample_umap

    # We keep per-epitope sample_pca for clustering, but the temporary list of all
    # PCA blocks is no longer needed once joint UMAP has been fitted.
    del sample_pca_blocks, sample_umap_blocks
    gc.collect()

    logging.info("Stage 3/4 done: joint UMAP prepared for %d epitopes", len(epitope_inputs))
    return epitope_inputs, bg_umap


def _stage_run_per_epitope_analysis(
    epitope_inputs: list[dict[str, object]],
    *,
    args,
    paths,
    bg_pca: np.ndarray,
    bg_reps: pd.DataFrame,
    bg_ids: pd.Series,
    bg_index_path: Path,
    bg_umap,
) -> tuple[list[pd.DataFrame], list[pd.DataFrame], list[pd.DataFrame], list[pd.DataFrame], list[dict[str, object]]]:
    """Stage 4: run per-epitope clustering and collect exported tables."""
    logging.info("Stage 4/4: running per-epitope clustering analysis")
    clustered_tables: list[pd.DataFrame] = []
    cluster_members_tables: list[pd.DataFrame] = []
    all_cluster_members_tables: list[pd.DataFrame] = []
    clonotype_coords_tables: list[pd.DataFrame] = []
    parameter_rows: list[dict[str, object]] = []
    for epitope_data in epitope_inputs:
        clustered_df, cluster_members_df, all_cluster_members_df, clonotype_coords_df, parameter_row = process_epitope(
            epitope_data,
            args=args,
            paths=paths,
            bg_pca=bg_pca,
            bg_reps=bg_reps,
            bg_ids=bg_ids,
            bg_index_path=bg_index_path,
            bg_umap=bg_umap,
        )
        clustered_tables.append(clustered_df)
        cluster_members_tables.append(cluster_members_df)
        all_cluster_members_tables.append(all_cluster_members_df)
        clonotype_coords_tables.append(clonotype_coords_df)
        parameter_rows.append(parameter_row)
    logging.info("Stage 4/4 done")
    return clustered_tables, cluster_members_tables, all_cluster_members_tables, clonotype_coords_tables, parameter_rows


def _compute_summary_and_significance(summary_df: pd.DataFrame, sample_ids: pd.Series, bg_ids: pd.Series):
    """Compute summary statistics and significance."""
    from .compat import add_log_fold_change, add_z_binom_pvalues

    summary_df = add_z_binom_pvalues(summary_df, total_sample=len(sample_ids), total_background=len(bg_ids))
    summary_df = add_log_fold_change(summary_df, total_sample=len(sample_ids), total_background=len(bg_ids))
    summary_df["significant"] = (
        (summary_df["enrichment_fdr_zbinom"] < 0.05) & (summary_df["log_fold_change"] > 0)
    )
    significant_cluster_ids = set(summary_df.loc[summary_df["significant"], "cluster_id"].astype(int))
    return summary_df, significant_cluster_ids


def _save_cluster_results(
    cluster_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    sample_cluster_df: pd.DataFrame,
    cluster_members_df: pd.DataFrame,
    paths,
    prefix: str,
):
    """Save clustering results to files."""
    cluster_df.to_csv(paths.tcrempnet_dir / f"{prefix}_tcremp_clusters.tsv", sep="\t", index=False)
    summary_df[
        [
            "cluster_id",
            "cluster_size",
            "sample",
            "background",
            "enrichment_pvalue_zbinom",
            "enrichment_fdr_zbinom",
            "log_fold_change",
            "significant",
        ]
    ].to_csv(paths.tcrempnet_dir / f"{prefix}_summary_tcrempnet.tsv", sep="\t", index=False)
    sample_cluster_df.to_csv(paths.tcrempnet_dir / f"{prefix}_clustered_sample_clonotypes.tsv", sep="\t", index=False)
    cluster_members_df.to_csv(paths.tcrempnet_dir / f"{prefix}_cluster_members.tsv", sep="\t", index=False)


def _build_epitope_clustering_artifacts(
    *,
    labels,
    sample_ids: pd.Series,
    bg_ids: pd.Series,
    sample_reps: pd.DataFrame,
    bg_reps: pd.DataFrame,
) -> EpitopeClusteringArtifacts:
    """Build the per-epitope clustering tables.

    This mirrors the shared RedCEA pattern of:
    1. creating a joint cluster table,
    2. computing summary statistics,
    3. extracting enriched/significant sample clusters,
    while preserving the exact vdjdb_redcea selection logic.
    """
    from .compat import compute_cluster_summary

    joint_ids = pd.concat([sample_ids, bg_ids], ignore_index=True)
    joint_reps = pd.concat([sample_reps, bg_reps], ignore_index=True)
    cluster_df = pd.DataFrame({"clone_id": joint_ids, "cluster_id": labels}).merge(
        joint_reps, on="clone_id", how="left"
    )

    summary_df = compute_cluster_summary(cluster_df.copy(), sample_ids)
    summary_df, significant_cluster_ids = _compute_summary_and_significance(summary_df, sample_ids, bg_ids)

    sample_cluster_df = cluster_df[
        (cluster_df["clone_id"].isin(set(sample_ids))) & (cluster_df["cluster_id"] != -1)
    ].copy()
    enriched_sample_cluster_df = sample_cluster_df[
        sample_cluster_df["cluster_id"].isin(significant_cluster_ids)
    ].copy()

    return EpitopeClusteringArtifacts(
        cluster_df=cluster_df,
        summary_df=summary_df,
        significant_cluster_ids=significant_cluster_ids,
        sample_cluster_df=sample_cluster_df,
        enriched_sample_cluster_df=enriched_sample_cluster_df,
    )


def process_epitope(
    epitope_data: dict[str, object],
    *,
    args,  # type: ignore
    paths,
    bg_pca: np.ndarray,
    bg_reps: pd.DataFrame,
    bg_ids: pd.Series,
    bg_index_path: Path,
    bg_umap,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    """Process a single epitope for clustering.

    Args:
        epitope_data: Prepared per-epitope inputs.
        args: Arguments object.
        paths: OutputPaths.
        bg_pca: Background PCA.
        bg_reps: Background representations.
        bg_ids: Background IDs.
        bg_index_path: Background index path.
        bg_umap: Background UMAP.

    Returns:
        Cluster outputs, 2D clonotype coordinates, and effective parameters.
    """
    from .io import build_sample_members_table, sanitize_filename_token
    from .plotting import save_cluster_plot_html

    epitope = epitope_data["epitope"]
    ep_df = epitope_data["ep_df"]
    prefix = epitope_data["prefix"]
    sample_reps = epitope_data["sample_reps"]
    sample_ids = epitope_data["sample_ids"]
    sample_index_path = epitope_data["sample_index_path"]
    sample_pca = epitope_data["sample_pca"]
    sample_umap = epitope_data.get("sample_umap")

    chain = args.chain
    clustering_params = _resolve_epitope_clustering_params(args, str(epitope))

    logging.info(
        "Processing epitope %s with %d clonotypes using algo=%s, k=%d, eps_k=%d, resolution=%s, min_samples=%d, eps_mode=%s, sym_rule=%s",
        epitope,
        len(ep_df),
        clustering_params.cluster_algo,
        clustering_params.k_neighbors,
        clustering_params.eps_k_neighbors,
        clustering_params.leiden_resolution,
        clustering_params.cluster_min_samples,
        clustering_params.eps_estimation_based_on,
        clustering_params.vdbscan_sym_rule,
    )
    if clustering_params.cluster_algo == "leiden":
        sample_pca, labels = _perform_clustering(
            sample_pca,
            bg_pca,
            clustering_params,
            args.nproc,
            sample_index_path,
            bg_index_path,
            sample_ids,
            bg_ids,
        )
        artifacts = _build_epitope_clustering_artifacts(
            labels=labels,
            sample_ids=sample_ids,
            bg_ids=bg_ids,
            sample_reps=sample_reps,
            bg_reps=bg_reps,
        )
        sample_labels = np.asarray(labels[: len(sample_ids)], dtype=np.int32)
    elif clustering_params.cluster_algo in {"vdbscan_leiden", "vdbscan"}:
        artifacts = _run_external_vdbscan_clustering(
            epitope_data=epitope_data,
            args=args,
            paths=paths,
            clustering_params=clustering_params,
            sample_ids=sample_ids,
        )
        sample_labels = _build_sample_labels(sample_ids, artifacts.sample_cluster_df)
    else:
        raise ValueError(f"Unsupported per-epitope cluster_algo: {clustering_params.cluster_algo}")
    sample_cluster_count = artifacts.sample_cluster_df["cluster_id"].nunique()
    enriched_cluster_count = artifacts.enriched_sample_cluster_df["cluster_id"].nunique()

    # Keep cluster_members aligned with the HTML visualization: export only
    # statistically significant/enriched sample clusters.
    cluster_members_df = build_sample_members_table(
        artifacts.enriched_sample_cluster_df, artifacts.summary_df, chain, ep_df, epitope, sample_ids, sample_umap
    )
    lfc_positive_cluster_ids = set(
        artifacts.summary_df.loc[artifacts.summary_df["log_fold_change"] > 0, "cluster_id"].astype(int)
    )
    lfc_positive_sample_cluster_df = artifacts.sample_cluster_df[
        artifacts.sample_cluster_df["cluster_id"].isin(lfc_positive_cluster_ids)
    ].copy()
    all_cluster_members_df = build_sample_members_table(
        lfc_positive_sample_cluster_df, artifacts.summary_df, chain, ep_df, epitope, sample_ids, sample_umap
    )
    clonotype_coords_df = _build_clonotype_coords_table(
        epitope=epitope,
        chain=chain,
        sample_reps=sample_reps,
        sample_ids=sample_ids,
        sample_labels=sample_labels,
        summary_df=artifacts.summary_df,
        sample_umap=sample_umap,
        clustering_params=clustering_params,
    )

    _save_cluster_results(
        artifacts.cluster_df,
        artifacts.summary_df,
        artifacts.sample_cluster_df,
        cluster_members_df,
        paths,
        prefix,
    )

    if not args.skip_umap:
        viz_filename = (
            f"{sanitize_filename_token(args.species)}_"
            f"{sanitize_filename_token(epitope)}_"
            f"{sanitize_filename_token(chain)}.html"
        )
        save_cluster_plot_html(
            epitope=epitope,
            chain=chain,
            sample_reps=sample_reps,
            sample_ids=sample_ids,
            sample_labels=sample_labels,
            significant_cluster_ids=artifacts.significant_cluster_ids,
            summary_df=artifacts.summary_df,
            sample_umap=sample_umap,
            bg_umap=bg_umap,
            output_path=paths.viz_dir / viz_filename,
        )

    logging.info(
        "Done %s: clustered_sample points=%d/%d, clusters=%d; enriched points=%d, clusters=%d",
        epitope,
        len(artifacts.sample_cluster_df),
        len(sample_reps),
        sample_cluster_count,
        len(artifacts.enriched_sample_cluster_df),
        enriched_cluster_count,
    )

    sample_cluster_df = artifacts.sample_cluster_df
    del sample_reps, sample_ids, sample_pca, sample_umap, sample_labels, artifacts
    gc.collect()
    return (
        sample_cluster_df,
        cluster_members_df,
        all_cluster_members_df,
        clonotype_coords_df,
        {
            "chain": chain,
            "epitope": epitope,
            "cluster_algo": clustering_params.cluster_algo,
            "k_neighbors": clustering_params.k_neighbors,
            "eps_k_neighbors": clustering_params.eps_k_neighbors,
            "leiden_resolution": clustering_params.leiden_resolution,
            "cluster_min_samples": clustering_params.cluster_min_samples,
            "eps_estimation_based_on": clustering_params.eps_estimation_based_on,
            "vdbscan_sym_rule": clustering_params.vdbscan_sym_rule,
            "leiden_sub_resolution": clustering_params.leiden_sub_resolution,
        },
    )


def main():
    """Main function for VDJdb epitope clustering."""
    try:
        from .compat import (
            SegmentLibrary,
            configure_logging,
            load_embedding_artifacts,
            load_prototype_repertoire,
            normalize_config,
            prepare_output_path,
            resolve_prototype_file,
            subsample_repertoire,
        )
        from .epitope_params import load_epitope_clustering_overrides
        from .io import filter_canonical_cdr3_rows, prepare_output_dirs

        args = get_arguments_vdjdb_clusters()
        args.epitope_clustering_overrides = (
            load_epitope_clustering_overrides(args.epitope_config) if args.epitope_config else {}
        )
        if args.epitope_clustering_overrides:
            logging.info(
                "Loaded per-epitope clustering overrides for %d epitope(s) from %s",
                len(args.epitope_clustering_overrides),
                args.epitope_config,
            )

        output_root = prepare_output_path(args.output)
        tcremp_cache_dir = Path(args.tcremp_cache_dir).resolve() if args.tcremp_cache_dir else None
        paths = prepare_output_dirs(output_root, tcremp_cache_dir=tcremp_cache_dir)

        configure_logging(Path(args.vdjdb), output_root, f"{args.chain.lower()}_vdjdb_clusters")
        import faiss

        faiss.omp_set_num_threads(args.nproc)

        genes = [args.chain]
        locus = CHAIN_COLS[args.chain]["locus"]
        lib = SegmentLibrary.load_default(genes=genes, organisms=args.species)

        vdjdb_df = pd.read_csv(args.vdjdb, sep="\t")
        vdjdb_df = _normalize_vdjdb_columns(vdjdb_df, args.chain)
        chain_cfg = CHAIN_COLS[args.chain]
        vdjdb_df = filter_canonical_cdr3_rows(
            vdjdb_df,
            chain_cfg["cdr3"],
            context=f"pipeline input ({args.chain})",
        )
        required_cols = ["antigen.epitope", chain_cfg["cdr3"]]
        before_chain_cleanup = len(vdjdb_df)
        vdjdb_df = vdjdb_df.dropna(subset=required_cols).copy()
        logging.info(
            "Filtered VDJdb rows for chain %s after dropping NaNs in %s: %d -> %d",
            args.chain,
            required_cols,
            before_chain_cleanup,
            len(vdjdb_df),
        )

        if args.epitopes is not None:
            vdjdb_df = vdjdb_df[vdjdb_df["antigen.epitope"].isin(args.epitopes)].copy()
        if args.min_epitope_clonotypes is not None:
            epitope_sizes = vdjdb_df.groupby("antigen.epitope").size()
            eligible_epitopes = epitope_sizes[epitope_sizes >= args.min_epitope_clonotypes].index
            skipped_epitopes = int((epitope_sizes < args.min_epitope_clonotypes).sum())
            vdjdb_df = vdjdb_df[vdjdb_df["antigen.epitope"].isin(eligible_epitopes)].copy()
            logging.info(
                "Filtered epitopes by minimum clonotype count >= %d: kept %d epitopes, skipped %d",
                args.min_epitope_clonotypes,
                len(eligible_epitopes),
                skipped_epitopes,
            )

        args.background = str(Path(args.background_airr).resolve())
        args.output = str(output_root)
        args.background_embedding = str(Path(args.background_embedding).resolve())

        proto_path = resolve_prototype_file(args.prototypes_path, chain=args.chain)

        logging.info("Loading prototypes")
        proto = load_prototype_repertoire(proto_path, lib, locus, args.index_col)
        proto = subsample_repertoire(proto, args.n_prototypes, args.sample_random_clonotypes, args.random_seed)
        logging.info("Loaded %d prototypes", len(proto))

        logging.info("Loading background embeddings")
        config = normalize_config(args)
        bg_artifacts = load_embedding_artifacts(
            path=config.background,
            args=config,
            is_sample=False,
            lib=lib,
            locus=locus,
            prefix=f"{args.chain.lower()}_background",
            output_path=paths.tcremp_dir,
        )
        bg_emb = bg_artifacts.embeddings
        bg_reps = bg_artifacts.representations
        bg_ids = bg_artifacts.ids
        bg_index_path = bg_artifacts.cache_path

        transform, transform_path, bg_pca = _stage_fit_background_transform(
            args=args,
            output_root=output_root,
            bg_emb=bg_emb,
        )

        # Fit the background-only transform before the heavy per-epitope
        # embedding stage. This avoids the cold-start path where a long
        # embedding precompute run can leave the subsequent transform step
        # appearing stalled on the very first launch.
        epitope_infos = _stage_precompute_epitope_embeddings(
            vdjdb_df,
            args=args,
            genes=genes,
            locus=locus,
            lib=lib,
            proto=proto,
            paths=paths,
        )

        if args.skip_umap:
            logging.info("Skipping Stage 3/4 joint plotting UMAP because --skip-umap was requested")
            epitope_inputs = []
            for epitope_info in epitope_infos:
                epitope_data = _load_epitope_input(
                    epitope_info,
                    args=args,
                    lib=lib,
                    locus=locus,
                    paths=paths,
                )
                sample_pca = transform.transform_pca(epitope_data["sample_emb"])
                epitope_data["sample_pca"] = sample_pca
                epitope_inputs.append(epitope_data)
                del epitope_data["sample_emb"]
            bg_umap = None
        else:
            epitope_inputs, bg_umap = _stage_prepare_joint_umap(
                epitope_infos,
                args=args,
                lib=lib,
                locus=locus,
                paths=paths,
                transform=transform,
                bg_pca=bg_pca,
                transform_path=transform_path,
            )

        (
            clustered_tables,
            cluster_members_tables,
            all_cluster_members_tables,
            clonotype_coords_tables,
            parameter_rows,
        ) = _stage_run_per_epitope_analysis(
            epitope_inputs,
            args=args,
            paths=paths,
            bg_pca=bg_pca,
            bg_reps=bg_reps,
            bg_ids=bg_ids,
            bg_index_path=bg_index_path,
            bg_umap=bg_umap,
        )
        del epitope_inputs, epitope_infos
        gc.collect()

        chain_lower = args.chain.lower()
        if clustered_tables:
            pd.concat(clustered_tables, ignore_index=True).to_csv(
                paths.tcrempnet_dir / f"{chain_lower}_vdjdb_clustered_clonotypes.tsv", sep="\t", index=False
            )
        if cluster_members_tables:
            pd.concat(cluster_members_tables, ignore_index=True).to_csv(
                output_root / _cluster_members_filename(args.chain, args.output_tag), sep="\t", index=False
            )
        if all_cluster_members_tables:
            pd.concat(all_cluster_members_tables, ignore_index=True).to_csv(
                output_root / _cluster_members_filename(args.chain, args.output_tag, prefix="cluster_members_all"),
                sep="\t",
                index=False,
            )
        if clonotype_coords_tables:
            pd.concat(clonotype_coords_tables, ignore_index=True).to_csv(
                output_root / f"{chain_lower}_vdjdb_clonotype_coords_2d.tsv",
                sep="\t",
                index=False,
            )
        if parameter_rows:
            pd.DataFrame(parameter_rows).to_csv(
                output_root / f"{chain_lower}_vdjdb_epitope_clustering_params.tsv",
                sep="\t",
                index=False,
            )
        if bg_umap is not None:
            _build_background_coords_table(chain=args.chain, bg_umap=bg_umap).to_csv(
                output_root / f"{chain_lower}_background_coords_2d.tsv",
                sep="\t",
                index=False,
            )
        del (
            clustered_tables,
            cluster_members_tables,
            all_cluster_members_tables,
            clonotype_coords_tables,
            parameter_rows,
            bg_umap,
            bg_pca,
            bg_reps,
            bg_ids,
        )
        gc.collect()

        logging.info("Done")
    except Exception as e:
        logging.error("An error occurred: %s", e)
        raise


if __name__ == "__main__":
    main()
