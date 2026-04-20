#!/usr/bin/env python
from __future__ import annotations

import gc
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .arguments import get_arguments_vdjdb_clusters
from .config import CHAIN_COLS, DEFAULT_PLOT_BG_POINTS


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


def _cluster_members_filename(chain: str, output_tag: str | None) -> str:
    if output_tag is None:
        return f"cluster_members_{chain}.txt"
    clean_tag = "".join(char if char.isalnum() or char in {"-", "_"} else "_" for char in output_tag).strip("_")
    if not clean_tag:
        return f"cluster_members_{chain}.txt"
    return f"cluster_members_{chain}_{clean_tag}.txt"


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
    args,
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
        args.k_neighbors,
        args.nproc,
        sample_index_path=sample_index_path,
        bg_index_path=bg_index_path,
    )
    labels = run_leiden_clustering(
        knn_indices=indices,
        knn_distances=distances,
        resolution=args.leiden_resolution,
        n_threads=args.nproc,
        min_cluster_size=args.cluster_min_samples,
        min_cluster_size_mask=np.arange(len(sample_pca) + len(bg_pca)) < len(sample_pca),
    )
    return sample_pca, labels


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
    from .io import build_airr_from_epitope

    chain = args.chain
    prefix = f"{chain.lower()}_vdjdb_{epitope}"
    airr_path = paths.airr_dir / f"{prefix}.tsv"

    logging.info("Preparing epitope %s with %d clonotypes", epitope, len(ep_df))
    build_airr_from_epitope(ep_df, chain).to_csv(airr_path, sep="\t", index=False)
    sample_artifacts = _compute_sample_embeddings(
        args, genes, locus, lib, proto, paths, chain, prefix, airr_path
    )
    return {
        "epitope": epitope,
        "ep_df": ep_df,
        "prefix": prefix,
        "airr_path": airr_path,
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
    for epitope, ep_df in vdjdb_df.groupby("antigen.epitope", sort=True):
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
) -> tuple[list[pd.DataFrame], list[pd.DataFrame]]:
    """Stage 4: run per-epitope clustering and save outputs."""
    logging.info("Stage 4/4: running per-epitope clustering analysis")
    clustered_tables: list[pd.DataFrame] = []
    cluster_members_tables: list[pd.DataFrame] = []
    for epitope_data in epitope_inputs:
        clustered_df, cluster_members_df = process_epitope(
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
    logging.info("Stage 4/4 done")
    return clustered_tables, cluster_members_tables


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
) -> tuple[pd.DataFrame, pd.DataFrame]:
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
        Tuple of clustered DataFrame and cluster members DataFrame.
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

    logging.info("Processing epitope %s with %d clonotypes", epitope, len(ep_df))
    sample_pca, labels = _perform_clustering(sample_pca, bg_pca, args, sample_index_path, bg_index_path, sample_ids, bg_ids)

    artifacts = _build_epitope_clustering_artifacts(
        labels=labels,
        sample_ids=sample_ids,
        bg_ids=bg_ids,
        sample_reps=sample_reps,
        bg_reps=bg_reps,
    )
    sample_cluster_count = artifacts.sample_cluster_df["cluster_id"].nunique()
    enriched_cluster_count = artifacts.enriched_sample_cluster_df["cluster_id"].nunique()

    # Keep cluster_members aligned with the HTML visualization: export only
    # statistically significant/enriched sample clusters.
    cluster_members_df = build_sample_members_table(
        artifacts.enriched_sample_cluster_df, artifacts.summary_df, chain, ep_df, epitope, sample_ids, sample_umap
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
        sample_labels = np.asarray(labels[: len(sample_ids)], dtype=np.int32)
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
    del sample_reps, sample_ids, sample_pca, sample_umap, labels, artifacts
    gc.collect()
    return sample_cluster_df, cluster_members_df


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
        from .io import prepare_output_dirs

        args = get_arguments_vdjdb_clusters()

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

        clustered_tables, cluster_members_tables = _stage_run_per_epitope_analysis(
            epitope_inputs,
            args=args,
            paths=paths,
            bg_pca=bg_pca,
            bg_reps=bg_reps,
            bg_ids=bg_ids,
            bg_index_path=bg_index_path,
            bg_umap=bg_umap,
        )
        del epitope_inputs, epitope_infos, bg_umap
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
        del clustered_tables, cluster_members_tables, bg_pca, bg_reps, bg_ids
        gc.collect()

        logging.info("Done")
    except Exception as e:
        logging.error("An error occurred: %s", e)
        raise


if __name__ == "__main__":
    main()
