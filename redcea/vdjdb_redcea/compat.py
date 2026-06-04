from __future__ import annotations

from mir.common.segments import SegmentLibrary

from redcea.analysis.background_transform import BackgroundTransform
from redcea.analysis.cluster_utils import compute_cluster_summary
from redcea.analysis.io import load_embedding_artifacts
from redcea.clustering import build_joint_knn_from_split, compute_blockwise_knn_merged, run_leiden_clustering
from redcea.config import PipelineConfig
from redcea.embeddings import compute_embeddings_if_needed
from redcea.utils.paths import resolve_prototype_file
from redcea.utils.stats import add_log_fold_change, add_z_binom_pvalues
from redcea.utils.tcremp import (
    configure_logging,
    load_prototype_repertoire,
    prepare_output_path,
    subsample_repertoire,
)


def normalize_config(args) -> PipelineConfig:
    """Build shared RedCEA config from the vdjdb_redcea argparse namespace."""
    return PipelineConfig.from_args(args)

__all__ = [
    "BackgroundTransform",
    "SegmentLibrary",
    "add_log_fold_change",
    "add_z_binom_pvalues",
    "build_joint_knn_from_split",
    "compute_blockwise_knn_merged",
    "compute_cluster_summary",
    "compute_embeddings_if_needed",
    "configure_logging",
    "load_embedding_artifacts",
    "load_prototype_repertoire",
    "normalize_config",
    "prepare_output_path",
    "resolve_prototype_file",
    "run_leiden_clustering",
    "subsample_repertoire",
]
