from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from .config import CHAIN_COLS


def _format_cluster_display_label(*, chain: str, epitope: str, cluster_id: int) -> str:
    """Format cluster labels exactly like ``cid`` in final cluster member tables."""
    chain_tag = "A" if chain == "TRA" else "B"
    return f"H.{chain_tag}.{epitope}.{int(cluster_id)}"


def build_cluster_plot(
    *,
    epitope: str,
    chain: str,
    sample_reps: pd.DataFrame,
    sample_ids: pd.Series,
    sample_labels: np.ndarray,
    significant_cluster_ids: set[int],
    summary_df: pd.DataFrame,
    sample_umap: np.ndarray,
    bg_umap: np.ndarray,
) -> go.Figure:
    """Build a Plotly figure for cluster visualization.

    Args:
        epitope: Epitope name.
        chain: Chain type.
        sample_reps: Sample representations.
        sample_ids: Sample IDs.
        sample_labels: Cluster labels for samples.
        significant_cluster_ids: Set of significant cluster IDs.
        summary_df: Summary DataFrame.
        sample_umap: Sample UMAP coordinates.
        bg_umap: Background UMAP coordinates.

    Returns:
        Plotly Figure object.
    """
    cfg = CHAIN_COLS[chain]
    sample_df = sample_reps.copy().reset_index(drop=True)
    sample_df["clone_id"] = sample_ids.to_numpy()
    sample_df["cluster_id"] = sample_labels
    sample_df["x"] = sample_umap[:, 0]
    sample_df["y"] = sample_umap[:, 1]
    sample_df["significant"] = sample_df["cluster_id"].isin(significant_cluster_ids)
    summary_by_cluster = summary_df.set_index("cluster_id")
    sample_df["cluster_size_sample"] = sample_df["cluster_id"].map(summary_by_cluster["sample"])
    sample_df["log_fold_change"] = sample_df["cluster_id"].map(summary_by_cluster["log_fold_change"])
    sample_df["cluster_label"] = sample_df["cluster_id"].map(
        lambda cluster_id: _format_cluster_display_label(chain=chain, epitope=epitope, cluster_id=cluster_id)
    )
    sample_df["cluster"] = np.where(sample_df["significant"], sample_df["cluster_label"], "unclustered")
    ordered_clusters = sorted(sample_df.loc[sample_df["significant"], "cluster_id"].unique().tolist())
    ordered_cluster_labels = [
        _format_cluster_display_label(chain=chain, epitope=epitope, cluster_id=cluster_id)
        for cluster_id in ordered_clusters
    ]
    category_orders = {"cluster": ["unclustered"] + ordered_cluster_labels}

    color_discrete_map = {"unclustered": "lightgrey"}
    enriched_palette = px.colors.qualitative.Plotly
    for i, cluster_label in enumerate(ordered_cluster_labels):
        color_discrete_map[cluster_label] = enriched_palette[i % len(enriched_palette)]

    hover_cols = {
        col: True
        for col in (
            f"cdr3aa_{cfg['gene']}",
            f"v_{cfg['gene']}",
            f"j_{cfg['gene']}",
            "clone_id",
            "cluster",
            "cluster_label",
            "cluster_id",
            "cluster_size_sample",
            "significant",
        )
        if col in sample_df.columns
    }
    if "log_fold_change" in sample_df.columns:
        hover_cols["log_fold_change"] = ":.2f"

    fig_scatter = px.scatter(
        sample_df,
        x="x",
        y="y",
        color="cluster",
        category_orders=category_orders,
        color_discrete_map=color_discrete_map,
        hover_data=hover_cols,
    )

    fig = go.Figure()
    fig.add_trace(
        go.Histogram2dContour(
            x=bg_umap[:, 0],
            y=bg_umap[:, 1],
            ncontours=20,
            contours=dict(coloring="fill", showlines=False),
            colorscale="Greys",
            showscale=False,
            hoverinfo="skip",
            opacity=0.35,
        )
    )
    for trace in fig_scatter.data:
        fig.add_trace(trace)

    fig.update_layout(
        title=f"TCR clustering with background density ({epitope})",
        width=1000,
        height=700,
        template="plotly_white",
    )
    return fig


def save_cluster_plot_html(
    *,
    epitope: str,
    chain: str,
    sample_reps: pd.DataFrame,
    sample_ids: pd.Series,
    sample_labels: np.ndarray,
    significant_cluster_ids: set[int],
    summary_df: pd.DataFrame,
    sample_umap: np.ndarray,
    bg_umap: np.ndarray,
    output_path: Path,
) -> None:
    """Save cluster plot as HTML file.

    Args:
        epitope: Epitope name.
        chain: Chain type.
        sample_reps: Sample representations.
        sample_ids: Sample IDs.
        sample_labels: Cluster labels.
        significant_cluster_ids: Significant cluster IDs.
        summary_df: Summary DataFrame.
        sample_umap: Sample UMAP coordinates.
        bg_umap: Background UMAP.
        output_path: Path to save HTML.
    """
    fig = build_cluster_plot(
        epitope=epitope,
        chain=chain,
        sample_reps=sample_reps,
        sample_ids=sample_ids,
        sample_labels=sample_labels,
        significant_cluster_ids=significant_cluster_ids,
        summary_df=summary_df,
        sample_umap=sample_umap,
        bg_umap=bg_umap,
    )
    output_path.write_text(fig.to_html(full_html=True, include_plotlyjs="cdn"), encoding="utf-8")
