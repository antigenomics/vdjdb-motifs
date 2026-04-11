from __future__ import annotations

import colorsys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.colors import hex_to_rgb, sample_colorscale, unlabel_rgb

from .config import CHAIN_COLS


def _format_cluster_display_label(*, chain: str, epitope: str, cluster_id: int) -> str:
    """Format cluster labels exactly like ``cid`` in final cluster member tables."""
    chain_tag = "A" if chain == "TRA" else "B"
    return f"H.{chain_tag}.{epitope}.{int(cluster_id)}"


def _parse_plotly_color(color: str) -> tuple[int, int, int]:
    """Parse Plotly ``#hex`` and ``rgb(...)`` colors into integer RGB."""
    if color.startswith("#"):
        return hex_to_rgb(color)
    if color.startswith("rgb"):
        return tuple(int(round(channel)) for channel in unlabel_rgb(color))
    raise ValueError(f"Unsupported Plotly color format: {color}")


def _is_non_gray_color(color: str, *, min_saturation: float = 0.28) -> bool:
    """Keep only sufficiently saturated colors to avoid gray-ish cluster markers."""
    red, green, blue = _parse_plotly_color(color)
    hue, lightness, saturation = colorsys.rgb_to_hls(red / 255, green / 255, blue / 255)
    _ = hue, lightness
    return saturation >= min_saturation


def _build_cluster_palette(n_colors: int) -> list[str]:
    """Return a large non-gray palette sized to the number of displayed clusters."""
    if n_colors <= 0:
        return []

    palette_names = ("Alphabet", "Light24", "Dark24", "Safe", "Vivid", "Prism")
    palette_pool: list[str] = []
    seen_colors: set[str] = set()

    for palette_name in palette_names:
        for color in getattr(px.colors.qualitative, palette_name):
            if not _is_non_gray_color(color):
                continue
            if color in seen_colors:
                continue
            palette_pool.append(color)
            seen_colors.add(color)

    if n_colors <= len(palette_pool):
        return palette_pool[:n_colors]

    extra_needed = n_colors - len(palette_pool)
    if extra_needed == 1:
        extra_positions = [0.5]
    else:
        extra_positions = np.linspace(0.0, 1.0, extra_needed, endpoint=False).tolist()
    extra_colors = [
        color
        for color in sample_colorscale("Turbo", extra_positions, colortype="rgb")
        if _is_non_gray_color(color) and color not in seen_colors
    ]
    return (palette_pool + extra_colors)[:n_colors]


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
    enriched_palette = _build_cluster_palette(len(ordered_cluster_labels))
    for cluster_label, color in zip(ordered_cluster_labels, enriched_palette):
        color_discrete_map[cluster_label] = color

    hover_cols = {
        col: True
        for col in (
            f"cdr3aa_{cfg['gene']}",
            f"v_{cfg['gene']}",
            f"j_{cfg['gene']}",
            "clone_id",
            "cluster",
            "cluster_id",
            "cluster_size_sample",
        )
        if col in sample_df.columns
    }
    hover_cols["x"] = False
    hover_cols["y"] = False
    hover_cols["cluster_label"] = False
    hover_cols["significant"] = False
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
