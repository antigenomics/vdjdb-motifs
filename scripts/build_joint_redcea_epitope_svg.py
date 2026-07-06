from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
MPL_CONFIG_DIR = REPO_ROOT / ".tmp" / "matplotlib"
MPL_CONFIG_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CONFIG_DIR))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from matplotlib import colors as mcolors
from matplotlib import lines as mlines
from scipy.ndimage import gaussian_filter
from scipy.stats import gaussian_kde


@dataclass(frozen=True)
class EpitopeStyle:
    outline: str
    fill_low: str
    fill_high: str
    density_cmap: mcolors.Colormap
    kde_bw: float
    fill_alpha: float
    line_alpha: float
    min_island_points: int


DEFAULT_STYLES: dict[str, EpitopeStyle] = {
    "GLCTLVAML": EpitopeStyle(
        outline="#6baed6",
        fill_low="#deebf7",
        fill_high="#2166ac",
        density_cmap=mcolors.LinearSegmentedColormap.from_list(
            "glc_blue",
            ["#deebf7", "#9ecae1", "#4292c6", "#2166ac"],
        ),
        kde_bw=0.08,
        fill_alpha=0.42,
        line_alpha=0.72,
        min_island_points=2,
    ),
    "YLQPRTFLL": EpitopeStyle(
        outline="#ef8a62",
        fill_low="#fddbc7",
        fill_high="#b2182b",
        density_cmap=mcolors.LinearSegmentedColormap.from_list(
            "ylq_orange_red",
            ["#fddbc7", "#f4a582", "#d6604d", "#b2182b"],
        ),
        kde_bw=0.12,
        fill_alpha=0.68,
        line_alpha=0.88,
        min_island_points=2,
    ),
}


def parse_args() -> argparse.Namespace:
    default_input_root = REPO_ROOT / "results" / "redcea_ylq_glc_joint_umap"
    default_output_path = default_input_root / "viz" / "joint_redcea_epitope_overlay.svg"
    default_metadata_path = default_input_root / "viz" / "joint_redcea_epitope_overlay.metadata.json"

    parser = argparse.ArgumentParser(
        description=(
            "Build one SVG overlay for the REDCEA joint UMAP run with background "
            "contours and epitope-specific point/cluster contours."
        )
    )
    parser.add_argument("--input-root", type=Path, default=default_input_root)
    parser.add_argument("--coords-file", type=Path, default=None)
    parser.add_argument("--background-file", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=default_output_path)
    parser.add_argument("--metadata-output", type=Path, default=default_metadata_path)
    parser.add_argument("--bins", type=int, default=240, help="Histogram bins per axis for contour estimation.")
    parser.add_argument("--background-sigma", type=float, default=2.4)
    parser.add_argument("--background-levels", type=float, nargs="+", default=[0.10, 0.18, 0.30, 0.46, 0.66])
    parser.add_argument("--sample-levels", type=float, nargs="+", default=[0.025, 0.05, 0.09, 0.15, 0.24, 0.36, 0.52])
    parser.add_argument("--pvalue-threshold", type=float, default=0.05)
    parser.add_argument("--min-log-fold-change", type=float, default=0.0)
    return parser.parse_args()


def resolve_input_paths(args: argparse.Namespace) -> tuple[Path, Path]:
    coords_file = args.coords_file or (args.input_root / "trb_vdjdb_clonotype_coords_2d.tsv")
    background_file = args.background_file or (args.input_root / "trb_background_coords_2d.tsv")
    if not coords_file.exists():
        raise FileNotFoundError(f"Coordinate file not found: {coords_file}")
    if not background_file.exists():
        raise FileNotFoundError(f"Background file not found: {background_file}")
    return coords_file, background_file


def find_summary_file(input_root: Path, epitope: str) -> Path:
    pattern = f"per_epitope_cluster_runs/trb_vdjdb_{epitope}/**/*_summary_tcrempnet.tsv"
    matches = sorted(input_root.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No summary_tcrempnet file found for epitope {epitope!r} under {input_root}")
    return matches[0]


def load_inputs(input_root: Path, coords_file: Path, background_file: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, pd.DataFrame]]:
    coords = pd.read_csv(coords_file, sep="\t")
    background = pd.read_csv(background_file, sep="\t")
    summaries: dict[str, pd.DataFrame] = {}
    for epitope in sorted(coords["epitope"].dropna().unique()):
        summary = pd.read_csv(find_summary_file(input_root, epitope), sep="\t")
        summaries[epitope] = summary.rename(columns={"cluster_id": "cluster_id_summary"})
    return coords, background, summaries


def attach_cluster_statistics(coords: pd.DataFrame, summaries: dict[str, pd.DataFrame]) -> pd.DataFrame:
    merged_parts: list[pd.DataFrame] = []
    for epitope, ep_coords in coords.groupby("epitope", sort=True):
        summary = summaries[epitope].rename(columns={"cluster_id_summary": "cluster_id"})
        summary = summary[
            [
                "cluster_id",
                "cluster_size",
                "sample",
                "background",
                "log_fold_change",
                "enrichment_pvalue_zbinom",
                "enrichment_fdr_zbinom",
            ]
        ].copy()
        merged = ep_coords.merge(summary, on="cluster_id", how="left", suffixes=("", "_summary"))
        if "log_fold_change_summary" in merged.columns:
            merged["log_fold_change"] = merged["log_fold_change"].fillna(merged["log_fold_change_summary"])
            merged = merged.drop(columns=["log_fold_change_summary"])
        merged_parts.append(merged)
    return pd.concat(merged_parts, ignore_index=True)


def flag_enriched_clusters(coords: pd.DataFrame, pvalue_threshold: float, min_log_fold_change: float) -> pd.DataFrame:
    coords = coords.copy()
    coords["is_enriched_cluster"] = (
        coords["cluster_id"].ge(0)
        & coords["enrichment_pvalue_zbinom"].fillna(np.inf).lt(pvalue_threshold)
        & coords["log_fold_change"].fillna(-np.inf).gt(min_log_fold_change)
    )
    return coords


def build_grid(
    coords: pd.DataFrame,
    background: pd.DataFrame,
    bins: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    combined_x = np.concatenate([coords["x"].to_numpy(dtype=float), background["x"].to_numpy(dtype=float)])
    combined_y = np.concatenate([coords["y"].to_numpy(dtype=float), background["y"].to_numpy(dtype=float)])
    x_pad = (combined_x.max() - combined_x.min()) * 0.06 or 1.0
    y_pad = (combined_y.max() - combined_y.min()) * 0.06 or 1.0
    x_edges = np.linspace(combined_x.min() - x_pad, combined_x.max() + x_pad, bins + 1)
    y_edges = np.linspace(combined_y.min() - y_pad, combined_y.max() + y_pad, bins + 1)
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2.0
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2.0
    xx, yy = np.meshgrid(x_centers, y_centers)
    return x_edges, y_edges, xx, yy


def smooth_histogram(
    x: np.ndarray,
    y: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    sigma: float,
) -> np.ndarray:
    hist, _, _ = np.histogram2d(x, y, bins=[x_edges, y_edges])
    return gaussian_filter(hist.T, sigma=sigma, mode="nearest")


def scaled_levels(density: np.ndarray, fractions: list[float]) -> list[float]:
    if density.size == 0:
        return []
    max_value = float(np.nanmax(density))
    if max_value <= 0:
        return []
    levels = sorted({max_value * fraction for fraction in fractions if 0 < fraction < 1})
    return [level for level in levels if level < max_value]


def kde_density(
    x: np.ndarray,
    y: np.ndarray,
    xx: np.ndarray,
    yy: np.ndarray,
    bandwidth_factor: float,
) -> np.ndarray:
    values = np.vstack([x, y])
    if values.shape[1] < 2:
        return np.zeros_like(xx, dtype=float)
    if np.allclose(np.var(x), 0.0) and np.allclose(np.var(y), 0.0):
        return np.zeros_like(xx, dtype=float)
    kde = gaussian_kde(values, bw_method=bandwidth_factor)
    density = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
    return density


def connected_components(points: np.ndarray, radius: float) -> list[np.ndarray]:
    if len(points) == 0:
        return []
    visited = np.zeros(len(points), dtype=bool)
    components: list[np.ndarray] = []
    for start_idx in range(len(points)):
        if visited[start_idx]:
            continue
        stack = [start_idx]
        visited[start_idx] = True
        component: list[int] = []
        while stack:
            idx = stack.pop()
            component.append(idx)
            deltas = points - points[idx]
            neighbors = np.where((deltas[:, 0] ** 2 + deltas[:, 1] ** 2) <= radius**2)[0]
            for neighbor_idx in neighbors:
                if not visited[neighbor_idx]:
                    visited[neighbor_idx] = True
                    stack.append(int(neighbor_idx))
        components.append(np.asarray(component, dtype=int))
    return components


def draw_component_fallback_contours(
    ax: plt.Axes,
    points: np.ndarray,
    style: EpitopeStyle,
) -> int:
    if len(points) == 0:
        return 0
    span_x = float(points[:, 0].max() - points[:, 0].min()) if len(points) > 1 else 0.0
    span_y = float(points[:, 1].max() - points[:, 1].min()) if len(points) > 1 else 0.0
    patch = mpatches.Ellipse(
        (float(points[:, 0].mean()), float(points[:, 1].mean())),
        width=max(0.7, span_x + 0.8),
        height=max(0.7, span_y + 0.8),
        facecolor=mcolors.to_rgba(style.fill_low, alpha=style.fill_alpha * 0.45),
        edgecolor=mcolors.to_rgba(style.fill_high, alpha=style.line_alpha),
        linewidth=1.0,
        zorder=3.05,
    )
    ax.add_patch(patch)
    return 1


def build_background_cmap() -> mcolors.Colormap:
    return mcolors.LinearSegmentedColormap.from_list(
        "soft_background_grey",
        ["#ffffff", "#ffffff", "#e0e0e0", "#b8b8b8", "#7a7a7a"],
    )


def draw_background_contours(
    ax: plt.Axes,
    background: pd.DataFrame,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    xx: np.ndarray,
    yy: np.ndarray,
    sigma: float,
    level_fractions: list[float],
) -> None:
    density = smooth_histogram(
        background["x"].to_numpy(dtype=float),
        background["y"].to_numpy(dtype=float),
        x_edges,
        y_edges,
        sigma=sigma,
    )
    levels = scaled_levels(density, level_fractions)
    if levels:
        ax.contourf(
            xx,
            yy,
            density,
            levels=[0.0] + levels + [float(np.nanmax(density))],
            cmap=build_background_cmap(),
            alpha=1.0,
            antialiased=True,
            zorder=0,
        )
        ax.contour(
            xx,
            yy,
            density,
            levels=levels,
            colors=["#b3b3b3"] * len(levels),
            linewidths=np.linspace(0.75, 1.25, len(levels)),
            alpha=0.95,
            zorder=1,
        )


def draw_epitope_points(
    ax: plt.Axes,
    epitope_df: pd.DataFrame,
    style: EpitopeStyle,
) -> None:
    non_clustered = epitope_df[~epitope_df["is_enriched_cluster"]].copy()
    if not non_clustered.empty:
        ax.scatter(
            non_clustered["x"],
            non_clustered["y"],
            s=66,
            facecolors=mcolors.to_rgba(style.fill_low, alpha=0.62),
            edgecolors=mcolors.to_rgba(style.outline, alpha=0.82),
            linewidths=0.95,
            zorder=2,
        )

    clustered = epitope_df[epitope_df["is_enriched_cluster"]].copy()
    if clustered.empty:
        return
    ax.scatter(
        clustered["x"],
        clustered["y"],
        s=60,
        c=style.fill_high,
        edgecolors=mcolors.to_rgba("#202020", alpha=0.78),
        linewidths=0.7,
        alpha=0.98,
        zorder=4,
    )


def draw_sample_density_contours(
    ax: plt.Axes,
    epitope_df: pd.DataFrame,
    style: EpitopeStyle,
    xx: np.ndarray,
    yy: np.ndarray,
    level_fractions: list[float],
) -> dict[str, float]:
    clustered_df = epitope_df[epitope_df["is_enriched_cluster"]].copy()
    if clustered_df.empty:
        return {"density_max": 0.0, "drawn": 0, "points_used": 0, "fallback_islands": 0}

    sample_density = kde_density(
        clustered_df["x"].to_numpy(dtype=float),
        clustered_df["y"].to_numpy(dtype=float),
        xx,
        yy,
        bandwidth_factor=style.kde_bw,
    )
    density_max = float(np.nanmax(sample_density))
    if density_max <= 0:
        return {"density_max": density_max, "drawn": 0, "points_used": int(len(clustered_df)), "fallback_islands": 0}
    levels = sorted({density_max * fraction for fraction in level_fractions if 0 < fraction < 1})
    levels = [level for level in levels if 0 < level < density_max]
    if not levels:
        return {"density_max": density_max, "drawn": 0, "points_used": int(len(clustered_df)), "fallback_islands": 0}

    fill_levels = levels + [density_max]
    ax.contourf(
        xx,
        yy,
        sample_density,
        levels=fill_levels,
        cmap=style.density_cmap,
        alpha=style.fill_alpha,
        antialiased=True,
        zorder=2.8,
    )
    ax.contour(
        xx,
        yy,
        sample_density,
        levels=levels,
        colors=[mcolors.to_rgba(style.fill_high, alpha=style.line_alpha)] * len(levels),
        linewidths=np.linspace(0.7, 1.2, len(levels)),
        zorder=3.2,
    )
    grid_x = xx[0]
    grid_y = yy[:, 0]
    x_idx = np.clip(np.searchsorted(grid_x, clustered_df["x"].to_numpy(dtype=float)), 0, len(grid_x) - 1)
    y_idx = np.clip(np.searchsorted(grid_y, clustered_df["y"].to_numpy(dtype=float)), 0, len(grid_y) - 1)
    point_density = sample_density[y_idx, x_idx]
    missing_points = clustered_df.loc[point_density < levels[0], ["x", "y"]].to_numpy(dtype=float)
    fallback_islands = 0
    for component in connected_components(missing_points, radius=1.2):
        component_points = missing_points[component]
        if len(component_points) >= style.min_island_points:
            fallback_islands += draw_component_fallback_contours(ax, component_points, style)
    return {
        "density_max": density_max,
        "drawn": 1,
        "points_used": int(len(clustered_df)),
        "fallback_islands": int(fallback_islands),
    }


def add_legends(
    fig: plt.Figure,
    ax: plt.Axes,
    epitope_order: list[str],
    styles: dict[str, EpitopeStyle],
) -> None:
    handles = [
        mlines.Line2D([], [], color="#d0d0d0", linewidth=1.0, label="Background contours"),
    ]
    for epitope in epitope_order:
        style = styles[epitope]
        handles.append(
            mlines.Line2D(
                [],
                [],
                color=style.fill_low,
                marker="o",
                markerfacecolor=style.fill_low,
                markeredgecolor=style.fill_low,
                markeredgewidth=0.0,
                linewidth=0.0,
                label=f"{epitope} sample outside clusters",
            )
        )
        handles.append(
            mlines.Line2D(
                [],
                [],
                color=style.fill_high,
                marker="o",
                markerfacecolor=style.fill_high,
                markeredgecolor=style.fill_high,
                markeredgewidth=0.0,
                linewidth=0.0,
                label=f"{epitope} sample in clusters",
            )
        )
        handles.append(
            mlines.Line2D(
                [],
                [],
                color=style.fill_high,
                linewidth=1.4,
                label=f"{epitope} sample density",
            )
        )
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=10)


def save_metadata(
    metadata_path: Path,
    *,
    args: argparse.Namespace,
    styles: dict[str, EpitopeStyle],
    coords_file: Path,
    background_file: Path,
    summary_files: dict[str, Path],
    contour_stats: dict[str, dict[str, int]],
    epitope_counts: dict[str, dict[str, int]],
    sample_density_stats: dict[str, dict[str, float]],
) -> None:
    payload = {
        "input_root": str(args.input_root.resolve()),
        "coords_file": str(coords_file.resolve()),
        "background_file": str(background_file.resolve()),
        "summary_files": {epitope: str(path.resolve()) for epitope, path in summary_files.items()},
        "output": str(args.output.resolve()),
        "bins": args.bins,
        "background_sigma": args.background_sigma,
        "epitope_kde_bw": {epitope: styles[epitope].kde_bw for epitope in styles},
        "epitope_fill_alpha": {epitope: styles[epitope].fill_alpha for epitope in styles},
        "epitope_line_alpha": {epitope: styles[epitope].line_alpha for epitope in styles},
        "background_levels": args.background_levels,
        "sample_levels": args.sample_levels,
        "pvalue_threshold": args.pvalue_threshold,
        "min_log_fold_change": args.min_log_fold_change,
        "contour_stats": contour_stats,
        "sample_density_stats": sample_density_stats,
        "epitope_counts": epitope_counts,
    }
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main() -> None:
    args = parse_args()
    coords_file, background_file = resolve_input_paths(args)
    coords_raw, background, summaries = load_inputs(args.input_root, coords_file, background_file)
    coords = attach_cluster_statistics(coords_raw, summaries)
    coords = flag_enriched_clusters(coords, args.pvalue_threshold, args.min_log_fold_change)

    x_edges, y_edges, xx, yy = build_grid(coords, background, bins=args.bins)
    epitope_order = sorted(coords["epitope"].dropna().unique().tolist())
    styles = {epitope: DEFAULT_STYLES.get(epitope, DEFAULT_STYLES["GLCTLVAML"]) for epitope in epitope_order}

    fig, ax = plt.subplots(figsize=(14, 11))
    fig.subplots_adjust(left=0.05, right=0.80, top=0.94, bottom=0.05)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    draw_background_contours(
        ax,
        background,
        x_edges,
        y_edges,
        xx,
        yy,
        sigma=args.background_sigma,
        level_fractions=list(args.background_levels),
    )

    contour_stats: dict[str, dict[str, int]] = {}
    epitope_counts: dict[str, dict[str, int]] = {}
    sample_density_stats: dict[str, dict[str, float]] = {}
    summary_files: dict[str, Path] = {}
    for epitope in epitope_order:
        epitope_df = coords[coords["epitope"].eq(epitope)].copy()
        summary_files[epitope] = find_summary_file(args.input_root, epitope)
        sample_density_stats[epitope] = draw_sample_density_contours(
            ax,
            epitope_df,
            styles[epitope],
            xx,
            yy,
            level_fractions=list(args.sample_levels),
        )
        draw_epitope_points(ax, epitope_df, styles[epitope])
        contour_stats[epitope] = {"cluster_points_enriched": int(epitope_df["is_enriched_cluster"].sum())}
        epitope_counts[epitope] = {
            "sample_points": int(len(epitope_df)),
            "cluster_points_any": int(epitope_df["cluster_id"].ge(0).sum()),
            "cluster_points_enriched": int(epitope_df["is_enriched_cluster"].sum()),
            "cluster_ids_any": int(epitope_df.loc[epitope_df["cluster_id"].ge(0), "cluster_id"].nunique()),
            "cluster_ids_enriched": int(epitope_df.loc[epitope_df["is_enriched_cluster"], "cluster_id"].nunique()),
        }

    add_legends(fig, ax, epitope_order, styles)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal", adjustable="box")
    for spine in ax.spines.values():
        spine.set_visible(False)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, format="svg", dpi=300, bbox_inches="tight", transparent=False)
    plt.close(fig)

    save_metadata(
        args.metadata_output,
        args=args,
        styles=styles,
        coords_file=coords_file,
        background_file=background_file,
        summary_files=summary_files,
        contour_stats=contour_stats,
        epitope_counts=epitope_counts,
        sample_density_stats=sample_density_stats,
    )
    print(f"Saved SVG: {args.output}")
    print(f"Saved metadata: {args.metadata_output}")


if __name__ == "__main__":
    main()
