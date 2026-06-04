#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import logomaker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import pairwise2
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = ROOT / "figures"
VDJDB_SLIM_PATH = ROOT / "vdjdb_release" / "vdjdb.slim.txt"
TCRNET_CLUSTER_PATH = ROOT / "results" / "tcrnet" / "cluster_members.txt"
REDCEA_TRA_PATH = ROOT / "results" / "redcea" / "cluster_members_TRA.txt"
REDCEA_TRB_PATH = ROOT / "results" / "redcea" / "cluster_members_TRB.txt"
OUTPUT_BASENAME = FIGURE_DIR / "tcrnet_tcrempnet_redcea_final_one_panel"

POINT_KEY = ["species", "gene", "antigen.epitope", "cdr3aa", "v.segm", "j.segm"]
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")
CHAIN = "TRB"
EPITOPE = "YLQPRTFLL"
PREFIX = "CSAR"
COLORS = {
    "TCRnet": "#4c78a8",
    "RedCEA": "#f58518",
    "VDJdb": "#9e9e9e",
}
TCRNET_CLUSTER_COLORS = {
    "C21": "#72c7b0",
    "C24": "#fc946a",
}
FONT_DELTA = 3


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build the final TCRnet vs RedCEA one-panel figure from scratch."
    )
    parser.add_argument(
        "--output-basename",
        type=Path,
        default=OUTPUT_BASENAME,
        help="Output path without extension. PNG and PDF will be created.",
    )
    return parser.parse_args()


def apply_publication_style() -> None:
    rc = {
        "font.family": "DejaVu Serif",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "font.size": 12 + FONT_DELTA,
        "axes.labelsize": 14 + FONT_DELTA,
        "xtick.labelsize": 12 + FONT_DELTA,
        "ytick.labelsize": 12 + FONT_DELTA,
        "legend.fontsize": 12 + FONT_DELTA,
        "savefig.transparent": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    }
    if hasattr(sns, "set_theme"):
        sns.set_theme(style="whitegrid", context="paper", rc=rc)
    else:
        sns.set(style="whitegrid", context="paper", rc=rc)


def fmt_int(value: int) -> str:
    return f"{int(value):,}".replace(",", " ")


def normalize_segment(value: object) -> str:
    if pd.isna(value):
        return ""
    value = str(value).split(",")[0].strip()
    return value.split("*")[0]


def load_cluster_table(path: Path, method: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t").copy()
    if "species" in df.columns:
        df = df[df["species"].eq("HomoSapiens")].copy()
    for col in ["v.segm", "j.segm"]:
        if col in df.columns:
            df[col] = df[col].map(normalize_segment)
    df["method"] = method
    df["cdr3aa"] = df["cdr3aa"].fillna("").astype(str)
    df["cid"] = df["cid"].fillna("").astype(str)
    return df.drop_duplicates(subset=POINT_KEY + ["cid"]).reset_index(drop=True)


def cluster_sizes_from_table(df: pd.DataFrame, method: str) -> pd.DataFrame:
    sizes = (
        df.drop_duplicates(subset=POINT_KEY + ["cid"])
        .groupby(["gene", "antigen.epitope", "cid"])
        .size()
        .reset_index(name="cluster_size")
    )
    sizes["method"] = method
    return sizes


def per_epitope_counts(df: pd.DataFrame) -> pd.DataFrame:
    points = df.drop_duplicates(subset=POINT_KEY)
    return (
        points.groupby(["gene", "antigen.epitope"])
        .size()
        .reset_index(name="clustered_points")
    )


def build_coverage_table(tcrnet_df: pd.DataFrame, redcea_df: pd.DataFrame) -> pd.DataFrame:
    tcr_counts = per_epitope_counts(tcrnet_df).rename(
        columns={"clustered_points": "tcrnet_clustered_points"}
    )
    red_counts = per_epitope_counts(redcea_df).rename(
        columns={"clustered_points": "redcea_clustered_points"}
    )
    coverage = tcr_counts.merge(red_counts, on=["gene", "antigen.epitope"], how="inner")
    coverage["redcea_minus_tcrnet"] = (
        coverage["redcea_clustered_points"] - coverage["tcrnet_clustered_points"]
    )
    coverage["redcea_to_tcrnet_ratio"] = (
        coverage["redcea_clustered_points"] / coverage["tcrnet_clustered_points"]
    )
    return coverage


def compute_multi_stats(
    gene: str,
    vdjdb_slim: pd.DataFrame,
    redcea_for_multi: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    redcea_multi_counts = (
        redcea_for_multi[gene][["cdr3aa", "antigen.epitope"]]
        .drop_duplicates()
        .groupby("cdr3aa")
        .size()
    )
    gene_vdjdb = vdjdb_slim[vdjdb_slim["gene"].eq(gene)].copy()
    initial_multi = pd.DataFrame(
        gene_vdjdb[["cdr3", "antigen.epitope"]].drop_duplicates().groupby("cdr3").size()[
            lambda x: x > 1
        ]
    )
    out = (
        initial_multi.reset_index()
        .rename(columns={0: "vdjdb_count", "cdr3": "cdr3aa"})
        .merge(
            pd.DataFrame(redcea_multi_counts, columns=["redcea_count"]).reset_index(),
            how="left",
            on="cdr3aa",
        )
        .fillna(0)
        .astype({"redcea_count": int})
    )
    return out[out["cdr3aa"].apply(lambda x: x.startswith("C") and x.endswith("F"))].copy()


def load_ylq_csar(path: Path, method: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    sub = df[
        df["gene"].eq(CHAIN)
        & df["antigen.epitope"].eq(EPITOPE)
        & df["cdr3aa"].astype(str).str.startswith(PREFIX)
    ].copy()
    sub["method"] = method
    sub["cdr3_len"] = sub["cdr3aa"].str.len()
    sub["cluster_short"] = sub["cid"].astype(str).str.extract(r"([^.]+)$")[0]
    sub["cluster_label"] = "C" + sub["cluster_short"]
    sub["csz"] = pd.to_numeric(sub["csz"], errors="coerce")
    return sub.sort_values(["cluster_short", "cdr3aa"]).reset_index(drop=True)


def prettify_axis(ax: plt.Axes, *, grid_axis: str = "y") -> None:
    ax.grid(axis=grid_axis, color="#d9d9d9", linewidth=0.7, alpha=0.75)
    ax.set_axisbelow(True)


def add_panel_label(ax: plt.Axes, label: str, *, x: float = -0.02, y: float = 1.02) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=22 + FONT_DELTA,
        fontweight="normal",
        ha="left",
        va="bottom",
        color="#222222",
    )


def redcea_length_palette(redcea_df: pd.DataFrame) -> dict[int, tuple[float, float, float]]:
    length_order = sorted(redcea_df["cdr3_len"].unique())
    greens = sns.color_palette("Greens", n_colors=len(length_order) + 3)[2:]
    return dict(zip(length_order, greens))


def align_to_longest(seqs: list[str]) -> list[str]:
    seqs = list(seqs)
    if not seqs:
        return []
    if len({len(seq) for seq in seqs}) == 1:
        return seqs

    ref = max(seqs, key=lambda seq: (len(seq), seq))
    aligned = []
    for seq in seqs:
        aln_ref, aln_seq, *_ = pairwise2.align.globalms(
            ref,
            seq,
            2,
            -1,
            -5,
            -0.5,
            one_alignment_only=True,
        )[0]
        if "-" in aln_ref:
            aln_seq = "".join(s for r, s in zip(aln_ref, aln_seq) if r != "-")
        aligned.append(aln_seq)
    return aligned


def padded_alignment(aligned_seqs: list[str]) -> list[str]:
    max_len = max(len(seq) for seq in aligned_seqs)
    return [seq + "-" * (max_len - len(seq)) for seq in aligned_seqs]


def frequency_matrix(aligned_seqs: list[str]) -> pd.DataFrame:
    padded = padded_alignment(aligned_seqs)
    rows = []
    for pos in range(len(padded[0])):
        counts = {aa: 0 for aa in AA_ORDER}
        for seq in padded:
            aa = seq[pos]
            if aa in counts:
                counts[aa] += 1
        rows.append({aa: counts[aa] / len(padded) for aa in AA_ORDER})
    mat = pd.DataFrame(rows)
    mat.index = range(1, len(mat) + 1)
    return mat


def gap_fractions(aligned_seqs: list[str]) -> pd.Series:
    padded = padded_alignment(aligned_seqs)
    return pd.Series(
        [(sum(seq[pos] == "-" for seq in padded) / len(padded)) for pos in range(len(padded[0]))],
        index=range(1, len(padded[0]) + 1),
    )


def draw_gap_marks(ax: plt.Axes, gaps: pd.Series, *, y: float = 1.08) -> None:
    nonzero = gaps[gaps.gt(0)]
    if nonzero.empty:
        return
    for pos, frac in nonzero.items():
        ax.hlines(
            y,
            pos - 0.36,
            pos + 0.36,
            color="#1f1f1f",
            linewidth=2.5 + 3.0 * frac,
            clip_on=False,
        )
        ax.text(
            pos + 0.08,
            y + 0.13,
            f"{frac * 100:.0f}%",
            ha="center",
            va="bottom",
            rotation=45,
            rotation_mode="anchor",
            fontsize=10 + FONT_DELTA,
            color="#1f1f1f",
            clip_on=False,
        )


def plot_logo(
    ax: plt.Axes,
    seqs: list[str],
    *,
    show_gap_marks: bool = False,
    show_ylabel: bool = True,
    show_xlabel: bool = False,
) -> None:
    aligned = align_to_longest(seqs)
    mat = frequency_matrix(aligned)
    gaps = gap_fractions(aligned)
    logo = logomaker.Logo(
        mat,
        ax=ax,
        color_scheme="NajafabadiEtAl2017",
        width=0.90,
        vpad=0.04,
    )
    logo.style_spines(visible=False)
    logo.style_spines(spines=["left", "bottom"], visible=True, linewidth=0.8)
    ax.set_ylim(0, 1.22 if show_gap_marks and gaps.gt(0).any() else 1.02)
    ax.set_xlim(0.5, len(mat) + 0.5)
    ax.set_ylabel("Frequency" if show_ylabel else "")
    xticks = list(range(1, len(mat), 2))
    ax.set_xticks(xticks)
    ax.tick_params(axis="x", labelsize=10 + FONT_DELTA, labelrotation=20)
    ax.tick_params(axis="y", labelsize=10 + FONT_DELTA)
    ax.grid(axis="y", alpha=0.18, linewidth=0.7)
    if show_gap_marks:
        draw_gap_marks(ax, gaps)
    if show_xlabel:
        ax.set_xlabel("Aligned CDR3 position")
    else:
        ax.set_xlabel("")


def plot_counts_panel(ax: plt.Axes, chain_totals: pd.DataFrame) -> None:
    x = np.arange(2)
    width = 0.24
    count_methods = ["TCRnet", "RedCEA", "VDJdb"]
    for midx, method in enumerate(count_methods):
        vals = [
            int(
                chain_totals.loc[
                    (chain_totals["method"].eq(method)) & (chain_totals["gene"].eq(gene)),
                    "clustered_points",
                ].iloc[0]
            )
            for gene in ["TRA", "TRB"]
        ]
        bars = ax.bar(
            x + (midx - 1) * width,
            vals,
            width=width,
            color=COLORS[method],
            label=method,
        )
        for bar_idx, (bar, value) in enumerate(zip(bars, vals)):
            x_offset = (-0.03, 0.0, 0.03)[midx]
            y_mult = 1.02
            if bar_idx == 0:
                # The three TRA labels are packed tightly, so pin them to
                # deliberately separated positions instead of relying on rotation.
                x_offset = (-0.055, 0.030, 0.045)[midx]
                y_mult = (1.005, 1.135, 1.075)[midx]
            if bar_idx == 0 and midx == 1:
                y_mult = 0.9
            ax.text(
                bar.get_x() + bar.get_width() / 2 + x_offset,
                value * y_mult,
                fmt_int(value),
                ha="center",
                va="bottom",
                fontsize=10 + FONT_DELTA,
                rotation=15,
                rotation_mode="anchor",
            )
    ax.set_xticks(x)
    ax.set_xticklabels(["TRA", "TRB"])
    ax.set_ylabel("Rows / unique clustered points")
    ax.set_ylim(0, chain_totals["clustered_points"].max() * 1.20)
    ax.legend(frameon=False, ncol=1, loc="upper left")
    prettify_axis(ax)
    add_panel_label(ax, "A")


def plot_cluster_sizes_panel(ax: plt.Axes, cluster_size_df: pd.DataFrame) -> None:
    swarm_parts = []
    for (_gene, _method), sub in cluster_size_df.groupby(["gene", "method"]):
        sub = sub[sub["cluster_size"].gt(0)].copy()
        cutoff = sub["cluster_size"].quantile(0.95)
        sub = sub[sub["cluster_size"].le(cutoff)].copy()
        sub["log_cluster_size"] = np.log10(sub["cluster_size"])
        swarm_parts.append(sub)
    cluster_swarm_df = pd.concat(swarm_parts, ignore_index=True)

    sns.swarmplot(
        data=cluster_swarm_df,
        x="gene",
        y="log_cluster_size",
        hue="method",
        order=["TRA", "TRB"],
        hue_order=["TCRnet", "RedCEA"],
        dodge=True,
        size=1.8,
        alpha=0.58,
        palette={"TCRnet": COLORS["TCRnet"], "RedCEA": COLORS["RedCEA"]},
        ax=ax,
    )
    for xpos, gene in enumerate(["TRA", "TRB"]):
        for offset, method in [(-0.2, "TCRnet"), (0.2, "RedCEA")]:
            vals = cluster_swarm_df[
                (cluster_swarm_df["gene"].eq(gene)) & (cluster_swarm_df["method"].eq(method))
            ]["log_cluster_size"]
            if len(vals):
                ax.scatter(
                    xpos + offset,
                    vals.mean(),
                    s=32,
                    color="white",
                    edgecolor="black",
                    zorder=4,
                )
                ax.plot(
                    [xpos + offset - 0.09, xpos + offset + 0.09],
                    [vals.median(), vals.median()],
                    color="black",
                    linewidth=1.1,
                    zorder=4,
                )
    ax.set_yticks([0, np.log10(3), 1, np.log10(30), 2, np.log10(300)])
    ax.set_yticklabels(["1", "3", "10", "30", "100", "300"])
    ax.set_ylabel("Members per cluster, log scale")
    ax.set_xlabel("")
    ax.set_ylim(0, np.log10(300))
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], frameon=False, loc="upper left", ncol=2)
    prettify_axis(ax)
    add_panel_label(ax, "B")


def plot_multi_cdr3(ax: plt.Axes, stats: pd.DataFrame, label: str, seed: int) -> None:
    rng = np.random.default_rng(seed)
    plot_df = stats.copy()
    plot_df["vdjdb_jitter"] = plot_df["vdjdb_count"] + rng.uniform(-0.08, 0.08, len(plot_df))
    plot_df["redcea_jitter"] = plot_df["redcea_count"] + rng.uniform(-0.06, 0.06, len(plot_df))
    ax.scatter(
        plot_df["vdjdb_jitter"],
        plot_df["redcea_jitter"],
        alpha=0.35,
        s=35,
        color=COLORS["TCRnet"],
        edgecolors="none",
    )
    ax.set_xscale("symlog")
    ax.set_xlim(left=1.8)
    xticks = [2, 3, 5, 10, 20, 50, 70]
    xticks = [x for x in xticks if x <= stats["vdjdb_count"].max()]
    ax.set_xticks(xticks)
    tick_labels = [str(x) for x in xticks]
    if 50 in xticks and 70 in xticks:
        tick_labels[xticks.index(50)] = "50 "
        tick_labels[xticks.index(70)] = "   70"
    ax.set_xticklabels(tick_labels, rotation=20, ha="right", rotation_mode="anchor")
    y_counts = stats["redcea_count"].value_counts().sort_index()
    ax.set_yticks(y_counts.index)
    x_text = ax.get_xlim()[1]
    for y, count in y_counts.items():
        ax.text(
            x_text,
            y,
            f"  n={count}",
            va="center",
            ha="left",
            fontsize=10 + FONT_DELTA,
            color="dimgray",
        )
    ax.set_xlim(ax.get_xlim()[0], ax.get_xlim()[1] * 1.25)
    ax.set_xlabel("vdjdb_count")
    ax.set_ylabel("redcea_count")
    prettify_axis(ax, grid_axis="both")
    add_panel_label(ax, label)


def plot_umap(ax: plt.Axes, df: pd.DataFrame, method: str, panel_label: str) -> None:
    plot_df = df.copy()
    if method == "RedCEA":
        center = plot_df[["x", "y"]].median()
        dx = plot_df["x"] - center["x"]
        dy = plot_df["y"] - center["y"]
        dist = np.sqrt(dx**2 + dy**2)
        outlier_cutoff = dist.quantile(0.90)
        shrink = np.where(dist.gt(outlier_cutoff), 0.72, 1.0)
        rng = np.random.default_rng(7)
        plot_df["plot_x"] = center["x"] + dx * shrink + rng.normal(0, 0.006, len(plot_df))
        plot_df["plot_y"] = center["y"] + dy * shrink + rng.normal(0, 0.006, len(plot_df))
        if len(plot_df) > 1:
            sorted_y = plot_df["plot_y"].sort_values()
            lowest_idx = sorted_y.index[0]
            next_lowest_y = sorted_y.iloc[1]
            plot_df.loc[lowest_idx, "plot_y"] = next_lowest_y
    else:
        plot_df["plot_x"] = plot_df["x"]
        plot_df["plot_y"] = plot_df["y"]

    min_len = plot_df["cdr3_len"].min()
    max_len = plot_df["cdr3_len"].max()
    denom = max(max_len - min_len, 1)
    sizes = 42 + 95 * (plot_df["cdr3_len"] - min_len) / denom

    if method == "RedCEA":
        length_palette = redcea_length_palette(plot_df)
        for length, sub in plot_df.groupby("cdr3_len"):
            ax.scatter(
                sub["plot_x"],
                sub["plot_y"],
                s=sizes.loc[sub.index],
                color=length_palette[length],
                edgecolors="white",
                linewidths=0.65,
                alpha=0.90,
            )
    else:
        for cluster_label, sub in plot_df.groupby("cluster_label"):
            ax.scatter(
                sub["plot_x"],
                sub["plot_y"],
                s=sizes.loc[sub.index],
                color=TCRNET_CLUSTER_COLORS[cluster_label],
                edgecolors="white",
                linewidths=0.65,
                alpha=0.88,
            )

    for cluster_label, sub in plot_df.groupby("cluster_label"):
        cx = sub["plot_x"].median()
        cy = sub["plot_y"].median()
        n = len(sub)
        ax.annotate(
            f"{cluster_label}\nn={n}",
            xy=(cx, cy),
            xytext=(18, -18) if method == "RedCEA" else (8, 8),
            textcoords="offset points",
            fontsize=12 + FONT_DELTA,
            fontweight="bold",
            color="#222222",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.82),
        )

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    if method == "RedCEA":
        y_min = plot_df["plot_y"].min()
        y_max = plot_df["plot_y"].max()
        y_span = max(y_max - y_min, 1e-6)
        ax.set_ylim(y_min - 0.75 * y_span, y_max + 0.06 * y_span)
    else:
        ax.set_ylim(bottom=-1600, top=-600)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(True, alpha=0.18, linewidth=0.7)
    add_panel_label(ax, panel_label, x=-0.01, y=1.01)
    sns.despine(ax=ax)

    if method == "RedCEA":
        palette = redcea_length_palette(plot_df)
        handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="None",
                markersize=8,
                markerfacecolor=palette[length],
                markeredgecolor="white",
                label=str(length),
            )
            for length in sorted(palette)
        ]
        ax.legend(
            handles=handles,
            frameon=True,
            fancybox=True,
            title="CDR3 length",
            loc="upper center",
            bbox_to_anchor=(0.50, 0.98),
            ncol=4,
            borderpad=0.35,
            handletextpad=0.45,
            columnspacing=0.85,
        )
    else:
        ax.legend(
            handles=[
                Line2D([0], [0], marker="o", linestyle="None", markersize=9, color=color, label=label)
                for label, color in TCRNET_CLUSTER_COLORS.items()
            ],
            frameon=True,
            fancybox=True,
            title="Cluster",
            loc="upper right",
            borderpad=0.35,
            handletextpad=0.45,
        )


def add_tcrnet_logo_insets(ax: plt.Axes, tcrnet_df: pd.DataFrame) -> None:
    inset_specs = [
        ("C21", [0.035, 0.055, 0.42, 0.19]),
        ("C24", [0.520, 0.055, 0.42, 0.19]),
    ]
    for cluster_label, bounds in inset_specs:
        sub = tcrnet_df[tcrnet_df["cluster_label"].eq(cluster_label)]
        inset = ax.inset_axes(bounds)
        plot_logo(
            inset,
            sub["cdr3aa"].tolist(),
            show_gap_marks=False,
            show_ylabel=False,
            show_xlabel=False,
        )
        inset.set_xlabel("")
        inset.set_facecolor("white")
        inset.patch.set_alpha(0.86)
        for spine in inset.spines.values():
            spine.set_visible(False)
        inset.set_yticks([])
        inset.tick_params(axis="x", labelsize=11 + FONT_DELTA, length=2, pad=1)
        inset.tick_params(axis="y", length=0)


def add_redcea_summary_insets(ax: plt.Axes, redcea_df: pd.DataFrame) -> None:
    logo_inset = ax.inset_axes([0.365, 0.070, 0.62, 0.29])
    plot_logo(
        logo_inset,
        redcea_df["cdr3aa"].tolist(),
        show_gap_marks=True,
        show_ylabel=False,
        show_xlabel=False,
    )
    logo_inset.set_xlabel("")
    logo_inset.set_facecolor("white")
    logo_inset.patch.set_alpha(0.88)
    for spine in logo_inset.spines.values():
        spine.set_visible(False)
    logo_inset.set_yticks([])
    logo_inset.tick_params(axis="x", labelsize=11 + FONT_DELTA, length=2, pad=1)
    logo_inset.tick_params(axis="y", length=0)

    length_counts = redcea_df["cdr3_len"].value_counts().sort_index()
    palette = redcea_length_palette(redcea_df)
    dist_inset = ax.inset_axes([0.055, 0.105, 0.26, 0.20])
    x = np.arange(len(length_counts))
    dist_inset.bar(
        x,
        length_counts.values,
        color=[palette[length] for length in length_counts.index],
        edgecolor="white",
        linewidth=0.7,
        alpha=0.94,
    )
    dist_inset.set_xticks(x)
    dist_inset.set_xticklabels(length_counts.index.astype(str), rotation=20, ha="right", rotation_mode="anchor")
    for xpos, value in zip(x, length_counts.values):
        dist_inset.text(
            xpos,
            value + max(length_counts.values) * 0.035,
            str(int(value)),
            ha="center",
            va="bottom",
            fontsize=11 + FONT_DELTA,
        )
    dist_inset.set_xlabel("aa", fontsize=11 + FONT_DELTA, labelpad=1)
    dist_inset.set_ylim(0, max(length_counts.values) * 1.30)
    dist_inset.set_yticks([])
    dist_inset.tick_params(axis="x", labelsize=11 + FONT_DELTA, length=2, pad=1)
    dist_inset.tick_params(axis="y", length=0)
    dist_inset.grid(axis="y", alpha=0.22, linewidth=0.5)
    dist_inset.set_facecolor("white")
    dist_inset.patch.set_alpha(0.88)
    sns.despine(ax=dist_inset)


def build_figure(output_basename: Path) -> None:
    apply_publication_style()

    vdjdb_slim = pd.read_csv(VDJDB_SLIM_PATH, sep="\t")
    redcea_for_multi = {
        "TRA": pd.read_csv(REDCEA_TRA_PATH, sep="\t"),
        "TRB": pd.read_csv(REDCEA_TRB_PATH, sep="\t"),
    }

    multi_stats_tra = compute_multi_stats("TRA", vdjdb_slim, redcea_for_multi)
    multi_stats_trb = compute_multi_stats("TRB", vdjdb_slim, redcea_for_multi)

    tcrnet_all = load_cluster_table(TCRNET_CLUSTER_PATH, "TCRnet")
    redcea_all = pd.concat(
        [
            load_cluster_table(REDCEA_TRA_PATH, "RedCEA"),
            load_cluster_table(REDCEA_TRB_PATH, "RedCEA"),
        ],
        ignore_index=True,
    )
    coverage = build_coverage_table(tcrnet_all, redcea_all)
    coverage = coverage[coverage["gene"].isin(["TRA", "TRB"])].copy()
    cluster_size_df = pd.concat(
        [
            cluster_sizes_from_table(tcrnet_all, "TCRnet"),
            cluster_sizes_from_table(redcea_all, "RedCEA"),
        ],
        ignore_index=True,
    )

    chain_totals = []
    vdjdb_raw_counts = vdjdb_slim[vdjdb_slim["species"].eq("HomoSapiens")].groupby("gene").size()
    for gene in ["TRA", "TRB"]:
        chain_totals.append(
            {
                "method": "VDJdb",
                "gene": gene,
                "clustered_points": int(vdjdb_raw_counts.get(gene, 0)),
            }
        )
    for method, table in [("TCRnet", tcrnet_all), ("RedCEA", redcea_all)]:
        by_chain = table.drop_duplicates(subset=POINT_KEY).groupby("gene").size()
        for gene in ["TRA", "TRB"]:
            chain_totals.append(
                {
                    "method": method,
                    "gene": gene,
                    "clustered_points": int(by_chain.get(gene, 0)),
                }
            )
    chain_totals_df = pd.DataFrame(chain_totals)

    tcrnet_csar = load_ylq_csar(TCRNET_CLUSTER_PATH, "TCRnet")
    redcea_csar = load_ylq_csar(REDCEA_TRB_PATH, "RedCEA")

    fig = plt.figure(figsize=(17.5, 10.5))
    ax_a = fig.add_axes([0.06, 0.66, 0.17, 0.28])
    ax_b = fig.add_axes([0.29, 0.66, 0.24, 0.28])
    ax_c = fig.add_axes([0.59, 0.68, 0.145, 0.24])
    ax_d = fig.add_axes([0.825, 0.68, 0.145, 0.24])
    ax_e = fig.add_axes([0.06, 0.08, 0.42, 0.48])
    ax_f = fig.add_axes([0.56, 0.08, 0.40, 0.48])

    plot_counts_panel(ax_a, chain_totals_df)
    plot_cluster_sizes_panel(ax_b, cluster_size_df)
    plot_multi_cdr3(ax_c, multi_stats_tra, "C", 17)
    plot_multi_cdr3(ax_d, multi_stats_trb, "D", 19)
    plot_umap(ax_e, tcrnet_csar, "TCRnet", "E")
    add_tcrnet_logo_insets(ax_e, tcrnet_csar)
    plot_umap(ax_f, redcea_csar, "RedCEA", "F")
    add_redcea_summary_insets(ax_f, redcea_csar)

    output_basename.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_basename.with_suffix(".png"), dpi=300)
    fig.savefig(output_basename.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    args = parse_args()
    build_figure(args.output_basename)
    print("Saved figure to:")
    print(f"  {args.output_basename.with_suffix('.png')}")
    print(f"  {args.output_basename.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
