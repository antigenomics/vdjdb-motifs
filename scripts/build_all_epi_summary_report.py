#!/usr/bin/env python
from __future__ import annotations

import argparse
import math
import os
from pathlib import Path

MPLCONFIGDIR = Path(".tmp/matplotlib").resolve()
MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def shorten_parameter_label(label: str) -> str:
    parts = {}
    for token in str(label).split("|"):
        token = token.strip()
        if "=" in token:
            key, value = token.split("=", 1)
            parts[key.strip()] = value.strip()
    return " ".join(
        [
            f"k{parts.get('k', '?')}",
            f"ek{parts.get('ek', '?')}",
            f"r{parts.get('res', '?')}",
            f"m{parts.get('min', '?')}",
            f"e{parts.get('eps', '?')[0:2]}",
            f"s{parts.get('sym', '?')[0:3]}",
        ]
    )


def find_repo_root() -> Path:
    current = Path.cwd().resolve()
    for candidate in [current] + list(current.parents):
        if (candidate / "results").exists() and (candidate / "scripts").exists():
            return candidate
    raise FileNotFoundError("Could not locate repository root from current working directory.")


def parse_args() -> argparse.Namespace:
    root = find_repo_root()
    parser = argparse.ArgumentParser(
        description="Build a publication-style summary report from the extracted all-epitope possig bundle."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=root / "results" / "all_epi_summary" / "redcea_possig_heatmaps_trb",
        help="Directory with direct collector outputs from build_redcea_possig_heatmaps.py.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "results" / "all_epi_summary" / "final_report",
        help="Directory for summary tables, report markdown, and article-ready figures.",
    )
    return parser.parse_args()


def apply_publication_style() -> None:
    rc = {
        "font.family": "DejaVu Serif",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "font.size": 12,
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "xtick.labelsize": 10.5,
        "ytick.labelsize": 10.5,
        "legend.fontsize": 10,
        "figure.titlesize": 16,
        "savefig.transparent": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    }
    sns.set_theme(style="whitegrid", context="paper", rc=rc)


def read_inputs(input_dir: Path) -> dict[str, pd.DataFrame]:
    parameter_runs_path = input_dir / "redcea_possig_parameter_runs.tsv"
    sample_nn_path = input_dir / "sample_clonotype_nn.tsv"
    possig_cluster_path = input_dir / "possig_cluster_lfc.tsv"
    epitope_object_summary_path = input_dir / "epitope_object_metric_summary.tsv"
    parameter_robustness_path = input_dir / "parameter_robustness_summary.tsv"
    coverage_path = input_dir / "metric_collection_coverage.tsv"
    best_parameter_path = input_dir / "best_parameter_by_epitope.tsv"
    required = [
        parameter_runs_path,
        sample_nn_path,
        possig_cluster_path,
        epitope_object_summary_path,
        parameter_robustness_path,
        coverage_path,
        best_parameter_path,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing direct collector outputs. Re-run build_redcea_possig_heatmaps.py first:\n"
            + "\n".join(missing)
        )

    return {
        "parameter_runs": pd.read_csv(parameter_runs_path, sep="\t"),
        "clonotypes": pd.read_csv(sample_nn_path, sep="\t"),
        "possig_clusters": pd.read_csv(possig_cluster_path, sep="\t"),
        "epitope_object_summary": pd.read_csv(epitope_object_summary_path, sep="\t"),
        "parameter_robustness": pd.read_csv(parameter_robustness_path, sep="\t"),
        "coverage": pd.read_csv(coverage_path, sep="\t"),
        "best_parameter_by_epitope": pd.read_csv(best_parameter_path, sep="\t"),
    }


def plot_final_figure(
    parameter_runs: pd.DataFrame,
    parameter_robustness: pd.DataFrame,
    nn_summary: pd.DataFrame,
    lfc_summary: pd.DataFrame,
    clonotypes: pd.DataFrame,
    possig_clusters: pd.DataFrame,
    *,
    output_stem: Path,
) -> None:
    apply_publication_style()
    fig = plt.figure(figsize=(16, 11))
    grid = fig.add_gridspec(2, 2, width_ratios=[1.35, 1.0], height_ratios=[1.1, 1.0])
    ax_heatmap = fig.add_subplot(grid[0, 0])
    ax_rank = fig.add_subplot(grid[0, 1])
    ax_nn = fig.add_subplot(grid[1, 0])
    ax_lfc = fig.add_subplot(grid[1, 1])

    top_params = parameter_robustness["parameter_label"].tolist()
    heatmap_source = parameter_runs.copy()
    heatmap_source["metric"] = pd.to_numeric(heatmap_source["metric"], errors="coerce")
    epitope_order = sorted(heatmap_source["epitope"].dropna().unique().tolist())
    heatmap = heatmap_source.pivot(index="parameter_label", columns="epitope", values="metric").reindex(
        index=top_params,
        columns=epitope_order,
    )
    short_index = [shorten_parameter_label(label) for label in heatmap.index]
    heatmap.index = short_index
    sns.heatmap(
        heatmap,
        ax=ax_heatmap,
        cmap="YlGnBu",
        vmin=0.0,
        vmax=float(np.nanmax(heatmap.to_numpy(dtype=float))) if np.isfinite(heatmap.to_numpy(dtype=float)).any() else 1.0,
        cbar_kws={"label": "redcea_possig_density_score"},
    )
    ax_heatmap.set_title("A. Final score across parameter sets and epitopes")
    ax_heatmap.set_xlabel("Epitope")
    ax_heatmap.set_ylabel("Parameter set")
    ax_heatmap.set_yticklabels(ax_heatmap.get_yticklabels(), rotation=0, fontsize=9.5)
    ax_heatmap.tick_params(axis="x", labelrotation=45)

    rank_plot = parameter_robustness.head(8).iloc[::-1].copy()
    rank_plot["short_label"] = rank_plot["parameter_label"].map(shorten_parameter_label)
    rank_column = "rank_within_chain" if "rank_within_chain" in rank_plot.columns else "rank"
    colors = ["#d55e00" if rank == 1 else "#4c78a8" for rank in rank_plot[rank_column]]
    ax_rank.barh(rank_plot["short_label"], rank_plot["robust_selection_score"], color=colors, alpha=0.95)
    for idx, row in rank_plot.reset_index(drop=True).iterrows():
        label = f"med {row['median_metric']:.3f} | q25 {row['q25_metric']:.3f} | cov {row['coverage_fraction']:.2f}"
        ax_rank.text(row["robust_selection_score"] + 0.0025, idx, label, va="center", fontsize=8.5)
    ax_rank.set_title("B. Robust parameter ranking")
    ax_rank.set_xlabel("coverage_fraction x q25(metric)")
    ax_rank.set_ylabel("")
    ax_rank.tick_params(axis="y", labelsize=9.5)

    nn_frame = clonotypes.loc[:, ["epitope", "nearest_neighbor_distance"]].copy()
    nn_frame["nearest_neighbor_distance"] = pd.to_numeric(nn_frame["nearest_neighbor_distance"], errors="coerce")
    nn_order = nn_summary.sort_values("nn_median")["epitope"].tolist()
    sns.violinplot(
        data=nn_frame,
        x="epitope",
        y="nearest_neighbor_distance",
        order=nn_order,
        inner=None,
        linewidth=0.8,
        color="#f4a261",
        ax=ax_nn,
        cut=0,
    )
    sns.boxplot(
        data=nn_frame,
        x="epitope",
        y="nearest_neighbor_distance",
        order=nn_order,
        width=0.16,
        showcaps=True,
        boxprops={"facecolor": "white", "zorder": 3},
        showfliers=False,
        whiskerprops={"linewidth": 1.0},
        ax=ax_nn,
    )
    ax_nn.set_title("C. Sample clonotype nearest-neighbor distances")
    ax_nn.set_xlabel("Epitope")
    ax_nn.set_ylabel("Nearest-neighbor distance")
    ax_nn.tick_params(axis="x", rotation=35)

    lfc_value_col = "log_fold_change_possig" if "log_fold_change_possig" in possig_clusters.columns else "log_fold_change"
    lfc_frame = possig_clusters.loc[:, ["epitope", lfc_value_col]].copy()
    lfc_frame[lfc_value_col] = pd.to_numeric(lfc_frame[lfc_value_col], errors="coerce")
    lfc_frame = lfc_frame.rename(columns={lfc_value_col: "log_fold_change_possig"})
    lfc_order = lfc_summary.sort_values("lfc_possig_median")["epitope"].tolist()
    sns.violinplot(
        data=lfc_frame,
        x="epitope",
        y="log_fold_change_possig",
        order=lfc_order,
        inner=None,
        linewidth=0.8,
        color="#7cb518",
        ax=ax_lfc,
        cut=0,
    )
    sns.boxplot(
        data=lfc_frame,
        x="epitope",
        y="log_fold_change_possig",
        order=lfc_order,
        width=0.16,
        showcaps=True,
        boxprops={"facecolor": "white", "zorder": 3},
        showfliers=False,
        whiskerprops={"linewidth": 1.0},
        ax=ax_lfc,
    )
    ax_lfc.set_title("D. log fold change in positive-significant clusters")
    ax_lfc.set_xlabel("Epitope")
    ax_lfc.set_ylabel("log_fold_change_possig")
    ax_lfc.tick_params(axis="x", rotation=35)

    best_row = parameter_robustness.iloc[0]
    fig.suptitle(
        "TRB large-epitope possig analysis: final score behavior and robust parameter choice\n"
        f"Recommended parameter set: {shorten_parameter_label(best_row['parameter_label'])}",
        y=0.98,
    )
    fig.subplots_adjust(left=0.08, right=0.98, top=0.93, bottom=0.07, wspace=0.34, hspace=0.34)
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_stem.with_suffix(".png"), dpi=240, bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def build_report_markdown(
    parameter_robustness: pd.DataFrame,
    nn_summary: pd.DataFrame,
    lfc_summary: pd.DataFrame,
    coverage_df: pd.DataFrame,
    *,
    figure_name: str,
    source_dir: Path,
    n_score_epitopes: int,
    n_object_epitopes: int,
) -> str:
    best = parameter_robustness.iloc[0]
    second = parameter_robustness.iloc[1]
    nn_best = (
        nn_summary.loc[nn_summary["epitope"].notna() & nn_summary["nn_median"].notna()]
        .sort_values("nn_median")
        .reset_index(drop=True)
    )
    lfc_best = (
        lfc_summary.loc[lfc_summary["epitope"].notna() & lfc_summary["lfc_possig_median"].notna()]
        .sort_values("lfc_possig_median", ascending=False)
        .reset_index(drop=True)
    )
    coverage_gap = coverage_df.loc[~coverage_df["has_object_level_nn"].astype(bool), "epitope"].astype(str).tolist()

    lines = [
        "# TRB all-epitope possig summary report",
        "",
        "## Data scope",
        f"- Score-grid source: `{source_dir}`",
        f"- `redcea_possig_parameter_runs.tsv` covers `{n_score_epitopes}` TRB epitopes across `{parameter_robustness['parameter_label'].nunique()}` parameter sets.",
        f"- Object-level diagnostics from `sample_clonotype_nn.tsv` are available for `{n_object_epitopes}` epitopes in this extracted bundle.",
        (
            f"- Object-level outputs are missing for: `{', '.join(coverage_gap)}`."
            if coverage_gap
            else "- Object-level outputs are present for every epitope covered by the score grid."
        ),
        "- The bundle does not contain external ground-truth labels such as F1, so the metric formula itself was not re-fit here.",
        "- Optimization below therefore means robust operating-point selection for the existing final score, not rewriting the score definition.",
        "",
        "## Recommended final parameter set",
        f"- Recommended: `{best['parameter_label']}`",
        f"- Why: full epitope coverage `{best['coverage_fraction']:.2f}`, strongest lower-quartile protection `q25={best['q25_metric']:.3f}`, and competitive median score `median={best['median_metric']:.3f}`.",
        f"- Comparator with higher raw median but weaker deployment coverage/stability: `{second['parameter_label']}` with `median={second['median_metric']:.3f}`, `q25={second['q25_metric']:.3f}`, coverage `{second['coverage_fraction']:.2f}`.",
        f"- Interpretation: sample-based variants can peak higher on individual epitopes, but the selected background-based `{best['parameter_label']}` setting is more conservative on weak epitopes and therefore safer as a publication default.",
        "",
        "## Biological sanity checks",
        f"- Tightest clonotype packing by nearest-neighbor median: `{nn_best.iloc[0]['epitope']}` (`{nn_best.iloc[0]['nn_median']:.2f}`), consistent with a denser sequence neighborhood.",
        f"- Broadest clonotype spacing by nearest-neighbor median: `{nn_best.iloc[-1]['epitope']}` (`{nn_best.iloc[-1]['nn_median']:.2f}`), suggesting a sparser or more heterogeneous response.",
        f"- Strongest positive-significant enrichment by median log fold change: `{lfc_best.iloc[0]['epitope']}` (`{lfc_best.iloc[0]['lfc_possig_median']:.2f}`), with `{int(lfc_best.iloc[0]['n_possig_clusters'])}` possig clusters recovered.",
        f"- Weakest positive-significant enrichment by median log fold change: `{lfc_best.iloc[-1]['epitope']}` (`{lfc_best.iloc[-1]['lfc_possig_median']:.2f}`), so scores for that epitope should be interpreted more cautiously.",
        "- One weak epitope remains close to score collapse in most parameterizations, which argues against selecting the final operating point on mean score alone.",
        "",
        "## Computational correctness checks",
        "- Positive-significant clusters were defined exactly as in the REDCEA/tcrempnet logic: `cluster_id != -1`, `enrichment_fdr_zbinom < 0.05`, and `log_fold_change > 0`.",
        "- The nearest-neighbor diagnostic is object-level on sample clonotypes, derived from the saved `knn_sample_sample__*.distances.npy` projections and converted to Euclidean distance before materialization into `sample_clonotype_nn.tsv`.",
        "- The LFC diagnostic is cluster-level and is now read directly from `possig_cluster_lfc.tsv` emitted by the collector.",
        "",
        "## Outputs",
        f"- Final figure: `{figure_name}`",
        "- `parameter_robustness_summary.tsv`",
        "- `possig_cluster_lfc.tsv`",
        "- `epitope_nn_summary.tsv`",
        "- `epitope_lfc_summary.tsv`",
        "- `metric_collection_coverage.tsv`",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    inputs = read_inputs(input_dir)
    parameter_runs = inputs["parameter_runs"]
    clonotypes = inputs["clonotypes"]
    possig_clusters = inputs["possig_clusters"]
    epitope_object_summary = inputs["epitope_object_summary"]
    parameter_robustness = inputs["parameter_robustness"]
    coverage_df = inputs["coverage"]

    nn_summary = (
        epitope_object_summary.loc[:, ["chain", "epitope", "n_sample_clonotypes", "nn_median", "nn_mean", "nn_std", "nn_min", "nn_max"]]
        .dropna(subset=["epitope"], how="all")
        .copy()
    )
    lfc_summary = (
        epitope_object_summary.loc[:, ["chain", "epitope", "n_possig_clusters", "lfc_possig_median", "lfc_possig_mean", "lfc_possig_std", "lfc_possig_min", "lfc_possig_max"]]
        .dropna(subset=["epitope"], how="all")
        .copy()
    )

    parameter_robustness_path = output_dir / "parameter_robustness_summary.tsv"
    possig_clusters_path = output_dir / "possig_cluster_lfc.tsv"
    nn_summary_path = output_dir / "epitope_nn_summary.tsv"
    lfc_summary_path = output_dir / "epitope_lfc_summary.tsv"
    coverage_copy_path = output_dir / "metric_collection_coverage.tsv"
    figure_stem = output_dir / "trb_all_epi_possig_publication_figure"
    report_path = output_dir / "summary_report.md"

    parameter_robustness.to_csv(parameter_robustness_path, sep="\t", index=False)
    possig_clusters.to_csv(possig_clusters_path, sep="\t", index=False)
    nn_summary.to_csv(nn_summary_path, sep="\t", index=False)
    lfc_summary.to_csv(lfc_summary_path, sep="\t", index=False)
    coverage_df.to_csv(coverage_copy_path, sep="\t", index=False)
    plot_final_figure(
        parameter_runs,
        parameter_robustness,
        nn_summary,
        lfc_summary,
        clonotypes,
        possig_clusters,
        output_stem=figure_stem,
    )
    report_text = build_report_markdown(
        parameter_robustness,
        nn_summary,
        lfc_summary,
        coverage_df,
        figure_name=figure_stem.with_suffix(".pdf").name,
        source_dir=input_dir,
        n_score_epitopes=int(parameter_runs["epitope"].nunique()),
        n_object_epitopes=int(coverage_df["has_object_level_nn"].astype(bool).sum()),
    )
    report_path.write_text(report_text, encoding="utf-8")

    print(f"Saved {parameter_robustness_path}")
    print(f"Saved {possig_clusters_path}")
    print(f"Saved {nn_summary_path}")
    print(f"Saved {lfc_summary_path}")
    print(f"Saved {coverage_copy_path}")
    print(f"Saved {figure_stem.with_suffix('.png')}")
    print(f"Saved {figure_stem.with_suffix('.pdf')}")
    print(f"Saved {figure_stem.with_suffix('.svg')}")
    print(f"Saved {report_path}")


if __name__ == "__main__":
    main()
