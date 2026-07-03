#!/usr/bin/env python
from __future__ import annotations

import argparse
import math
import os
import re
from pathlib import Path

MPLCONFIGDIR = Path(".tmp/matplotlib").resolve()
MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SUMMARY_SUFFIX = "_summary_tcrempnet.tsv"
COMBO_RE = re.compile(
    r"vdbscan_leiden"
    r"_k(?P<k_neighbors>\d+)"
    r"_ek(?P<eps_k_neighbors>\d+)"
    r"_res(?P<resolution_tag>[A-Za-z0-9p.-]+)"
    r"_min(?P<cluster_min_samples>\d+)"
    r"_eps(?P<eps_estimation_based_on>[A-Za-z0-9_.-]+)"
    r"_sym(?P<vdbscan_sym_rule>[A-Za-z0-9_.-]+)$"
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
        description=(
            "Compute repo-native redcea_possig_density_score from per-run summary tables "
            "and build publication-style parameter heatmaps."
        )
    )
    parser.add_argument(
        "--runs-root",
        type=Path,
        required=True,
        help="Root directory with run folders, typically <...>/TRB/<EPITOPE>/<COMBO>/run_metadata.tsv.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "results" / "redcea_possig_heatmaps",
        help="Directory for TSV outputs and heatmap figures.",
    )
    parser.add_argument(
        "--glc-knn",
        type=Path,
        default=None,
        help="Optional reference GLC sample-sample kNN distance array used in d_ref.",
    )
    parser.add_argument(
        "--ylq-knn",
        type=Path,
        default=None,
        help="Optional reference YLQ sample-sample kNN distance array used in d_ref.",
    )
    parser.add_argument(
        "--d-ref-mode",
        choices=["dataset_median", "dataset_mean", "glc_ylq"],
        default="dataset_median",
        help=(
            "How to compute d_ref: median or mean of current-run d_epi values, "
            "or the original GLC/YLQ reference mode."
        ),
    )
    parser.add_argument(
        "--comparison-group",
        choices=["epitope", "chain_epitope"],
        default="chain_epitope",
        help="Grouping used for coverage_rank_pct and effective_size_penalty thresholds.",
    )
    return parser.parse_args()


def read_metadata_tsv(path: Path) -> dict[str, str]:
    frame = pd.read_csv(path, sep="\t")
    if not {"key", "value"}.issubset(frame.columns):
        raise ValueError(f"Expected key/value TSV in {path}")
    return {
        str(key): "" if pd.isna(value) else str(value)
        for key, value in zip(frame["key"], frame["value"], strict=False)
    }


def parse_resolution(tag: str) -> float:
    return float(str(tag).replace("p", "."))


def parse_combo_tag(combo_tag: str) -> dict[str, object]:
    match = COMBO_RE.match(combo_tag.strip())
    if not match:
        raise ValueError(f"Could not parse combo tag: {combo_tag}")
    values = match.groupdict()
    return {
        "algorithm": "vdbscan_leiden",
        "param_k_neighbors": int(values["k_neighbors"]),
        "param_eps_k_neighbors": int(values["eps_k_neighbors"]),
        "param_leiden_resolution": parse_resolution(values["resolution_tag"]),
        "param_resolution_tag": values["resolution_tag"],
        "param_cluster_min_samples": int(values["cluster_min_samples"]),
        "param_eps_estimation_based_on": values["eps_estimation_based_on"],
        "param_vdbscan_sym_rule": values["vdbscan_sym_rule"],
    }


def parameter_label_for_row(row: pd.Series) -> str:
    return " | ".join(
        [
            f"k={int(row['param_k_neighbors'])}",
            f"ek={int(row['param_eps_k_neighbors'])}",
            f"res={row['param_leiden_resolution']:g}",
            f"min={int(row['param_cluster_min_samples'])}",
            f"eps={row['param_eps_estimation_based_on']}",
            f"sym={row['param_vdbscan_sym_rule']}",
        ]
    )


def locate_single_file(run_dir: Path, pattern: str) -> Path:
    matches = sorted(run_dir.rglob(pattern))
    if not matches:
        raise FileNotFoundError(f"No files matching {pattern} under {run_dir}")
    return matches[0]


def load_reference_distance(path: Path) -> float:
    if not path.exists():
        raise FileNotFoundError(f"Reference distance file not found: {path}")
    return float(np.median(np.load(path)))


def resolve_reference_distance(
    run_df: pd.DataFrame,
    *,
    d_ref_mode: str,
    glc_knn: Path | None,
    ylq_knn: Path | None,
) -> tuple[float, str]:
    d_epi_values = pd.to_numeric(run_df["d_epi"], errors="coerce").dropna()
    if d_epi_values.empty:
        raise RuntimeError("Could not compute d_ref because no valid d_epi values were collected.")

    if d_ref_mode == "glc_ylq":
        if glc_knn is None or ylq_knn is None:
            raise ValueError("--glc-knn and --ylq-knn are required when --d-ref-mode glc_ylq is used.")
        d_glc = load_reference_distance(glc_knn)
        d_ylq = load_reference_distance(ylq_knn)
        d_ref = (d_glc + d_ylq) / 2.0
        return d_ref, f"glc_ylq (d_glc={d_glc:.6f}, d_ylq={d_ylq:.6f})"

    if d_ref_mode == "dataset_mean":
        return float(d_epi_values.mean()), "dataset_mean(d_epi)"

    return float(d_epi_values.median()), "dataset_median(d_epi)"


def possig_mask(summary_df: pd.DataFrame) -> pd.Series:
    return (
        (pd.to_numeric(summary_df["cluster_id"], errors="coerce") != -1)
        & (pd.to_numeric(summary_df["enrichment_fdr_zbinom"], errors="coerce") < 0.05)
        & (pd.to_numeric(summary_df["log_fold_change"], errors="coerce") > 0)
    )


def comparison_key(row: pd.Series, mode: str) -> tuple[str, ...]:
    if mode == "epitope":
        return (str(row["epitope"]),)
    return (str(row["chain"]), str(row["epitope"]))


def collect_run_rows(runs_root: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    metadata_paths = sorted(runs_root.rglob("run_metadata.tsv"))
    if not metadata_paths:
        raise FileNotFoundError(f"No run_metadata.tsv files found under {runs_root}")

    for metadata_path in metadata_paths:
        run_dir = metadata_path.parent
        metadata = read_metadata_tsv(metadata_path)
        chain = metadata.get("chain") or run_dir.parent.parent.name
        epitope = metadata.get("epitope") or run_dir.parent.name
        combo_tag = metadata.get("output_tag") or run_dir.name
        params = parse_combo_tag(combo_tag)

        summary_path = locate_single_file(run_dir, f"*{SUMMARY_SUFFIX}")
        knn_path = locate_single_file(run_dir, "knn_sample_sample__*.distances.npy")
        summary_df = pd.read_csv(summary_path, sep="\t")
        mask = possig_mask(summary_df)

        sample_values = pd.to_numeric(summary_df.get("sample", pd.Series(dtype=float)), errors="coerce")
        cluster_size_values = pd.to_numeric(summary_df.get("cluster_size", pd.Series(dtype=float)), errors="coerce")
        lfc_values = pd.to_numeric(summary_df.get("log_fold_change", pd.Series(dtype=float)), errors="coerce")

        sample_mass_enriched = float(sample_values.loc[mask].fillna(0.0).sum()) if len(summary_df) else 0.0
        cluster_mass_enriched = float(cluster_size_values.loc[mask].fillna(0.0).sum()) if len(summary_df) else 0.0
        background_mass_enriched = (
            float(cluster_mass_enriched - sample_mass_enriched)
            if pd.notna(cluster_mass_enriched) and pd.notna(sample_mass_enriched)
            else np.nan
        )
        total_sample = float(sample_values.fillna(0.0).sum()) if len(summary_df) else np.nan
        sample_coverage_enriched = (
            float(sample_mass_enriched / total_sample) if pd.notna(total_sample) and total_sample > 0 else np.nan
        )

        effective_n = np.nan
        effective_sample_cluster_size = np.nan
        if sample_mass_enriched > 0:
            cluster_weights = sample_values.loc[mask].fillna(0.0)
            weight_fraction = cluster_weights / sample_mass_enriched
            simpson_denom = float((weight_fraction**2).sum())
            if simpson_denom > 0:
                effective_n = float(1.0 / simpson_denom)
                effective_sample_cluster_size = float(sample_mass_enriched / effective_n)

        lfc_possig_mean = float(lfc_values.loc[mask].mean()) if mask.any() else np.nan
        d_epi = float(np.median(np.load(knn_path)))

        row = {
            "run_id": str(run_dir.relative_to(runs_root)).replace("\\", "/"),
            "run_dir": str(run_dir),
            "summary_path": str(summary_path),
            "knn_sample_sample_path": str(knn_path),
            "chain": str(chain),
            "epitope": str(epitope),
            "output_tag": combo_tag,
            "parameter_label": "",
            "n_clusters_in_summary": int(len(summary_df)),
            "n_clusters_enriched": int(mask.sum()),
            "sample_mass_enriched": sample_mass_enriched,
            "cluster_mass_enriched": cluster_mass_enriched,
            "background_mass_enriched": background_mass_enriched,
            "sample_coverage_enriched": sample_coverage_enriched,
            "effective_n_enriched_clusters_sample_weighted": effective_n,
            "effective_sample_cluster_size": effective_sample_cluster_size,
            "log_fold_change_possig__mean": lfc_possig_mean,
            "d_epi": d_epi,
            **params,
        }
        row["parameter_label"] = parameter_label_for_row(pd.Series(row))
        rows.append(row)

    frame = pd.DataFrame(rows)
    if frame.empty:
        raise RuntimeError(f"No run rows collected from {runs_root}")
    return frame


def append_repo_native_scores(run_df: pd.DataFrame, *, comparison_group: str, d_ref: float) -> pd.DataFrame:
    out = run_df.copy()
    out["comparison_group"] = out.apply(lambda row: " | ".join(comparison_key(row, comparison_group)), axis=1)
    out["coverage_rank_pct"] = np.nan
    out["eff_low_threshold"] = np.nan
    out["eff_high_threshold"] = np.nan
    out["effective_size_penalty"] = np.nan

    for group_value, group_frame in out.groupby("comparison_group", sort=True):
        group_index = group_frame.index
        coverage_source = np.log1p(pd.to_numeric(group_frame["sample_mass_enriched"], errors="coerce"))
        out.loc[group_index, "coverage_rank_pct"] = coverage_source.rank(method="average", pct=True)

        effective_sizes = pd.to_numeric(group_frame["effective_sample_cluster_size"], errors="coerce")
        if effective_sizes.notna().any():
            eff_low = max(float(effective_sizes.quantile(0.25)), 2.0)
            eff_high = float(effective_sizes.quantile(0.75))
        else:
            eff_low = np.nan
            eff_high = np.nan

        out.loc[group_index, "eff_low_threshold"] = eff_low
        out.loc[group_index, "eff_high_threshold"] = eff_high

        penalty = pd.Series(np.nan, index=group_index, dtype="float64")
        valid_eff_mask = effective_sizes.notna() & (effective_sizes > 0) & pd.notna(eff_low) & pd.notna(eff_high)
        if valid_eff_mask.any():
            valid_eff = effective_sizes.loc[valid_eff_mask]
            penalty.loc[valid_eff.index] = np.minimum(1.0, valid_eff / eff_low) * np.minimum(1.0, eff_high / valid_eff)
        out.loc[group_index, "effective_size_penalty"] = penalty

    zero_enriched_mask = out["sample_mass_enriched"].eq(0)
    out.loc[zero_enriched_mask, "effective_size_penalty"] = 0.0
    out["redcea_dense_score_base"] = out["coverage_rank_pct"] * out["effective_size_penalty"]

    out["final_gamma"] = np.power(d_ref / pd.to_numeric(out["d_epi"], errors="coerce"), 2.0)
    lfc_numeric = pd.to_numeric(out["log_fold_change_possig__mean"], errors="coerce")
    out["final_lfc_component"] = np.power(1.0 - np.exp(-lfc_numeric / 8.0), 0.5)
    out["redcea_possig_density_score"] = out["redcea_dense_score_base"] * np.power(
        out["final_lfc_component"],
        out["final_gamma"],
    )
    out.loc[zero_enriched_mask, "final_lfc_component"] = 0.0
    out.loc[zero_enriched_mask, "redcea_possig_density_score"] = 0.0
    return out


def summarize_parameter_heatmap(run_df: pd.DataFrame) -> pd.DataFrame:
    group_cols = [
        "chain",
        "parameter_label",
        "epitope",
        "param_cluster_min_samples",
        "param_eps_k_neighbors",
        "param_k_neighbors",
        "param_leiden_resolution",
        "param_eps_estimation_based_on",
        "param_vdbscan_sym_rule",
    ]
    summary = (
        run_df.groupby(group_cols, dropna=False)
        .agg(
            metric=("redcea_possig_density_score", "median"),
            count=("run_id", "size"),
        )
        .reset_index()
    )
    sort_cols = [
        "param_k_neighbors",
        "param_eps_k_neighbors",
        "param_leiden_resolution",
        "param_cluster_min_samples",
        "param_eps_estimation_based_on",
        "param_vdbscan_sym_rule",
        "parameter_label",
        "epitope",
    ]
    return summary.sort_values(sort_cols).reset_index(drop=True)


def summarize_distances(run_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    run_distance_df = (
        run_df.loc[:, ["chain", "epitope", "run_id", "parameter_label", "d_epi"]]
        .sort_values(["chain", "epitope", "parameter_label", "run_id"])
        .reset_index(drop=True)
    )
    epitope_distance_df = (
        run_df.groupby(["chain", "epitope"], dropna=False)
        .agg(
            d_epi_mean=("d_epi", "mean"),
            d_epi_median=("d_epi", "median"),
            d_epi_std=("d_epi", "std"),
            d_epi_min=("d_epi", "min"),
            d_epi_max=("d_epi", "max"),
            n_runs=("run_id", "size"),
        )
        .reset_index()
        .sort_values(["chain", "epitope"])
        .reset_index(drop=True)
    )
    return run_distance_df, epitope_distance_df


def plot_chain_heatmap(summary_df: pd.DataFrame, *, chain: str, output_stem: Path) -> None:
    chain_df = summary_df.loc[summary_df["chain"].astype(str) == str(chain)].copy()
    if chain_df.empty:
        return

    row_order = (
        chain_df.loc[
            :,
            [
                "parameter_label",
                "param_k_neighbors",
                "param_eps_k_neighbors",
                "param_leiden_resolution",
                "param_cluster_min_samples",
                "param_eps_estimation_based_on",
                "param_vdbscan_sym_rule",
            ],
        ]
        .drop_duplicates()
        .sort_values(
            [
                "param_k_neighbors",
                "param_eps_k_neighbors",
                "param_leiden_resolution",
                "param_cluster_min_samples",
                "param_eps_estimation_based_on",
                "param_vdbscan_sym_rule",
                "parameter_label",
            ]
        )["parameter_label"]
        .tolist()
    )
    column_order = sorted(chain_df["epitope"].astype(str).unique().tolist())
    heatmap_df = chain_df.pivot(index="parameter_label", columns="epitope", values="metric").reindex(
        index=row_order,
        columns=column_order,
    )

    values = heatmap_df.to_numpy(dtype=float)
    finite_values = values[np.isfinite(values)]
    vmin = float(finite_values.min()) if finite_values.size else 0.0
    vmax = float(finite_values.max()) if finite_values.size else 1.0
    if math.isclose(vmin, vmax):
        vmax = vmin + 1e-6

    fig_width = max(10.0, 0.8 * len(column_order) + 6.0)
    fig_height = max(8.0, 0.42 * len(row_order) + 2.5)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)
    image = ax.imshow(values, cmap="YlGnBu", vmin=vmin, vmax=vmax, aspect="auto")

    ax.set_xticks(range(len(column_order)))
    ax.set_xticklabels(column_order, rotation=45, ha="right")
    ax.set_yticks(range(len(row_order)))
    ax.set_yticklabels(row_order)
    ax.set_xlabel("Epitope")
    ax.set_ylabel("Parameter set")
    ax.set_title(f"redcea_possig_density_score heatmap ({chain})")

    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            value = values[i, j]
            if not np.isfinite(value):
                continue
            text_color = "white" if value >= (vmin + vmax) / 2.0 else "black"
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=8, color=text_color)

    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label("redcea_possig_density_score")

    output_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_stem.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def plot_epitope_distance_histograms(run_df: pd.DataFrame, *, output_dir: Path) -> None:
    distance_dir = output_dir / "distance_histograms"
    distance_dir.mkdir(parents=True, exist_ok=True)

    for (chain, epitope), group_df in run_df.groupby(["chain", "epitope"], sort=True):
        values = pd.to_numeric(group_df["d_epi"], errors="coerce").dropna().to_numpy(dtype=float)
        if values.size == 0:
            continue

        bins = min(12, max(5, int(np.ceil(np.sqrt(values.size)))))
        fig, ax = plt.subplots(figsize=(6.5, 4.5), constrained_layout=True)
        ax.hist(values, bins=bins, color="#4c78a8", edgecolor="white")
        ax.axvline(values.mean(), color="#f58518", linestyle="--", linewidth=2, label=f"mean={values.mean():.3f}")
        ax.axvline(np.median(values), color="#54a24b", linestyle="-.", linewidth=2, label=f"median={np.median(values):.3f}")
        ax.set_xlabel("d_epi")
        ax.set_ylabel("Run count")
        ax.set_title(f"{chain} {epitope}: distribution of run-level d_epi")
        ax.legend(frameon=False)

        stem = distance_dir / f"{str(chain).lower()}_{epitope}_d_epi_hist"
        fig.savefig(stem.with_suffix(".png"), dpi=220, bbox_inches="tight")
        fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
        plt.close(fig)


def plot_epitope_distance_summary(epitope_distance_df: pd.DataFrame, *, output_dir: Path) -> None:
    if epitope_distance_df.empty:
        return

    fig_height = max(5.0, 0.38 * len(epitope_distance_df) + 1.5)
    fig, ax = plt.subplots(figsize=(9.0, fig_height), constrained_layout=True)
    labels = [f"{row.chain}:{row.epitope}" for row in epitope_distance_df.itertuples(index=False)]
    y = np.arange(len(epitope_distance_df))
    means = epitope_distance_df["d_epi_mean"].to_numpy(dtype=float)
    medians = epitope_distance_df["d_epi_median"].to_numpy(dtype=float)
    stds = epitope_distance_df["d_epi_std"].fillna(0.0).to_numpy(dtype=float)

    ax.barh(y, means, xerr=stds, color="#72b7b2", alpha=0.9, ecolor="#4c4c4c", capsize=3)
    ax.scatter(medians, y, color="#e45756", s=26, zorder=3, label="median")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("d_epi")
    ax.set_title("Mean d_epi by epitope with run-to-run spread")
    ax.legend(frameon=False)

    stem = output_dir / "epitope_mean_distance_summary"
    fig.savefig(stem.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    runs_root = args.runs_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    run_df = collect_run_rows(runs_root)
    d_ref, d_ref_source = resolve_reference_distance(
        run_df,
        d_ref_mode=args.d_ref_mode,
        glc_knn=args.glc_knn.resolve() if args.glc_knn is not None else None,
        ylq_knn=args.ylq_knn.resolve() if args.ylq_knn is not None else None,
    )
    scored_df = append_repo_native_scores(run_df, comparison_group=args.comparison_group, d_ref=d_ref)
    summary_df = summarize_parameter_heatmap(scored_df)
    run_distance_df, epitope_distance_df = summarize_distances(scored_df)

    run_level_path = output_dir / "run_level_metric_breakdown.tsv"
    summary_path = output_dir / "redcea_possig_parameter_runs.tsv"
    run_distance_path = output_dir / "run_level_distances.tsv"
    epitope_distance_path = output_dir / "epitope_distance_summary.tsv"
    scored_df.to_csv(run_level_path, sep="\t", index=False)
    summary_df.to_csv(summary_path, sep="\t", index=False)
    run_distance_df.to_csv(run_distance_path, sep="\t", index=False)
    epitope_distance_df.to_csv(epitope_distance_path, sep="\t", index=False)

    for chain in sorted(summary_df["chain"].astype(str).unique().tolist()):
        plot_chain_heatmap(summary_df, chain=chain, output_stem=output_dir / f"redcea_possig_density_score_{chain.lower()}")
    plot_epitope_distance_histograms(scored_df, output_dir=output_dir)
    plot_epitope_distance_summary(epitope_distance_df, output_dir=output_dir)

    print(f"Discovered {len(scored_df)} run(s) from {runs_root}")
    print(f"Reference distance: d_ref={d_ref:.6f} from {d_ref_source}")
    print(f"Saved {run_level_path}")
    print(f"Saved {summary_path}")
    print(f"Saved {run_distance_path}")
    print(f"Saved {epitope_distance_path}")
    for chain in sorted(summary_df['chain'].astype(str).unique().tolist()):
        stem = output_dir / f"redcea_possig_density_score_{chain.lower()}"
        print(f"Saved {stem.with_suffix('.png')}")
        print(f"Saved {stem.with_suffix('.pdf')}")
        print(f"Saved {stem.with_suffix('.svg')}")
    print(f"Saved {(output_dir / 'epitope_mean_distance_summary.png')}")
    print(f"Saved {(output_dir / 'epitope_mean_distance_summary.pdf')}")
    print(f"Saved {(output_dir / 'epitope_mean_distance_summary.svg')}")
    print(f"Saved per-epitope histograms under {output_dir / 'distance_histograms'}")


if __name__ == "__main__":
    main()
