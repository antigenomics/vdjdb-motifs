#!/usr/bin/env python
from __future__ import annotations

import argparse
import math
import time
from pathlib import Path

import pandas as pd


LEGACY_COLUMNS = [
    "species",
    "antigen.epitope",
    "gene",
    "aa",
    "pos",
    "len",
    "v.segm.repr",
    "j.segm.repr",
    "cid",
    "csz",
    "count",
    "count.bg",
    "total.bg",
    "count.bg.i",
    "total.bg.i",
    "need.impute",
    "freq",
    "freq.bg",
    "I",
    "I.norm",
    "height.I",
    "height.I.norm",
    "antigen.gene",
    "antigen.species",
    "mhc.a",
    "mhc.b",
    "mhc.class",
]

REQUIRED_CLUSTER_COLUMNS = {
    "species",
    "antigen.epitope",
    "antigen.gene",
    "antigen.species",
    "mhc.a",
    "mhc.b",
    "mhc.class",
    "gene",
    "cdr3aa",
    "cid",
    "v.segm.repr",
    "j.segm.repr",
}

START_TIME = time.perf_counter()


def log(message: str) -> None:
    elapsed = time.perf_counter() - START_TIME
    print(f"[motif_pwms +{elapsed:8.1f}s] {message}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute legacy-style motif_pwms.txt tables from one or more "
            "cluster_members tables plus AIRR-like background repertoires."
        )
    )
    parser.add_argument(
        "--cluster-members",
        nargs="+",
        required=True,
        help="One or more cluster_members TSV files.",
    )
    parser.add_argument(
        "--background",
        nargs=3,
        action="append",
        metavar=("SPECIES", "GENE", "FILE"),
        required=True,
        help=(
            "Precomputed background PWM TSV for a species/chain pair in legacy "
            "VDJtools format. Example: --background HomoSapiens TRA "
            "tcrnet/pwms/human.tra.aa.aa_cdr3_pwm.txt"
        ),
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Destination motif_pwms.txt path.",
    )
    return parser.parse_args()


def read_cluster_members(paths: list[str]) -> pd.DataFrame:
    tables = []
    for path_str in paths:
        path = Path(path_str)
        log(f"Reading cluster_members from {path}")
        table = pd.read_csv(path, sep="\t")
        missing = REQUIRED_CLUSTER_COLUMNS - set(table.columns)
        if missing:
            raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
        log(f"Loaded {len(table):,} rows from {path}")
        tables.append(table)

    df = pd.concat(tables, ignore_index=True)
    df = df.loc[df["cdr3aa"].notna() & df["cid"].notna()].copy()
    df["cdr3aa"] = df["cdr3aa"].astype(str)
    df = df.loc[df["cdr3aa"].str.len() > 0].copy()
    log(
        "Prepared cluster_members table with "
        f"{len(df):,} rows, {df['cid'].nunique():,} clusters, "
        f"{df['cdr3aa'].nunique():,} unique CDR3s"
    )
    return df


def build_cluster_pwm(df: pd.DataFrame) -> pd.DataFrame:
    t0 = time.perf_counter()
    cluster_meta = (
        df.groupby("cid", as_index=False)
        .agg(
            {
                "species": "first",
                "antigen.epitope": "first",
                "antigen.gene": "first",
                "antigen.species": "first",
                "mhc.a": "first",
                "mhc.b": "first",
                "mhc.class": "first",
                "gene": "first",
                "v.segm.repr": "first",
                "j.segm.repr": "first",
            }
        )
        .rename(columns={"cid": "_cid"})
    )
    cluster_sizes = df.groupby("cid").size().rename("csz").reset_index().rename(columns={"cid": "_cid"})
    cluster_meta = cluster_meta.merge(cluster_sizes, on="_cid", how="left").rename(columns={"_cid": "cid"})

    seqs = df[["cid", "cdr3aa"]].drop_duplicates().copy()
    log(f"Expanding {len(seqs):,} unique cluster-member sequences into position rows")
    aa_rows: list[tuple[str, int, int, str]] = []
    for cid, cdr3 in seqs.itertuples(index=False):
        seq_len = len(cdr3)
        for pos, aa in enumerate(cdr3):
            aa_rows.append((cid, seq_len, pos, aa))
    aa_df = pd.DataFrame(aa_rows, columns=["cid", "len", "pos", "aa"])
    log(f"Expanded cluster members to {len(aa_df):,} position rows")

    counts = aa_df.groupby(["cid", "len", "pos", "aa"], as_index=False).size().rename(columns={"size": "count"})
    out = counts.merge(cluster_meta, on="cid", how="left")
    log(
        "Built cluster PWM counts with "
        f"{len(out):,} rows in {time.perf_counter() - t0:.1f}s"
    )
    return out


def build_background_tables(background_specs: list[list[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    exact_tables = []
    imputed_tables = []

    for species, gene, path_str in background_specs:
        path = Path(path_str)
        t0 = time.perf_counter()
        log(f"Reading precomputed background PWM for {species} {gene} from {path}")
        bg = pd.read_csv(path, sep="\t")
        required = {"v", "j", "len", "pos", "aa", "value"}
        missing = required - set(bg.columns)
        if missing:
            raise ValueError(f"{path} is missing required background columns: {sorted(missing)}")

        exact = bg[["v", "j", "len", "pos", "aa", "value"]].copy()
        exact["species"] = species
        exact["gene"] = gene
        exact = exact.rename(columns={"v": "v.segm.repr", "j": "j.segm.repr", "value": "count.bg"})
        exact["count.bg"] = exact["count.bg"].astype(float)
        log(f"Background {species} {gene}: loaded {len(exact):,} PWM rows")

        exact_totals = (
            exact.loc[exact["pos"] == 0]
            .groupby(["species", "gene", "v.segm.repr", "j.segm.repr", "len"], as_index=False)["count.bg"]
            .sum()
            .rename(columns={"count.bg": "total.bg"})
        )
        exact_tables.append(
            exact.merge(exact_totals, on=["species", "gene", "v.segm.repr", "j.segm.repr", "len"], how="left")
        )

        imputed = (
            exact.groupby(["species", "gene", "len", "pos", "aa"], as_index=False)["count.bg"]
            .sum()
            .rename(columns={"count.bg": "count.bg.i"})
        )
        imputed_totals = (
            imputed.loc[imputed["pos"] == 0]
            .groupby(["species", "gene", "len"], as_index=False)["count.bg.i"]
            .sum()
            .rename(columns={"count.bg.i": "total.bg.i"})
        )
        imputed_tables.append(imputed.merge(imputed_totals, on=["species", "gene", "len"], how="left"))
        log(
            f"Background {species} {gene}: exact rows={len(exact):,}, "
            f"imputed rows={len(imputed):,}, done in {time.perf_counter() - t0:.1f}s"
        )

    exact_df = pd.concat(exact_tables, ignore_index=True)
    imputed_df = pd.concat(imputed_tables, ignore_index=True)
    log(
        f"Combined background tables: exact={len(exact_df):,} rows, "
        f"imputed={len(imputed_df):,} rows"
    )
    return exact_df, imputed_df


def compute_information(df: pd.DataFrame) -> pd.DataFrame:
    log20 = math.log(20.0)
    group_cols = [
        "species",
        "antigen.epitope",
        "gene",
        "cid",
        "csz",
        "v.segm.repr",
        "j.segm.repr",
        "pos",
        "len",
    ]
    log(
        "Computing information metrics for "
        f"{len(df):,} merged rows across {df[group_cols].drop_duplicates().shape[0]:,} position groups"
    )

    def per_position(group: pd.DataFrame) -> pd.DataFrame:
        freq = group["count"] / group["csz"]
        if bool(group["need.impute"].iloc[0]):
            freq_bg = (group["count.bg.i"] + 1.0) / (group["total.bg.i"] + 1.0)
        else:
            freq_bg = (group["count.bg"] + 1.0) / (group["total.bg"] + 1.0)

        info = 1.0 + (freq * freq.map(math.log)).sum() / log20
        info_norm = -((freq * freq_bg.map(math.log)).sum()) / log20 / 2.0

        out = group.copy()
        out["freq"] = freq
        out["freq.bg"] = freq_bg
        out["I"] = info
        out["I.norm"] = info_norm
        out["height.I"] = freq * info
        out["height.I.norm"] = freq * info_norm
        return out

    t0 = time.perf_counter()
    out = df.groupby(group_cols, group_keys=False).apply(per_position).reset_index(drop=True)
    log(f"Computed information metrics in {time.perf_counter() - t0:.1f}s")
    return out


def main() -> None:
    args = parse_args()

    log("Starting motif PWM computation")
    cluster_df = read_cluster_members(args.cluster_members)
    cluster_pwm = build_cluster_pwm(cluster_df)

    bg_exact, bg_imputed = build_background_tables(args.background)
    supported_pairs = set(bg_imputed[["species", "gene"]].drop_duplicates().itertuples(index=False, name=None))
    available_pairs = set(cluster_pwm[["species", "gene"]].drop_duplicates().itertuples(index=False, name=None))
    unsupported_pairs = sorted(available_pairs - supported_pairs)
    if unsupported_pairs:
        print(f"Skipping clusters without matching background: {unsupported_pairs}")

    log(f"Supported species/gene pairs: {sorted(supported_pairs)}")
    cluster_pwm = cluster_pwm.loc[
        cluster_pwm.apply(lambda row: (row["species"], row["gene"]) in supported_pairs, axis=1)
    ].copy()
    if cluster_pwm.empty:
        raise RuntimeError("No cluster rows remain after background species/gene filtering.")
    log(f"{len(cluster_pwm):,} cluster PWM rows remain after species/gene filtering")

    t0 = time.perf_counter()
    merged = cluster_pwm.merge(
        bg_exact,
        on=["species", "gene", "v.segm.repr", "j.segm.repr", "len", "pos", "aa"],
        how="left",
    ).merge(
        bg_imputed,
        on=["species", "gene", "len", "pos", "aa"],
        how="left",
    )
    log(f"Merged cluster/background tables into {len(merged):,} rows in {time.perf_counter() - t0:.1f}s")

    # Keep only rows with an exact V/J/len background, matching the legacy Rmd:
    # df.pwms %>% merge(df.bg.pwms, all.x = T) %>% ... %>% filter(total.bg > 0)
    merged = merged.loc[merged["total.bg"].fillna(0) > 0].copy()
    merged["need.impute"] = False
    merged["count.bg"] = merged["count.bg"].fillna(0).astype(int)
    merged["count.bg.i"] = merged["count.bg.i"].fillna(0).astype(int)
    merged["total.bg"] = merged["total.bg"].fillna(0).astype(int)
    merged["total.bg.i"] = merged["total.bg.i"].fillna(0).astype(int)
    log(f"{len(merged):,} merged rows remain after exact-background filtering")

    pwm = compute_information(merged)
    pwm = pwm.sort_values(["species", "antigen.epitope", "gene", "cid", "len", "pos", "aa"]).copy()
    pwm = pwm[LEGACY_COLUMNS]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pwm.to_csv(output_path, sep="\t", index=False)
    log(f"Wrote {len(pwm):,} motif PWM rows to {output_path}")


if __name__ == "__main__":
    main()
