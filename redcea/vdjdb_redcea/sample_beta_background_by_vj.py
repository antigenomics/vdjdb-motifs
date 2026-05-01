from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


SOURCE_CDR3_COLUMNS = ("junction_aa", "cdr3aa", "cdr3", "cdr3_beta_aa", "cdr3.beta")
SOURCE_V_COLUMNS_BY_CHAIN = {
    "TRA": ("v_call", "v", "TRAV", "v_alpha", "v.alpha", "v.segm"),
    "TRB": ("v_call", "v", "TRBV", "v_beta", "v.beta", "v.segm"),
}
SOURCE_J_COLUMNS_BY_CHAIN = {
    "TRA": ("j_call", "j", "TRAJ", "j_alpha", "j.alpha", "j.segm"),
    "TRB": ("j_call", "j", "TRBJ", "j_beta", "j.beta", "j.segm"),
}
SOURCE_LOCUS_COLUMNS = ("locus", "gene", "chain")
TARGET_GENE_BY_CHAIN = {"TRA": "TRA", "TRB": "TRB"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sample a chain-specific background repertoire so its V/J usage matches "
            "the selected subset of vdjdb.slim.txt."
        )
    )
    parser.add_argument(
        "--chain",
        choices=("TRA", "TRB"),
        default="TRB",
        help="TCR chain to prepare.",
    )
    parser.add_argument("--source", required=True, help="Path to the large source background table.")
    parser.add_argument("--target", required=True, help="Path to vdjdb.slim.txt.")
    parser.add_argument("--output", required=True, help="Where to save the sampled beta background TSV.")
    parser.add_argument(
        "--source-embedding",
        default=None,
        help="Optional path to source background embeddings parquet aligned row-wise with --source.",
    )
    parser.add_argument(
        "--output-embedding",
        default=None,
        help="Optional path where sampled background embeddings parquet should be written.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducible sampling.",
    )
    parser.add_argument(
        "--keep-alleles",
        action="store_true",
        help="Match V/J with alleles preserved. Default behavior strips allele suffixes on both sides.",
    )
    parser.add_argument(
        "--allow-shortfall",
        action="store_true",
        help=(
            "Allow output to be smaller than requested when some V/J pairs are missing "
            "or too small in the source background."
        ),
    )
    return parser.parse_args()


def _read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    return pd.read_csv(path, sep=None, engine="python")


def _pick_column(frame: pd.DataFrame, candidates: tuple[str, ...], label: str) -> str:
    for column in candidates:
        if column in frame.columns:
            return column
    raise ValueError(f"Could not find a {label} column. Tried: {', '.join(candidates)}")


def _strip_allele(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.replace(r"\*.*$", "", regex=True)


def _normalize_gene(series: pd.Series, keep_alleles: bool) -> pd.Series:
    cleaned = series.astype(str).str.strip()
    return cleaned if keep_alleles else _strip_allele(cleaned)


def _normalize_source_background(frame: pd.DataFrame, chain: str, keep_alleles: bool) -> pd.DataFrame:
    cdr3_col = _pick_column(frame, SOURCE_CDR3_COLUMNS, "CDR3")
    v_col = _pick_column(frame, SOURCE_V_COLUMNS_BY_CHAIN[chain], "V")
    j_col = _pick_column(frame, SOURCE_J_COLUMNS_BY_CHAIN[chain], "J")

    out = pd.DataFrame(
        {
            "source_row_index": frame.index.to_numpy(),
            "junction_aa": frame[cdr3_col].astype(str).str.strip(),
            "v_call": frame[v_col].astype(str).str.strip(),
            "j_call": frame[j_col].astype(str).str.strip(),
        }
    )

    locus_col = next((column for column in SOURCE_LOCUS_COLUMNS if column in frame.columns), None)
    if locus_col is None:
        out["locus"] = "beta"
    else:
        locus = frame[locus_col].astype(str).str.strip().str.upper()
        out["locus"] = locus

    out = out[out["junction_aa"].ne("") & out["v_call"].ne("") & out["j_call"].ne("")].copy()
    out = out[~out["junction_aa"].str.lower().eq("nan")].copy()
    allowed_locus = {"TRA": {"TRA", "ALPHA", "TCRA"}, "TRB": {"TRB", "BETA", "TCRB"}}[chain]
    out = out[out["locus"].isin(allowed_locus)].copy()
    out["locus"] = "alpha" if chain == "TRA" else "beta"
    out["v_norm"] = _normalize_gene(out["v_call"], keep_alleles=keep_alleles)
    out["j_norm"] = _normalize_gene(out["j_call"], keep_alleles=keep_alleles)
    out = out.drop_duplicates(subset=["junction_aa", "v_call", "j_call"]).reset_index(drop=True)
    return out


def _load_target_vj_counts(frame: pd.DataFrame, chain: str, keep_alleles: bool) -> pd.DataFrame:
    required = {"gene", "v.segm", "j.segm"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Target table is missing required columns: {sorted(missing)}")

    target = frame.loc[
        frame["gene"].astype(str).str.upper().eq(TARGET_GENE_BY_CHAIN[chain]),
        ["v.segm", "j.segm"],
    ].copy()
    target = target.dropna().reset_index(drop=True)
    if target.empty:
        raise ValueError(
            f"Target vdjdb.slim.txt contains no {chain} rows with both V and J segments."
        )

    target["v_norm"] = _normalize_gene(target["v.segm"], keep_alleles=keep_alleles)
    target["j_norm"] = _normalize_gene(target["j.segm"], keep_alleles=keep_alleles)
    counts = (
        target.groupby(["v_norm", "j_norm"])
        .size()
        .reset_index(name="target_count")
        .sort_values(["target_count", "v_norm", "j_norm"], ascending=[False, True, True])
        .reset_index(drop=True)
    )
    return counts


def sample_background_by_vj(
    source_path: str | Path,
    target_path: str | Path,
    output_path: str | Path,
    chain: str = "TRB",
    source_embedding_path: str | Path | None = None,
    output_embedding_path: str | Path | None = None,
    seed: int = 42,
    keep_alleles: bool = False,
    allow_shortfall: bool = False,
) -> pd.DataFrame:
    source_raw = _read_table(source_path)
    target_raw = _read_table(target_path)

    source = _normalize_source_background(source_raw, chain=chain, keep_alleles=keep_alleles)
    target_counts = _load_target_vj_counts(target_raw, chain=chain, keep_alleles=keep_alleles)
    source_embeddings = None
    if source_embedding_path is not None:
        source_embeddings = pd.read_parquet(source_embedding_path)
        if len(source_embeddings) != len(source_raw):
            raise ValueError(
                "Source background embeddings must have the same number of rows as the source background table."
            )

    sampled_parts: list[pd.DataFrame] = []
    shortfalls: list[str] = []

    for offset, row in target_counts.iterrows():
        pair_pool = source[(source["v_norm"] == row["v_norm"]) & (source["j_norm"] == row["j_norm"])].copy()
        available = len(pair_pool)
        needed = int(row["target_count"])

        if available < needed:
            message = (
                f"{row['v_norm']} / {row['j_norm']}: requested {needed}, available {available}"
            )
            if not allow_shortfall:
                raise ValueError(f"Not enough source clonotypes for V/J pair {message}")
            shortfalls.append(message)
            needed = available

        if needed == 0:
            continue

        sampled_parts.append(
            pair_pool.sample(n=needed, replace=False, random_state=seed + int(offset)).copy()
        )

    if not sampled_parts:
        raise ValueError("Sampling produced no rows.")

    sampled = pd.concat(sampled_parts, ignore_index=True)
    sampled_embeddings = None
    if source_embeddings is not None:
        sampled_embeddings = source_embeddings.iloc[sampled["source_row_index"].to_numpy()].reset_index(drop=True)

    sampled = sampled[["junction_aa", "v_call", "j_call", "locus"]].reset_index(drop=True)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sampled.to_csv(output_path, sep="\t", index=False)
    if output_embedding_path is not None:
        if sampled_embeddings is None:
            raise ValueError("--output-embedding requires --source-embedding.")
        output_embedding_path = Path(output_embedding_path)
        output_embedding_path.parent.mkdir(parents=True, exist_ok=True)
        sampled_embeddings.to_parquet(output_embedding_path, index=False)

    print(f"Chain: {chain}")
    print(f"Source clonotypes: {len(source)}")
    print(f"Target {chain} rows: {int(target_counts['target_count'].sum())}")
    print(f"Sampled rows: {len(sampled)}")
    print(f"Unique target V/J pairs: {len(target_counts)}")
    if shortfalls:
        print("Shortfalls:")
        for item in shortfalls:
            print(f"  - {item}")
    print(f"Saved: {output_path}")
    if output_embedding_path is not None:
        print(f"Saved embeddings: {output_embedding_path}")
    return sampled


def main() -> None:
    args = parse_args()
    sample_background_by_vj(
        source_path=args.source,
        target_path=args.target,
        output_path=args.output,
        chain=args.chain,
        source_embedding_path=args.source_embedding,
        output_embedding_path=args.output_embedding,
        seed=args.seed,
        keep_alleles=args.keep_alleles,
        allow_shortfall=args.allow_shortfall,
    )


if __name__ == "__main__":
    main()
