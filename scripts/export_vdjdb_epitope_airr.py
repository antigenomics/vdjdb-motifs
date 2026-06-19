#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import pandas as pd


CHAIN_COLS = {
    "TRA": {"cdr3": "cdr3", "v": "v.segm", "j": "j.segm", "locus": "TRA"},
    "TRB": {"cdr3": "cdr3", "v": "v.segm", "j": "j.segm", "locus": "TRB"},
}


def configure_csv_field_limit() -> None:
    max_field_limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(max_field_limit)
            return
        except OverflowError:
            max_field_limit //= 10


def is_canonical_cdr3(value: object) -> bool:
    if pd.isna(value):
        return False
    sequence = str(value).strip().upper()
    return bool(sequence) and sequence.startswith("C") and sequence.endswith(("F", "W"))


def normalize_vdjdb_columns(vdjdb_df: pd.DataFrame, chain: str) -> pd.DataFrame:
    normalized = vdjdb_df.copy()
    if "gene" in normalized.columns:
        normalized = normalized[normalized["gene"].astype(str).str.upper() == chain].copy()

    generic_targets = {
        "cdr3": {"TRA": "cdr3.alpha", "TRB": "cdr3.beta"}[chain],
        "v.segm": {"TRA": "v.alpha", "TRB": "v.beta"}[chain],
        "j.segm": {"TRA": "j.alpha", "TRB": "j.beta"}[chain],
    }
    renamed_cols: dict[str, str] = {}
    for target_col, source_col in generic_targets.items():
        if target_col not in normalized.columns and source_col in normalized.columns:
            renamed_cols[source_col] = target_col
    if renamed_cols:
        normalized = normalized.rename(columns=renamed_cols)
    return normalized


def build_airr_table(ep_df: pd.DataFrame, chain: str, *, include_noncanonical: bool) -> pd.DataFrame:
    cfg = CHAIN_COLS[chain]
    required_cols = [cfg["cdr3"], cfg["v"], cfg["j"]]
    filtered = ep_df.dropna(subset=required_cols).copy()
    if not include_noncanonical:
        filtered = filtered[filtered[cfg["cdr3"]].map(is_canonical_cdr3)].copy()

    airr = filtered[[cfg["cdr3"], cfg["v"], cfg["j"]]].copy()
    airr.columns = ["junction_aa", "v_call", "j_call"]
    airr["locus"] = cfg["locus"]
    return airr.reset_index(drop=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export one VDJdb epitope into AIRR TSV format.")
    parser.add_argument("--vdjdb", required=True, help="Path to vdjdb.slim.txt or compatible VDJdb TSV.")
    parser.add_argument("--epitope", required=True, help="Exact epitope sequence to export.")
    parser.add_argument("--chain", required=True, choices=["TRA", "TRB"], help="Chain to export.")
    parser.add_argument(
        "--species",
        default="HomoSapiens",
        help="Species filter applied before export. Use empty string to disable species filtering.",
    )
    parser.add_argument("--output", required=True, help="Output AIRR TSV path.")
    parser.add_argument(
        "--include-noncanonical",
        action="store_true",
        help="Keep non-canonical CDR3 sequences instead of filtering them out.",
    )
    return parser.parse_args()


def main() -> int:
    configure_csv_field_limit()
    args = parse_args()

    vdjdb_path = Path(args.vdjdb)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    vdjdb_df = pd.read_csv(vdjdb_path, sep="\t")
    vdjdb_df = normalize_vdjdb_columns(vdjdb_df, args.chain)

    if args.species:
        vdjdb_df = vdjdb_df[vdjdb_df["species"].astype(str) == args.species].copy()

    ep_df = vdjdb_df[vdjdb_df["antigen.epitope"].astype(str) == args.epitope].copy()
    if ep_df.empty:
        raise SystemExit(
            f"No VDJdb rows found for epitope={args.epitope!r}, chain={args.chain}, species={args.species!r}"
        )

    airr_df = build_airr_table(ep_df, args.chain, include_noncanonical=args.include_noncanonical)
    if airr_df.empty:
        raise SystemExit(
            f"Epitope {args.epitope!r} for chain={args.chain} has no clonotypes after filtering/export"
        )

    airr_df.to_csv(output_path, sep="\t", index=False)
    print(str(output_path))
    print(f"rows\t{len(airr_df)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
