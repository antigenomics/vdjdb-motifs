#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Select the largest VDJdb epitopes by clonotype count for one or more chains."
    )
    parser.add_argument(
        "--vdjdb",
        type=Path,
        default=Path("vdjdb_release") / "vdjdb.slim.txt",
        help="Path to vdjdb.slim.txt with generic gene/cdr3/species/antigen.epitope columns.",
    )
    parser.add_argument(
        "--chains",
        nargs="+",
        default=["TRB"],
        help="One or more chains to analyze, e.g. TRB or TRA TRB.",
    )
    parser.add_argument(
        "--species",
        default="HomoSapiens",
        help="Species filter applied before counting epitopes.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=10,
        help="Return top-N epitopes per chain.",
    )
    parser.add_argument(
        "--min-count",
        type=int,
        default=1,
        help="Optional minimum clonotype count threshold.",
    )
    parser.add_argument(
        "--include-noncanonical",
        action="store_true",
        help="If set, keep non-canonical CDR3s instead of filtering to C...F/W sequences.",
    )
    parser.add_argument(
        "--format",
        choices=["epitopes", "tsv", "json"],
        default="epitopes",
        help="Output format. `epitopes` prints one epitope per line for a single chain.",
    )
    return parser.parse_args()


def is_canonical_cdr3(value: str) -> bool:
    text = str(value).strip().upper()
    return bool(text) and text.startswith("C") and text.endswith(("F", "W"))


def load_top_epitopes(
    path: Path,
    *,
    chains: list[str],
    species: str,
    top_n: int,
    min_count: int,
    include_noncanonical: bool,
) -> list[dict[str, object]]:
    max_field_limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(max_field_limit)
            break
        except OverflowError:
            max_field_limit = max_field_limit // 10
    chains_upper = [str(chain).upper() for chain in chains]
    counters = {chain: Counter() for chain in chains_upper}

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene", "cdr3", "species", "antigen.epitope"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise KeyError("Missing required VDJdb columns: {0}".format(", ".join(sorted(missing))))

        for row in reader:
            chain = str(row.get("gene", "")).strip().upper()
            if chain not in counters:
                continue
            if str(row.get("species", "")).strip() != species:
                continue
            epitope = str(row.get("antigen.epitope", "")).strip()
            cdr3 = str(row.get("cdr3", "")).strip()
            if not epitope or not cdr3:
                continue
            if (not include_noncanonical) and (not is_canonical_cdr3(cdr3)):
                continue
            counters[chain][epitope] += 1

    rows: list[dict[str, object]] = []
    for chain in chains_upper:
        ranked = sorted(counters[chain].items(), key=lambda item: (-item[1], item[0]))
        for epitope, count in ranked[:top_n]:
            if int(count) < int(min_count):
                continue
            rows.append({"chain": chain, "epitope": epitope, "count": int(count)})
    return rows


def main() -> None:
    args = parse_args()
    rows = load_top_epitopes(
        args.vdjdb,
        chains=args.chains,
        species=args.species,
        top_n=int(args.top_n),
        min_count=int(args.min_count),
        include_noncanonical=bool(args.include_noncanonical),
    )

    if args.format == "json":
        print(json.dumps(rows, indent=2))
        return

    if args.format == "tsv":
        print("chain\tepitope\tcount")
        for row in rows:
            print("{chain}\t{epitope}\t{count}".format(**row))
        return

    if len({row["chain"] for row in rows}) > 1:
        raise ValueError("`--format epitopes` is only valid for a single chain. Use --format tsv or json instead.")
    for row in rows:
        print(str(row["epitope"]))


if __name__ == "__main__":
    main()
