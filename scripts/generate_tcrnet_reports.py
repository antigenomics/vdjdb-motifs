#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path


PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.35.2.min.js"
GRAY_COLOR = "rgba(160, 160, 160, 0.38)"
CLUSTER_COLORS = [
    "#1f77b4",
    "#d62728",
    "#2ca02c",
    "#ff7f0e",
    "#17becf",
    "#e377c2",
    "#bcbd22",
    "#8c564b",
    "#9467bd",
    "#7f7f7f",
    "#aec7e8",
    "#ff9896",
    "#98df8a",
    "#ffbb78",
    "#9edae5",
    "#f7b6d2",
    "#dbdb8d",
    "#c49c94",
    "#c5b0d5",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate per-epitope TCRNET HTML reports.")
    parser.add_argument(
        "--input",
        default="results/tcrnet/cluster_members.txt",
        help="Path to the TCRNET cluster members TSV.",
    )
    parser.add_argument(
        "--output-dir",
        default="results/tcrnet/viz",
        help="Directory where per-epitope HTML reports will be written.",
    )
    return parser.parse_args()


def sanitize_filename_token(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    cleaned = cleaned.strip("._-")
    return cleaned or "item"


def load_rows(input_path: Path) -> list[dict[str, str]]:
    with input_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [{str(key).strip(): value for key, value in row.items()} for row in reader]
    for row in rows:
        row["x"] = float(row["x"])
        row["y"] = float(row["y"])
        row["csz"] = int(float(row["csz"]))
    return rows


def build_hover_text(row: dict[str, str]) -> str:
    fields = [
        ("Epitope", row["antigen.epitope"]),
        ("Cluster", row["cid"]),
        ("CDR3", row["cdr3aa"]),
        ("V", row.get("v.segm.repr") or row.get("v.segm") or ""),
        ("J", row.get("j.segm.repr") or row.get("j.segm") or ""),
        ("Cluster size", str(row["csz"])),
        ("Species", row["species"]),
        ("Chain", row["gene"]),
        ("Antigen gene", row.get("antigen.gene", "")),
        ("Antigen species", row.get("antigen.species", "")),
        ("MHC A", row.get("mhc.a", "")),
        ("MHC B", row.get("mhc.b", "")),
        ("MHC class", row.get("mhc.class", "")),
    ]
    return "<br>".join(f"<b>{label}:</b> {value}" for label, value in fields if value)


def build_trace(name: str, rows: list[dict[str, str]], color: str, *, size: int, opacity: float) -> dict[str, object]:
    return {
        "type": "scattergl",
        "mode": "markers",
        "name": name,
        "x": [row["x"] for row in rows],
        "y": [row["y"] for row in rows],
        "text": [build_hover_text(row) for row in rows],
        "hovertemplate": "%{text}<extra></extra>",
        "marker": {
            "color": color,
            "size": size,
            "opacity": opacity,
            "line": {"width": 0},
        },
    }


def build_payload(*, epitope: str, rows_same_group: list[dict[str, str]], rows_epitope: list[dict[str, str]]) -> dict[str, object]:
    traces: list[dict[str, object]] = []
    other_rows = [row for row in rows_same_group if row["antigen.epitope"] != epitope]
    if other_rows:
        traces.append(build_trace("Other TCRs", other_rows, GRAY_COLOR, size=7, opacity=0.55))

    clusters: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows_epitope:
        clusters[row["cid"]].append(row)
    for idx, cluster_id in enumerate(sorted(clusters)):
        traces.append(build_trace(cluster_id, clusters[cluster_id], CLUSTER_COLORS[idx % len(CLUSTER_COLORS)], size=9, opacity=0.92))

    return {
        "traces": traces,
        "layout": {
            "title": {"text": ""},
            "template": "plotly_white",
            "width": 1100,
            "height": 760,
            "hovermode": "closest",
            "legend": {"orientation": "v", "x": 1.02, "xanchor": "left", "y": 1, "yanchor": "top"},
            "margin": {"l": 40, "r": 20, "t": 30, "b": 40},
            "xaxis": {"title": "", "showgrid": False, "zeroline": False, "showticklabels": False},
            "yaxis": {"title": "", "showgrid": False, "zeroline": False, "showticklabels": False, "scaleanchor": "x", "scaleratio": 1},
        },
        "config": {"displaylogo": False, "responsive": True},
    }


def build_html(*, species: str, gene: str, epitope: str, rows_same_group: list[dict[str, str]], rows_epitope: list[dict[str, str]]) -> str:
    title = epitope
    div_id = f"plot-{sanitize_filename_token(species)}-{sanitize_filename_token(gene)}-{sanitize_filename_token(epitope)}"
    payload = build_payload(epitope=epitope, rows_same_group=rows_same_group, rows_epitope=rows_epitope)
    payload_json = json.dumps(payload, ensure_ascii=False)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{title}</title>
  <script src="{PLOTLY_CDN}"></script>
</head>
<body>
  <div id="{div_id}"></div>
  <script>
    const payload = {payload_json};
    Plotly.newPlot("{div_id}", payload.traces, payload.layout, payload.config);
  </script>
</body>
</html>
"""


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = load_rows(input_path)
    rows_by_group: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    rows_by_epitope: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)

    for row in rows:
        group_key = (row["species"], row["gene"])
        epitope_key = (row["species"], row["gene"], row["antigen.epitope"])
        rows_by_group[group_key].append(row)
        rows_by_epitope[epitope_key].append(row)

    for (species, gene, epitope), rows_epitope in sorted(rows_by_epitope.items()):
        rows_same_group = rows_by_group[(species, gene)]
        html = build_html(
            species=species,
            gene=gene,
            epitope=epitope,
            rows_same_group=rows_same_group,
            rows_epitope=rows_epitope,
        )
        filename = (
            f"{sanitize_filename_token(species)}_"
            f"{sanitize_filename_token(gene)}_"
            f"{sanitize_filename_token(epitope)}.html"
        )
        (output_dir / filename).write_text(html, encoding="utf-8")


if __name__ == "__main__":
    main()
