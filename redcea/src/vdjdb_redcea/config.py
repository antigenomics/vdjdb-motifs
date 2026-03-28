from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


CHAIN_COLS = {
    "TRA": {"cdr3": "cdr3.alpha", "v": "v.alpha", "j": "j.alpha", "locus": "alpha", "gene": "alpha"},
    "TRB": {"cdr3": "cdr3.beta", "v": "v.beta", "j": "j.beta", "locus": "beta", "gene": "beta"},
}

DEFAULT_PLOT_BG_POINTS = 100_000


@dataclass(frozen=True)
class OutputPaths:
    output_root: Path
    viz_dir: Path
    airr_dir: Path
    tcremp_dir: Path
    tcrempnet_dir: Path
