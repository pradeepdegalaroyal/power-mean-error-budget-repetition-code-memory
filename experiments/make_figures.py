"""Re-draw all figures from the saved ``results/*.csv`` and ``summary.json``."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_error_budget.plots import make_all_figures  # noqa: E402

NAMES = ("budget_accuracy", "validity_map", "decoder_ablation", "decoder_mismatch",
         "shot_pareto", "threshold_scan", "exponent_ablation")

if __name__ == "__main__":
    tables = {n: pd.read_csv(ROOT / "results" / f"{n}.csv") for n in NAMES}
    summary = json.loads((ROOT / "results" / "summary.json").read_text())
    make_all_figures(tables, summary, ROOT / "figures")
