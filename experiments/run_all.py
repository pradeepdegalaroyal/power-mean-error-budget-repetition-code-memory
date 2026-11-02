"""Reproduce every result and figure: ``python experiments/run_all.py [--quick]``."""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from qec_error_budget import studies  # noqa: E402
from qec_error_budget.plots import make_all_figures  # noqa: E402
from qec_error_budget.summary import summarise  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="smoke run (< 1 minute)")
    parser.add_argument("--out", default=None, help="output root (default: repository root)")
    args = parser.parse_args()

    cfg = studies.Config.quick() if args.quick else studies.Config()
    out = Path(args.out) if args.out else (Path(tempfile.mkdtemp(prefix="qeb_quick_")) if args.quick else ROOT)
    results, figures = out / "results", out / "figures"
    results.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    timings: dict[str, float] = {}
    tables = {}
    for name, fn in (
        ("budget_accuracy", studies.budget_accuracy),
        ("validity_map", studies.validity_map),
        ("decoder_ablation", studies.decoder_ablation),
        ("decoder_mismatch", studies.decoder_mismatch),
        ("shot_pareto", studies.shot_pareto),
        ("threshold_scan", studies.threshold_scan),
    ):
        t0 = time.time()
        tables[name] = fn(cfg)
        timings[name] = round(time.time() - t0, 1)
        tables[name].to_csv(results / f"{name}.csv", index=False)
        print(f"{name}: {timings[name]} s", flush=True)

    tables["exponent_ablation"] = studies.exponent_ablation(tables["budget_accuracy"])
    tables["exponent_ablation"].to_csv(results / "exponent_ablation.csv", index=False)

    summary = summarise(tables)
    summary["runtime_seconds"] = timings
    summary["quick"] = bool(args.quick)
    (results / "summary.json").write_text(json.dumps(summary, indent=2))
    make_all_figures(tables, summary, figures)
    print("total", round(sum(timings.values()), 1), "s; wrote", results, "and", figures)


if __name__ == "__main__":
    main()
