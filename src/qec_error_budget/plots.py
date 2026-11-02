"""Figures: one accent colour (burnt orange) on greys, 200 dpi."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ACCENT = "#c2410c"
DARK, MID, LIGHT = "#3a3a3a", "#8c8c8c", "#cfcfcf"
DPI = 200
TOL = 0.25


def _style() -> None:
    plt.rcParams.update({
        "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
        "axes.titlesize": 10, "axes.labelsize": 9, "figure.dpi": 100,
    })


def _save(fig: plt.Figure, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def fig_accuracy_violin(acc: pd.DataFrame, path: Path) -> None:
    """Per-seed log10(prediction / Monte Carlo) for three budget predictors."""
    methods = [("additive", "additive (one source at a time)", DARK),
               ("power_mean", "power-mean (this work)", ACCENT),
               ("classical_majority", "lumped classical majority vote", MID)]
    dists = sorted(acc["distance"].unique())
    fig, axes = plt.subplots(1, len(dists), figsize=(9.5, 4), sharey=True,
                             gridspec_kw={"width_ratios": [acc[acc.distance == d]["scale"].nunique() for d in dists]})
    axes = np.atleast_1d(axes)
    rng = np.random.default_rng(0)
    for ax, d in zip(axes, dists):
        sub = acc[acc.distance == d]
        scales = sorted(sub["scale"].unique())
        for gi, s in enumerate(scales):
            g = sub[sub.scale == s]
            for mi, (key, _, col) in enumerate(methods):
                vals = np.log10(g[f"pred_{key}"] / g["eps_all"].replace(0, np.nan)).replace(
                    [np.inf, -np.inf], np.nan).dropna().values
                pos = gi + (mi - 1) * 0.27
                if len(vals) >= 2 and np.ptp(vals) > 0:
                    v = ax.violinplot(vals, [pos], widths=0.25, showextrema=False)
                    for b in v["bodies"]:
                        b.set_facecolor(col)
                        b.set_alpha(0.35)
                        b.set_edgecolor(col)
                ax.scatter(pos + rng.uniform(-0.04, 0.04, len(vals)), vals, s=9, color=col, zorder=3,
                           label=methods[mi][1] if gi == 0 and d == dists[0] else None)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(range(len(scales)))
        ax.set_xticklabels([f"{100 * s:g}%" for s in scales])
        ax.set_xlabel("noise scale s (gate error, per location)")
        ax.set_title(f"repetition code, distance {d} ({d} rounds)")
    axes[0].set_ylabel("log10 (predicted / Monte Carlo logical error)")
    axes[0].legend(frameon=False, fontsize=7.5, loc="lower right")
    _save(fig, path)


def fig_validity_map(vm: pd.DataFrame, path: Path) -> None:
    """Filled-contour map of budget error over noise scale and readout ratio."""
    pooled = vm.groupby(["scale", "ratio"]).mean(numeric_only=True).reset_index()
    scales = np.sort(pooled["scale"].unique())
    ratios = np.sort(pooled["ratio"].unique())
    X, Y = np.meshgrid(np.log10(scales), np.log2(ratios))
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4), sharey=True)
    levels = np.linspace(0, 0.6, 13)
    cs = None
    for ax, key, title in zip(axes, ("additive", "power_mean"),
                              ("additive budget", "power-mean budget")):
        Z = np.zeros_like(X)
        for i, r in enumerate(ratios):
            for j, s in enumerate(scales):
                row = pooled[(pooled.scale == s) & (pooled.ratio == r)].iloc[0]
                Z[i, j] = abs(np.log10(row[f"pred_{key}"] / row["eps_all"]))
        cs = ax.contourf(X, Y, np.clip(Z, 0, 0.6), levels=levels, cmap="Greys", extend="max")
        if Z.min() < np.log10(1 + TOL) < Z.max():
            ax.contour(X, Y, Z, levels=[np.log10(1 + TOL)], colors=ACCENT, linewidths=2)
        ax.set_title(title + "  (orange line: 25 % error)")
        ax.set_xlabel("noise scale s")
        ax.set_xticks(np.log10(scales))
        ax.set_xticklabels([f"{100 * s:g}%" for s in scales], fontsize=7)
        ax.plot(np.log10(0.004), 1.0, marker="*", color=ACCENT, ms=13, mec="k", zorder=5)
    axes[0].set_ylabel("readout / gate error ratio")
    axes[0].set_yticks(np.log2(ratios))
    axes[0].set_yticklabels([f"{r:g}" for r in ratios])
    fig.colorbar(cs, ax=axes, label="|log10 (prediction / Monte Carlo)|", shrink=0.9)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def fig_pareto(par: pd.DataFrame, path: Path) -> None:
    """Prediction error versus Monte Carlo cost, Pareto front highlighted."""
    g = par.dropna(subset=["log10_error"]).groupby(["method", "total_shots"])["log10_error"].agg(
        lambda x: float(np.sqrt(np.mean(np.square(x))))).reset_index()
    g.columns = ["method", "shots", "rms"]
    marks = {"direct Monte Carlo": ("s", MID), "additive budget": ("^", DARK), "power-mean budget": ("o", ACCENT)}
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    pts = g.sort_values("shots").values
    front, best = [], np.inf
    for m, n, r in pts:
        if r < best:
            front.append((n, r))
            best = r
    for m, (mk, col) in marks.items():
        s = g[g.method == m]
        ax.scatter(s.shots, s.rms, marker=mk, s=45, color=col, edgecolor="k", linewidth=0.5, label=m, zorder=3)
    ax.step([f[0] for f in front], [f[1] for f in front], where="post", color=ACCENT, lw=2, zorder=2,
            label="Pareto front")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("total Monte Carlo shots spent")
    ax.set_ylabel("RMS log10 error of combined logical error")
    ax.set_title("Cost versus accuracy, distance 3, s = 0.8 %")
    ax.legend(frameon=False, fontsize=8)
    _save(fig, path)


def fig_ablation(abl: pd.DataFrame, mm: pd.DataFrame, path: Path) -> None:
    """Decoder-weight ablation (box + seeds) and decoder mis-calibration sweep."""
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4))
    ax = axes[0]
    sub = abl[abl.distance == 3]
    modes = [("per_run", "matched per run"), ("fixed", "one fixed model"), ("uniform", "unit weights")]
    rng = np.random.default_rng(1)
    for i, (mode, label) in enumerate(modes):
        for key, off, col in (("additive", -0.18, DARK), ("power_mean", 0.18, ACCENT)):
            g = sub[sub.decoder == mode]
            vals = (g[f"pred_{key}"] / g["eps_all"].replace(0, np.nan)).dropna().values
            ax.boxplot(vals, positions=[i + off], widths=0.28, showfliers=False,
                       medianprops={"color": "k"}, boxprops={"color": col}, whiskerprops={"color": col},
                       capprops={"color": col})
            ax.scatter(i + off + rng.uniform(-0.05, 0.05, len(vals)), vals, s=10, color=col, zorder=3,
                       label=("power-mean" if key == "power_mean" else "additive") if i == 0 else None)
    ax.axhline(1, color="k", lw=0.8)
    ax.set_xticks(range(3))
    ax.set_xticklabels([m[1] for m in modes])
    ax.set_ylabel("predicted / Monte Carlo logical error")
    ax.set_title("Decoder-weight ablation (distance 3, s = 0.8 %)")
    ax.legend(frameon=False, fontsize=8)
    ax = axes[1]
    g = mm.groupby("kappa")
    kap = np.array(sorted(mm.kappa.unique()))
    for key, col, lab in (("pred_additive", DARK, "additive"), ("pred_power_mean", ACCENT, "power-mean")):
        r = (mm[key] / mm["eps_all"]).groupby(mm.kappa)
        ax.errorbar(kap, r.mean().loc[kap], r.std().loc[kap], color=col, marker="o", capsize=3, label=lab)
    ax.set_xscale("log", base=2)
    ax.axhline(1, color="k", lw=0.8)
    ax2 = ax.twinx()
    ax2.plot(kap, g["eps_all"].mean().loc[kap], color=MID, ls="--", marker="s")
    ax2.set_ylim(0, 0.02)
    ax2.set_ylabel("Monte Carlo logical error (dashed)", color=MID)
    ax2.spines["right"].set_visible(True)
    ax.set_xlabel("decoder's assumed readout error / true readout error")
    ax.set_ylabel("predicted / Monte Carlo logical error")
    ax.set_title("Decoder mis-calibration sweep")
    ax.legend(frameon=False, fontsize=8, loc="center left")
    _save(fig, path)


def fig_threshold(acc: pd.DataFrame, thr: pd.DataFrame, path: Path, crossing: float | None = None) -> None:
    """Logical error versus noise scale for the two distances (threshold behaviour)."""
    both = pd.concat([acc[["distance", "scale", "seed", "eps_all"]], thr[["distance", "scale", "seed", "eps_all"]]])
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for d, col, mk in ((3, MID, "s"), (5, ACCENT, "o")):
        g = both[both.distance == d].groupby("scale")["eps_all"].agg(["mean", "std"]).reset_index()
        ax.errorbar(g["scale"], g["mean"], g["std"].fillna(0), color=col, marker=mk, capsize=3,
                    label=f"distance {d}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ticks = [0.002, 0.004, 0.008, 0.016, 0.03]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{100 * t:g}%" for t in ticks])
    ax.minorticks_off()
    if crossing:
        ax.axvline(crossing, color=DARK, ls=":", lw=1)
        ax.text(crossing * 1.03, 0.0012, f"curves cross at s = {100 * crossing:.1f}%", fontsize=8)
    ax.set_xlabel("noise scale s (gate error per location)")
    ax.set_ylabel("logical error probability (mean +- std over seeds)")
    ax.set_title("Distance-3 and distance-5 memories")
    ax.legend(frameon=False)
    _save(fig, path)


def make_all_figures(tables: dict[str, pd.DataFrame], summary: dict, out: Path) -> None:
    """Write every figure into ``out``."""
    _style()
    fig_accuracy_violin(tables["budget_accuracy"], out / "fig1_budget_accuracy.png")
    fig_validity_map(tables["validity_map"], out / "fig2_validity_map.png")
    fig_pareto(tables["shot_pareto"], out / "fig3_cost_accuracy_pareto.png")
    fig_ablation(tables["decoder_ablation"], tables["decoder_mismatch"], out / "fig4_decoder_ablation.png")
    fig_threshold(tables["budget_accuracy"], tables["threshold_scan"], out / "fig5_logical_error_scaling.png",
                  summary.get("threshold_crossing_scale"))
