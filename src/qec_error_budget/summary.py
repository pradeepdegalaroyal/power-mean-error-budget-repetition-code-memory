"""Aggregate raw result tables into the numbers quoted in the README."""

from __future__ import annotations

import numpy as np
import pandas as pd

METHODS = ("additive", "power_mean", "classical_majority")
TOL = 0.25  # a budget is "accurate" when within +-25 % of the Monte Carlo value


def _ms(x: pd.Series) -> dict:
    return {"mean": float(np.mean(x)), "std": float(np.std(x, ddof=1)) if len(x) > 1 else 0.0}


def _accuracy(acc: pd.DataFrame) -> list[dict]:
    out = []
    for (d, s), g in acc.groupby(["distance", "scale"]):
        row = {"distance": int(d), "scale": float(s), "eps_all": _ms(g["eps_all"]),
               "n_seeds": int(len(g))}
        for m in METHODS:
            ratio = g[f"pred_{m}"] / g["eps_all"].replace(0, np.nan)
            row[f"ratio_{m}"] = _ms(ratio.dropna())
        out.append(row)
    return out


def _map(vm: pd.DataFrame) -> dict:
    pooled = vm.groupby(["scale", "ratio"]).mean(numeric_only=True).reset_index()
    res = {"n_cells": int(len(pooled))}
    for m in ("additive", "power_mean"):
        err = np.log10(pooled[f"pred_{m}"] / pooled["eps_all"])
        res[m] = {"mean_abs_log10_error": float(np.mean(np.abs(err))),
                  "fraction_within_tol": float(np.mean(np.abs(10**err - 1) <= TOL))}
    both = np.abs(np.log10(pooled["pred_additive"] / pooled["eps_all"])) > np.abs(
        np.log10(pooled["pred_power_mean"] / pooled["eps_all"]))
    res["fraction_power_mean_better"] = float(np.mean(both))
    # the same statistics restricted to scale <= 0.008 (well below threshold)
    low = pooled[pooled["scale"] <= 0.008]
    for m in ("additive", "power_mean"):
        err = np.log10(low[f"pred_{m}"] / low["eps_all"])
        res[m]["fraction_within_tol_low_scale"] = float(np.mean(np.abs(10**err - 1) <= TOL))
    return res


def _crossing(thr: pd.DataFrame) -> float | None:
    """Scale at which the d=5 and d=3 curves cross (log-linear interpolation)."""
    m = thr.groupby(["distance", "scale"])["eps_all"].mean().unstack(0)
    scales = m.index.values
    diff = (m[5] - m[3]).values
    for i in range(len(scales) - 1):
        if diff[i] < 0 <= diff[i + 1]:
            f = -diff[i] / (diff[i + 1] - diff[i])
            return float(scales[i] + f * (scales[i + 1] - scales[i]))
    return None


def summarise(tables: dict[str, pd.DataFrame]) -> dict:
    """Return a JSON-serialisable dictionary of all headline statistics."""
    s: dict = {"accuracy": _accuracy(tables["budget_accuracy"]), "map": _map(tables["validity_map"])}
    acc = tables["budget_accuracy"]
    s["accuracy_overall"] = {}
    for m in METHODS:
        err = np.log10(acc[f"pred_{m}"] / acc["eps_all"].replace(0, np.nan)).replace(
            [np.inf, -np.inf], np.nan).dropna()
        s["accuracy_overall"][m] = {"median_abs_log10_error": float(np.median(np.abs(err))),
                                    "mean_log10_error": float(np.mean(err))}
    abl = tables["decoder_ablation"]
    s["decoder_ablation"] = []
    for (d, mode), g in abl.groupby(["distance", "decoder"]):
        row = {"distance": int(d), "decoder": mode, "eps_all": _ms(g["eps_all"])}
        for m in ("additive", "power_mean"):
            row[f"ratio_{m}"] = _ms(g[f"pred_{m}"] / g["eps_all"].replace(0, np.nan)).copy()
        s["decoder_ablation"].append(row)
    s["exponent_ablation"] = tables["exponent_ablation"].to_dict("records")
    mm = tables["decoder_mismatch"]
    s["decoder_mismatch"] = [
        {"kappa": float(k), "eps_all": _ms(g["eps_all"]),
         "ratio_power_mean": _ms(g["pred_power_mean"] / g["eps_all"].replace(0, np.nan)),
         "ratio_additive": _ms(g["pred_additive"] / g["eps_all"].replace(0, np.nan))}
        for k, g in mm.groupby("kappa")]
    par = tables["shot_pareto"]
    s["shot_pareto"] = [
        {"method": m, "total_shots": int(n),
         "rms_log10_error": float(np.sqrt(np.nanmean(g["log10_error"] ** 2))),
         "n_nan": int(g["log10_error"].isna().sum())}
        for (m, n), g in par.groupby(["method", "total_shots"])]
    s["reference_eps_all_pareto"] = float(par["ref"].iloc[0])
    thr = tables["threshold_scan"]
    s["threshold_crossing_scale"] = _crossing(thr) if thr["scale"].nunique() > 1 else None
    s["threshold_table"] = [
        {"distance": int(d), "scale": float(sc), "eps_all": _ms(g["eps_all"])}
        for (d, sc), g in thr.groupby(["distance", "scale"])]
    return s
