"""Experiment drivers. Each function returns a :class:`pandas.DataFrame` of raw rows."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .budget import (
    additive_budget,
    classical_majority_budget,
    min_fault_weight,
    power_mean_budget,
    rescale_additive,
    rescale_power_mean,
)
from .dem import NoiseParams
from .simulate import logical_failures

#: Default noise mix: readout is twice as likely to flip as a CNOT output or idle qubit.
MIX: tuple[float, float, float] = (1.0, 2.0, 1.0)


def mix_noise(scale: float, mix: tuple[float, float, float] = MIX) -> NoiseParams:
    """Noise with source probabilities ``scale * mix``."""
    return NoiseParams(*(scale * m for m in mix))


@dataclass
class Config:
    """Sizes of all experiments; :meth:`quick` shrinks them for a smoke run."""

    seeds: int = 5
    map_seeds: int = 3
    acc_d3: tuple[float, ...] = (0.002, 0.004, 0.008, 0.016)
    acc_d5: tuple[float, ...] = (0.012, 0.016, 0.02)
    shots_d3: int = 100_000
    shots_d5: int = 60_000
    map_scales: tuple[float, ...] = (0.002, 0.004, 0.008, 0.012, 0.016, 0.024)
    map_ratios: tuple[float, ...] = (0.25, 0.5, 1.0, 2.0, 4.0, 8.0)
    map_ref_shots: int = 600_000
    map_shots: int = 100_000
    abl_shots: tuple[int, int] = (100_000, 60_000)
    kappas: tuple[float, ...] = (0.25, 0.5, 1.0, 2.0, 4.0)
    pareto_shots: tuple[int, ...] = (2_000, 5_000, 20_000, 80_000, 200_000)
    pareto_reps: int = 10
    pareto_ref_shots: int = 1_000_000
    thr_scales: tuple[float, ...] = (0.02, 0.03, 0.04)
    thr_shots: int = 10_000
    thr_seeds: int = 3
    extra: dict = field(default_factory=dict)

    @classmethod
    def quick(cls) -> "Config":
        """Tiny configuration for CI (well under one minute)."""
        return cls(
            seeds=2, map_seeds=2, acc_d3=(0.004, 0.016), acc_d5=(0.02,),
            shots_d3=4_000, shots_d5=2_000,
            map_scales=(0.004, 0.016), map_ratios=(0.5, 4.0),
            map_ref_shots=20_000, map_shots=4_000, abl_shots=(4_000, 2_000),
            kappas=(0.5, 2.0), pareto_shots=(1_000, 4_000), pareto_reps=3,
            pareto_ref_shots=20_000, thr_scales=(0.03,), thr_shots=1_000, thr_seeds=2,
        )


def _rng(*tags: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence(list(tags)))


def _single_and_all(
    d: int, noise: NoiseParams, shots: int, rng: np.random.Generator,
    decoder: str = "per_run",
) -> tuple[list[int], int]:
    """Failure counts of the three single-source runs and the all-sources run.

    ``decoder``: ``per_run`` rebuilds the matching weights from each run's own
    noise; ``fixed`` uses the full-noise model for all runs; ``uniform`` uses the
    full-noise graph with unit weights.
    """
    kw: dict = {}
    if decoder in ("fixed", "uniform"):
        kw = {"decoder_noise": noise, "uniform_weights": decoder == "uniform"}
    singles = [logical_failures(d, noise.only(k), shots, rng, **kw) for k in range(3)]
    return singles, logical_failures(d, noise, shots, rng, **kw)


def _predictions(d: int, eps: list[float], noise: NoiseParams) -> dict[str, float]:
    t = min_fault_weight(d)
    return {
        "additive": additive_budget(eps),
        "power_mean": power_mean_budget(eps, t),
        "classical_majority": classical_majority_budget(d, d, *noise.as_tuple()),
    }


def budget_accuracy(cfg: Config) -> pd.DataFrame:
    """Calibrate on single sources at the same point; compare with the full run."""
    rows = []
    for d, scales, shots in ((3, cfg.acc_d3, cfg.shots_d3), (5, cfg.acc_d5, cfg.shots_d5)):
        for s in scales:
            for seed in range(cfg.seeds):
                noise = mix_noise(s)
                rng = _rng(1, d, int(s * 1e4), seed)
                singles, full = _single_and_all(d, noise, shots, rng)
                eps = [c / shots for c in singles]
                row = {"distance": d, "scale": s, "seed": seed, "shots": shots,
                       "eps_gate": eps[0], "eps_readout": eps[1], "eps_idle": eps[2],
                       "eps_all": full / shots, "fail_all": full}
                row.update({f"pred_{k}": v for k, v in _predictions(d, eps, noise).items()})
                rows.append(row)
    return pd.DataFrame(rows)


def validity_map(cfg: Config) -> pd.DataFrame:
    """Extrapolate one calibration (d=3) over noise scale and readout ratio."""
    d, t = 3, min_fault_weight(3)
    s_ref, rho_ref = 0.004, MIX[1]
    rows = []
    for seed in range(cfg.map_seeds):
        rng = _rng(2, seed)
        ref = mix_noise(s_ref)
        eps = [logical_failures(d, ref.only(k), cfg.map_ref_shots, rng) / cfg.map_ref_shots
               for k in range(3)]
        for s in cfg.map_scales:
            for rho in cfg.map_ratios:
                noise = NoiseParams(s, rho * s, s)
                fac = (s / s_ref, rho * s / (rho_ref * s_ref), s / s_ref)
                fail = logical_failures(d, noise, cfg.map_shots, rng)
                rows.append({
                    "seed": seed, "scale": s, "ratio": rho, "shots": cfg.map_shots,
                    "eps_all": fail / cfg.map_shots,
                    "pred_additive": rescale_additive(eps, fac, t),
                    "pred_power_mean": rescale_power_mean(eps, fac, t),
                })
    return pd.DataFrame(rows)


def decoder_ablation(cfg: Config) -> pd.DataFrame:
    """Per-run matched decoder vs one fixed full-model decoder vs unit weights."""
    rows = []
    for (d, s), shots in zip(((3, 0.008), (5, 0.016)), cfg.abl_shots):
        for mode in ("per_run", "fixed", "uniform"):
            for seed in range(cfg.seeds):
                noise = mix_noise(s)
                rng = _rng(3, d, seed, ("per_run", "fixed", "uniform").index(mode))
                singles, full = _single_and_all(d, noise, shots, rng, decoder=mode)
                eps = [c / shots for c in singles]
                row = {"distance": d, "scale": s, "decoder": mode, "seed": seed,
                       "eps_all": full / shots}
                row.update({f"pred_{k}": v for k, v in _predictions(d, eps, noise).items()})
                rows.append(row)
    return pd.DataFrame(rows)


def exponent_ablation(acc: pd.DataFrame) -> pd.DataFrame:
    """Power-mean exponent: 1 (additive), theory ``t``, and a pooled fit per distance."""
    rows = []
    for d, grp in acc.groupby("distance"):
        pooled = grp.groupby("scale")[["eps_gate", "eps_readout", "eps_idle", "eps_all"]].mean()
        grid = np.linspace(1.0, 2.0 * min_fault_weight(d), 61)
        best_t, best_err = 1.0, np.inf
        for t in grid:
            err = np.mean([
                (np.log10(power_mean_budget(r[:3].values, t)) - np.log10(r["eps_all"])) ** 2
                for _, r in pooled.iterrows()
            ])
            if err < best_err:
                best_t, best_err = float(t), float(err)
        for label, t in (("t=1 (additive)", 1.0), ("t=(d+1)/2", float(min_fault_weight(d))),
                         ("t fitted", best_t)):
            errs = [np.log10(power_mean_budget(r[:3].values, t)) - np.log10(r["eps_all"])
                    for _, r in pooled.iterrows()]
            rows.append({"distance": d, "exponent": label, "t": t,
                         "rms_log10_error": float(np.sqrt(np.mean(np.square(errs))))})
    return pd.DataFrame(rows)


def decoder_mismatch(cfg: Config) -> pd.DataFrame:
    """A decoder that assumes readout noise ``kappa`` times too large (d=3)."""
    d, s = 3, 0.008
    shots = cfg.abl_shots[0]
    rows = []
    for kappa in cfg.kappas:
        for seed in range(cfg.seeds):
            noise = mix_noise(s)
            assumed = noise.scaled((1.0, kappa, 1.0))
            rng = _rng(4, int(kappa * 100), seed)
            kw = {"decoder_noise": assumed}
            singles = [logical_failures(d, noise.only(k), shots, rng, **kw) for k in range(3)]
            full = logical_failures(d, noise, shots, rng, **kw)
            eps = [c / shots for c in singles]
            rows.append({"kappa": kappa, "seed": seed, "eps_all": full / shots,
                         "pred_additive": additive_budget(eps),
                         "pred_power_mean": power_mean_budget(eps, min_fault_weight(d))})
    return pd.DataFrame(rows)


def shot_pareto(cfg: Config) -> pd.DataFrame:
    """Prediction error of the combined logical error versus Monte Carlo cost."""
    d, s = 3, 0.008
    t = min_fault_weight(d)
    noise = mix_noise(s)
    ref_rng = _rng(5, 0)
    ref = logical_failures(d, noise, cfg.pareto_ref_shots, ref_rng) / cfg.pareto_ref_shots
    rows = []
    for n in cfg.pareto_shots:
        for rep in range(cfg.pareto_reps):
            rng = _rng(5, n, rep)
            eps = [logical_failures(d, noise.only(k), n, rng) / n for k in range(3)]
            direct = logical_failures(d, noise, 3 * n, rng) / (3 * n)
            for method, pred in (("direct Monte Carlo", direct),
                                 ("power-mean budget", power_mean_budget(eps, t)),
                                 ("additive budget", additive_budget(eps))):
                err = np.log10(pred / ref) if pred > 0 else np.nan
                rows.append({"method": method, "total_shots": 3 * n, "rep": rep,
                             "pred": pred, "ref": ref, "log10_error": err})
    return pd.DataFrame(rows)


def threshold_scan(cfg: Config) -> pd.DataFrame:
    """Logical error of d=3 and d=5 near the threshold of the noise mix."""
    rows = []
    for d in (3, 5):
        for s in cfg.thr_scales:
            for seed in range(cfg.thr_seeds):
                rng = _rng(6, d, int(s * 1e4), seed)
                f = logical_failures(d, mix_noise(s), cfg.thr_shots, rng)
                rows.append({"distance": d, "scale": s, "seed": seed, "eps_all": f / cfg.thr_shots})
    return pd.DataFrame(rows)
