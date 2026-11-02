"""Error-budget predictors for the combined logical error probability.

All predictors take the logical failure probabilities ``eps_s`` measured with
*one* noise source switched on at a time and return a prediction for all sources
on together. ``t = (distance + 1) / 2`` is the minimum number of faults that can
defeat a distance-``distance`` repetition code.
"""

from __future__ import annotations

from math import comb
from typing import Sequence

import numpy as np


def min_fault_weight(distance: int) -> int:
    """Smallest number of faults that causes a logical error: ``(d + 1) / 2``."""
    return (distance + 1) // 2


def additive_budget(eps: Sequence[float]) -> float:
    """One-at-a-time budget: the sum of single-source logical error rates."""
    return float(np.sum(eps))


def power_mean_budget(eps: Sequence[float], t: float) -> float:
    """Power-mean budget ``(sum_s eps_s ** (1 / t)) ** t``.

    It is exact when logical failures are caused by ``t`` faults drawn from the
    sources in proportion to their effective rates ``w_s = eps_s ** (1 / t)``,
    and equals the additive budget for ``t = 1``.
    """
    w = np.power(np.asarray(eps, dtype=float), 1.0 / t)
    return float(np.sum(w) ** t)


def rescale_additive(eps: Sequence[float], factors: Sequence[float], t: float) -> float:
    """Additive budget after multiplying source rates by ``factors``."""
    return float(np.sum(np.asarray(eps) * np.power(np.asarray(factors, float), t)))


def rescale_power_mean(eps: Sequence[float], factors: Sequence[float], t: float) -> float:
    """Power-mean budget after multiplying source rates by ``factors``."""
    w = np.power(np.asarray(eps, dtype=float), 1.0 / t) * np.asarray(factors, float)
    return float(np.sum(w) ** t)


def classical_majority_budget(distance: int, rounds: int, gate: float, readout: float, idle: float) -> float:
    """Lumped classical baseline: majority vote over ``distance`` bits.

    Each data bit accumulates ``rounds * (idle + 2 * gate)`` flips (two CNOT
    layers, first-order) plus one readout flip; syndrome-measurement noise is
    ignored. The block fails when at least ``(distance + 1) / 2`` bits flip.
    """
    q = min(rounds * (idle + 2.0 * gate) + readout, 0.5)
    t = min_fault_weight(distance)
    return float(sum(comb(distance, k) * q**k * (1 - q) ** (distance - k) for k in range(t, distance + 1)))
