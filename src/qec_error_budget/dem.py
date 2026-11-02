"""Detector error model: merged independent error mechanisms and sampling."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .circuit import FaultTable


@dataclass(frozen=True)
class NoiseParams:
    """Per-location bit-flip probabilities of the three noise sources."""

    gate: float = 0.0
    readout: float = 0.0
    idle: float = 0.0

    def as_tuple(self) -> tuple[float, float, float]:
        """Return ``(gate, readout, idle)``."""
        return (self.gate, self.readout, self.idle)

    def only(self, source: int) -> "NoiseParams":
        """Keep a single source (0 gate, 1 readout, 2 idle) and zero the others."""
        v = [0.0, 0.0, 0.0]
        v[source] = self.as_tuple()[source]
        return NoiseParams(*v)

    def scaled(self, factors: tuple[float, float, float]) -> "NoiseParams":
        """Multiply each source by its own factor."""
        a, b, c = (p * f for p, f in zip(self.as_tuple(), factors))
        return NoiseParams(a, b, c)


class DetectorErrorModel:
    """Independent error mechanisms, each flipping a fixed set of detectors.

    Fault locations with identical effect are merged with the exclusive-or rule
    ``p <- p + q - 2 p q``.

    Args:
        table: Fault table of the circuit.
        noise: Source error probabilities.
    """

    def __init__(self, table: FaultTable, noise: NoiseParams) -> None:
        self.table = table
        self.noise = noise
        probs = noise.as_tuple()
        merged: dict[tuple[tuple[int, ...], bool], float] = {}
        for src, det, obs in zip(table.sources, table.detectors, table.flips_observable):
            q = probs[src]
            if q <= 0.0 or (not det and not obs):
                continue
            key = (det, obs)
            p = merged.get(key, 0.0)
            merged[key] = p + q - 2.0 * p * q
        self.mechanisms: list[tuple[tuple[int, ...], bool, float]] = [
            (d, o, p) for (d, o), p in merged.items()
        ]
        n = len(self.mechanisms)
        self._h = np.zeros((n, table.n_detectors), dtype=np.float32)
        self._obs = np.zeros(n, dtype=np.float32)
        self._p = np.zeros(n)
        for k, (d, o, p) in enumerate(self.mechanisms):
            self._h[k, list(d)] = 1.0
            self._obs[k] = float(o)
            self._p[k] = p

    def sample(self, shots: int, rng: np.random.Generator, chunk: int = 20000):
        """Draw ``shots`` samples; returns (detector events, observable flips)."""
        dets, obs = [], []
        done = 0
        while done < shots:
            m = min(chunk, shots - done)
            f = (rng.random((m, len(self._p))) < self._p).astype(np.float32)
            dets.append((f @ self._h).astype(np.int64) & 1)
            obs.append((f @ self._obs).astype(np.int64) & 1)
            done += m
        return np.concatenate(dets).astype(bool), np.concatenate(obs).astype(bool)
