"""Minimum-weight matching decoder on the detector graph of a DEM."""

from __future__ import annotations

import math

import networkx as nx
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

from .dem import DetectorErrorModel

_DP_LIMIT = 14  # exact subset DP up to this many defects, networkx matching beyond


class MatchingDecoder:
    """Decode by pairing detection events (or sending them to the boundary).

    Edge weights are ``ln((1 - p) / p)`` of the decoder's *assumed* model, or 1
    for every edge when ``uniform=True``. Weights may therefore differ from the
    noise that generated the data, which is how mis-calibration is studied.

    Args:
        model: Detector error model assumed by the decoder (graph-like).
        uniform: Use unit weights instead of log-likelihood weights.
    """

    def __init__(self, model: DetectorErrorModel, uniform: bool = False) -> None:
        n = model.table.n_detectors
        self.n_detectors = n
        self.boundary = n
        best: dict[tuple[int, int], tuple[float, bool]] = {}
        for det, obs, p in model.mechanisms:
            if len(det) == 1:
                key = (det[0], n)
            elif len(det) == 2:
                key = (det[0], det[1])
            else:
                raise ValueError(f"mechanism with {len(det)} detectors is not graph-like")
            w = 1.0 if uniform else math.log((1.0 - p) / p)
            w = max(w, 1e-9)
            if key not in best or w < best[key][0]:
                best[key] = (w, obs)
        rows, cols, vals = [], [], []
        for (u, v), (w, _) in best.items():
            rows += [u, v]
            cols += [v, u]
            vals += [w, w]
        graph = csr_matrix((vals, (rows, cols)), shape=(n + 1, n + 1))
        dist, pred = dijkstra(graph, directed=False, return_predecessors=True)
        self._dist = dist
        parity = np.zeros((n + 1, n + 1), dtype=np.int64)
        edge_obs = {}
        for (u, v), (_, o) in best.items():
            edge_obs[(u, v)] = edge_obs[(v, u)] = int(o)
        for s in range(n + 1):
            for t in range(n + 1):
                x, par = t, 0
                while x != s and pred[s, x] >= 0:
                    par ^= edge_obs[(int(pred[s, x]), x)]
                    x = int(pred[s, x])
                parity[s, t] = par
        self._parity = parity
        self._cache: dict[bytes, int] = {}

    def _solve(self, defects: list[int]) -> int:
        """Return the predicted observable flip for one defect set."""
        n = len(defects)
        b = self.boundary
        d, par = self._dist, self._parity
        if n > _DP_LIMIT:
            return self._solve_large(defects)
        memo: dict[int, tuple[float, int]] = {0: (0.0, 0)}

        def f(mask: int) -> tuple[float, int]:
            if mask in memo:
                return memo[mask]
            i = (mask & -mask).bit_length() - 1
            rest = mask & ~(1 << i)
            best = d[defects[i], b] + 0.0
            c0, p0 = f(rest)
            res = (best + c0, p0 ^ int(par[defects[i], b]))
            m = rest
            while m:
                j = (m & -m).bit_length() - 1
                m &= m - 1
                cj, pj = f(rest & ~(1 << j))
                cost = d[defects[i], defects[j]] + cj
                if cost < res[0] - 1e-12:
                    res = (cost, pj ^ int(par[defects[i], defects[j]]))
            memo[mask] = res
            return res

        return f((1 << n) - 1)[1]

    def _solve_large(self, defects: list[int]) -> int:
        """Matching via networkx for large defect sets (rare)."""
        g = nx.Graph()
        n = len(defects)
        for i in range(n):
            for j in range(i + 1, n):
                g.add_edge(i, j, weight=float(self._dist[defects[i], defects[j]]))
                g.add_edge(n + i, n + j, weight=0.0)
            g.add_edge(i, n + i, weight=float(self._dist[defects[i], self.boundary]))
        match = nx.min_weight_matching(g)
        par = 0
        for u, v in match:
            if u < n and v < n:
                par ^= int(self._parity[defects[u], defects[v]])
            elif u < n <= v and v - n == u:
                par ^= int(self._parity[defects[u], self.boundary])
            elif v < n <= u and u - n == v:
                par ^= int(self._parity[defects[v], self.boundary])
        return par

    def decode(self, events: np.ndarray) -> np.ndarray:
        """Predict logical flips for a ``(shots, n_detectors)`` boolean array."""
        uniq, inverse = np.unique(events, axis=0, return_inverse=True)
        inverse = inverse.reshape(-1)
        out = np.zeros(len(uniq), dtype=bool)
        for k, row in enumerate(uniq):
            key = np.packbits(row).tobytes()
            if key not in self._cache:
                self._cache[key] = self._solve([int(x) for x in np.flatnonzero(row)])
            out[k] = bool(self._cache[key])
        return out[inverse]
