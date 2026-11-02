"""Fault table of a bit-flip repetition-code memory experiment.

The circuit has ``rounds`` syndrome-extraction cycles. Each cycle applies two
layers of CNOT gates (data qubit -> check ancilla), an idle period on the data
qubits while the ancillas are read out, and an ancilla readout. After the last
cycle all data qubits are read out. Every noise location is a single-qubit
X (bit-flip) fault, so fault effects are linear and can be propagated exactly.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

#: Noise sources of the budget: two-qubit gate, readout, data idling.
SOURCES: tuple[str, ...] = ("gate", "readout", "idle")


@dataclass(frozen=True)
class FaultTable:
    """All single-fault locations of a memory experiment and their effects.

    Attributes:
        distance: Number of data qubits (odd).
        rounds: Number of syndrome-extraction cycles.
        n_detectors: Number of detectors, ``(distance - 1) * (rounds + 1)``.
        sources: Source index (into :data:`SOURCES`) of every fault location.
        detectors: Sorted tuple of detectors flipped by each fault location.
        flips_observable: Whether each fault location flips the logical observable.
    """

    distance: int
    rounds: int
    n_detectors: int
    sources: tuple[int, ...]
    detectors: tuple[tuple[int, ...], ...]
    flips_observable: tuple[bool, ...]

    @property
    def n_faults(self) -> int:
        """Number of fault locations."""
        return len(self.sources)


def _ops(distance: int, rounds: int) -> list[tuple]:
    """Build the operation list; ``("F", ...)`` entries are fault markers."""
    n_anc = distance - 1
    ops: list[tuple] = []
    for _ in range(rounds):
        for layer in (0, 1):
            for i in range(n_anc):
                ctrl = i + layer
                ops.append(("CX", ctrl, distance + i))
                ops.append(("F", "X", ctrl, 0))
                ops.append(("F", "X", distance + i, 0))
        for j in range(distance):
            ops.append(("F", "X", j, 2))
        for i in range(n_anc):
            ops.append(("F", "M", i, 1))
            ops.append(("MEAS", i, distance + i))
    for j in range(distance):
        ops.append(("F", "X", j, 1))
    return ops


def _propagate(ops: list[tuple], active: int, distance: int, rounds: int):
    """Run the circuit with only fault marker number ``active`` switched on."""
    n_anc = distance - 1
    flips = np.zeros(2 * distance - 1, dtype=np.uint8)
    meas = np.zeros((rounds, n_anc), dtype=np.uint8)
    t = 0
    marker = -1
    for op in ops:
        kind = op[0]
        if kind == "CX":
            flips[op[2]] ^= flips[op[1]]
        elif kind == "F":
            marker += 1
            if marker == active:
                if op[1] == "X":
                    flips[op[2]] ^= 1
                else:
                    meas[t, op[2]] ^= 1
        else:  # MEAS
            meas[t, op[1]] ^= flips[op[2]]
            flips[op[2]] = 0
            if op[1] == n_anc - 1:
                t += 1
    return meas, flips[:distance]


def build_fault_table(distance: int, rounds: int | None = None) -> FaultTable:
    """Enumerate every fault location and the detectors/observable it flips.

    Args:
        distance: Odd number of data qubits (>= 3).
        rounds: Syndrome-extraction cycles; defaults to ``distance``.

    Returns:
        The fault table. Detector ``t * (distance - 1) + i`` compares check ``i``
        in cycle ``t`` with cycle ``t - 1`` (the last detector layer compares the
        final data readout with the last cycle).
    """
    if distance < 3 or distance % 2 == 0:
        raise ValueError("distance must be an odd integer >= 3")
    rounds = distance if rounds is None else rounds
    ops = _ops(distance, rounds)
    markers = [op for op in ops if op[0] == "F"]
    n_anc = distance - 1
    sources, detectors, observables = [], [], []
    for k, op in enumerate(markers):
        meas, final = _propagate(ops, k, distance, rounds)
        det = np.zeros((rounds + 1, n_anc), dtype=np.uint8)
        det[0] = meas[0]
        det[1:rounds] = meas[1:] ^ meas[:-1]
        det[rounds] = final[:-1] ^ final[1:] ^ meas[-1]
        flat = np.flatnonzero(det.ravel())
        sources.append(op[3])
        detectors.append(tuple(int(x) for x in flat))
        observables.append(bool(final[0]))
    return FaultTable(
        distance=distance,
        rounds=rounds,
        n_detectors=(rounds + 1) * n_anc,
        sources=tuple(sources),
        detectors=tuple(detectors),
        flips_observable=tuple(observables),
    )
