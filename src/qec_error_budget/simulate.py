"""Monte Carlo estimation of logical failure probabilities."""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from .circuit import FaultTable, build_fault_table
from .decoder import MatchingDecoder
from .dem import DetectorErrorModel, NoiseParams


@lru_cache(maxsize=None)
def get_table(distance: int, rounds: int | None = None) -> FaultTable:
    """Cached :func:`build_fault_table`."""
    return build_fault_table(distance, rounds)


def logical_failures(
    distance: int,
    noise: NoiseParams,
    shots: int,
    rng: np.random.Generator,
    decoder_noise: NoiseParams | None = None,
    uniform_weights: bool = False,
    rounds: int | None = None,
) -> int:
    """Count logical failures of the matching decoder over ``shots`` samples.

    Args:
        distance: Repetition-code distance (odd).
        noise: True noise that generates the samples.
        shots: Number of Monte Carlo samples.
        rng: Random generator.
        decoder_noise: Noise assumed by the decoder; defaults to ``noise``. Its
            non-zero sources must include those of ``noise``.
        uniform_weights: Give every decoder edge unit weight.
        rounds: Syndrome cycles; defaults to ``distance``.

    Returns:
        Number of shots whose predicted logical flip differs from the true one.
    """
    table = get_table(distance, rounds)
    truth = DetectorErrorModel(table, noise)
    if not truth.mechanisms:
        return 0
    assumed = DetectorErrorModel(table, decoder_noise if decoder_noise is not None else noise)
    decoder = MatchingDecoder(assumed, uniform=uniform_weights)
    events, obs = truth.sample(shots, rng)
    return int(np.sum(decoder.decode(events) != obs))

