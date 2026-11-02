"""Fast tests of the fault table, decoder and budget predictors."""

from __future__ import annotations

import numpy as np
import pytest

from qec_error_budget import (
    DetectorErrorModel,
    MatchingDecoder,
    NoiseParams,
    additive_budget,
    build_fault_table,
    classical_majority_budget,
    logical_failures,
    min_fault_weight,
    power_mean_budget,
    rescale_additive,
    rescale_power_mean,
)


@pytest.mark.parametrize("distance", [3, 5])
def test_fault_table_is_graph_like(distance: int) -> None:
    table = build_fault_table(distance)
    assert table.n_detectors == (distance - 1) * (distance + 1)
    assert max(len(d) for d in table.detectors) <= 2
    # no fault flips the logical observable without being detected
    assert not any(o and not d for d, o in zip(table.detectors, table.flips_observable))


def test_invalid_distance() -> None:
    with pytest.raises(ValueError):
        build_fault_table(4)


@pytest.mark.parametrize("distance", [3, 5])
def test_decoder_corrects_every_single_fault(distance: int) -> None:
    table = build_fault_table(distance)
    model = DetectorErrorModel(table, NoiseParams(0.01, 0.02, 0.01))
    decoder = MatchingDecoder(model)
    for det, obs, _ in model.mechanisms:
        row = np.zeros((1, table.n_detectors), dtype=bool)
        row[0, list(det)] = True
        assert bool(decoder.decode(row)[0]) == obs


def test_sampler_marginals() -> None:
    table = build_fault_table(3)
    model = DetectorErrorModel(table, NoiseParams(0.0, 0.0, 0.05))
    events, _ = model.sample(20000, np.random.default_rng(0))
    assert events.shape == (20000, table.n_detectors)
    assert 0.0 < events.mean() < 0.2


def test_noiseless_never_fails() -> None:
    assert logical_failures(3, NoiseParams(), 100, np.random.default_rng(0)) == 0


def test_larger_distance_is_better_below_threshold() -> None:
    noise = NoiseParams(0.004, 0.008, 0.004)
    rng = np.random.default_rng(3)
    p3 = logical_failures(3, noise, 20000, rng) / 20000
    p5 = logical_failures(5, noise, 20000, rng) / 20000
    assert p5 < p3


def test_budget_identities() -> None:
    eps = [1e-3, 4e-4, 2e-4]
    assert additive_budget(eps) == pytest.approx(1.6e-3)
    assert power_mean_budget(eps, 1.0) == pytest.approx(additive_budget(eps))
    for t in (2, 3, 4):
        assert power_mean_budget(eps, t) >= additive_budget(eps)
        assert rescale_power_mean(eps, [1, 1, 1], t) == pytest.approx(power_mean_budget(eps, t))
        assert rescale_additive(eps, [1, 1, 1], t) == pytest.approx(additive_budget(eps))
    # halving every source divides a power-mean budget by 2**t
    assert rescale_power_mean(eps, [0.5] * 3, 3) == pytest.approx(power_mean_budget(eps, 3) / 8)
    assert min_fault_weight(5) == 3


def test_classical_baseline_monotone() -> None:
    lo = classical_majority_budget(3, 3, 0.001, 0.002, 0.001)
    hi = classical_majority_budget(3, 3, 0.002, 0.004, 0.002)
    assert 0.0 < lo < hi < 0.5
