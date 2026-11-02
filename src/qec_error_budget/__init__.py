"""Error budgets for a circuit-level repetition-code memory (simulation study)."""

from .budget import (
    additive_budget,
    classical_majority_budget,
    min_fault_weight,
    power_mean_budget,
    rescale_additive,
    rescale_power_mean,
)
from .circuit import SOURCES, FaultTable, build_fault_table
from .decoder import MatchingDecoder
from .dem import DetectorErrorModel, NoiseParams
from .simulate import get_table, logical_failures

__all__ = [
    "SOURCES",
    "FaultTable",
    "build_fault_table",
    "DetectorErrorModel",
    "NoiseParams",
    "MatchingDecoder",
    "get_table",
    "logical_failures",
    "additive_budget",
    "power_mean_budget",
    "rescale_additive",
    "rescale_power_mean",
    "classical_majority_budget",
    "min_fault_weight",
]

__version__ = "0.1.0"
