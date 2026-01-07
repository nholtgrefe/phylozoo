"""
Split module for PhyloZoo.

This submodule bundles classes for working with splits and split systems.
"""

from .base import Split, is_compatible, is_subsplit
from .splitsystem import SplitSystem
from .weighted_splitsystem import WeightedSplitSystem, to_weightedsplitsystem
from . import algorithms, classifications, io

__all__ = [
    "Split",
    "SplitSystem",
    "WeightedSplitSystem",
    "to_weightedsplitsystem",
    "is_compatible",
    "is_subsplit",
    "algorithms",
    "classifications",
    "io",
]


