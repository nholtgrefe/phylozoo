"""
Split module for PhyloZoo.

This submodule bundles classes for working with splits and split systems.
"""

from .base import Split
from .splitsystem import SplitSystem
from .weighted_splitsystem import WeightedSplitSystem, to_weightedsplitsystem

__all__ = ["Split", "SplitSystem", "WeightedSplitSystem", "to_weightedsplitsystem"]


