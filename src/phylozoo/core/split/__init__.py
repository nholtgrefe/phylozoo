"""
Split module for PhyloZoo.

This module provides classes for working with phylogenetic splits and split systems.
A split is a bipartition of a set of taxa, representing a division of the taxa into
two non-empty subsets. Split systems are collections of splits that can represent
phylogenetic trees or networks. The public API (Split, SplitSystem, WeightedSplitSystem,
to_weightedsplitsystem, is_compatible, is_subsplit, and the algorithms, classifications,
io submodules) is re-exported here; the implementation is split across the base,
splitsystem, weighted_splitsystem, algorithms, classifications, and io submodules.
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


