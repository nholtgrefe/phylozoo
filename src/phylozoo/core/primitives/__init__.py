"""
Fundamental data structures used throughout PhyloZoo.

This submodule provides the underlying graph and set-theoretic structures that
support the higher-level phylogenetic network classes: Partition, CircularOrdering,
DirectedMultiGraph, and MixedMultiGraph.
"""

from ..split import Split, SplitSystem
from .partition import Partition
from .circular_ordering import CircularOrdering, CircularSetOrdering
from .d_multigraph import DirectedMultiGraph
from .m_multigraph import MixedMultiGraph

__all__ = [
    "Split",
    "SplitSystem",
    "Partition",
    "CircularOrdering",
    "CircularSetOrdering",
    "MixedMultiGraph",
    "DirectedMultiGraph",
]
