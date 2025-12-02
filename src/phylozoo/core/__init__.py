"""
Core module for PhyloZoo.

This module contains fundamental data structures and classes used throughout
the package.
"""

from .network import DirectedNetwork, SemiDirectedNetwork, random_semi_directed_network, MixedGraph, MultiMixedGraph
from .structure import Split, SplitSystem, Partition, CircularOrdering, CircularSetOrdering
from .distance import DistanceMatrix
from .sequence import MSA

__all__ = [
    # Networks
    "DirectedNetwork",
    "SemiDirectedNetwork",
    "random_semi_directed_network",
    "MixedGraph",
    "MultiMixedGraph",
    # Structures
    "Split",
    "SplitSystem",
    "Partition",
    "CircularOrdering",
    "CircularSetOrdering",
    # Distance
    "DistanceMatrix",
    # Sequence
    "MSA",
]

