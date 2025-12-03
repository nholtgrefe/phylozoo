"""
Core module for PhyloZoo.

This module contains fundamental data structures and classes used throughout
the package.
"""

from .network import DirectedNetwork, SemiDirectedNetwork, random_semi_directed_network
from .structure import Partition, CircularOrdering, CircularSetOrdering, MixedMultiGraph, DirectedMultiGraph
from .split import Split, SplitSystem
from .distance import DistanceMatrix
from .sequence import MSA

__all__ = [
    # Networks
    "DirectedNetwork",
    "SemiDirectedNetwork",
    "random_semi_directed_network",
    "MixedMultiGraph",
    "DirectedMultiGraph",
    # Split structures
    "Split",
    "SplitSystem",
    # Other structures
    "Partition",
    "CircularOrdering",
    "CircularSetOrdering",
    # Distance
    "DistanceMatrix",
    # Sequence
    "MSA",
]

