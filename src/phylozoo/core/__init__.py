"""
Core module for PhyloZoo.

This module contains fundamental data structures and classes used throughout
the package.
"""

from .network import DirectedPhyNetwork, SemiDirectedPhyNetwork
from .primitives import Partition, CircularOrdering, CircularSetOrdering
from .split import Split, SplitSystem
from .distance import DistanceMatrix
from .sequence import MSA

__all__ = [
    # Networks
    "DirectedPhyNetwork",
    "SemiDirectedPhyNetwork",
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

