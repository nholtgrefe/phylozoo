"""
Core module for PhyloZoo.

This module contains fundamental data structures and classes used throughout
the package.
"""

from .network import DirectedPhyNetwork, SemiDirectedPhyNetwork
from .split import Split, SplitSystem, WeightedSplitSystem
from .distance import DistanceMatrix
from .sequence import MSA
from .quartet import Quartet, QuartetProfile, QuartetProfileSet

__all__ = [
    # Networks
    "DirectedPhyNetwork",
    "SemiDirectedPhyNetwork",
    # Split structures
    "Split",
    "SplitSystem",
    "WeightedSplitSystem",
    # Distance
    "DistanceMatrix",
    # Sequence
    "MSA",
    # Quartet
    "Quartet",
    "QuartetProfile",
    "QuartetProfileSet",
]

