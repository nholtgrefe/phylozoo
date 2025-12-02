"""
Structure module for PhyloZoo.

This submodule contains fundamental data structures used throughout the package.
"""

from .splits import Split
from .splitsystems import SplitSystem
from .partition import Partition
from .circular import CircularOrdering, CircularSetOrdering

__all__ = [
    "Split",
    "SplitSystem",
    "Partition",
    "CircularOrdering",
    "CircularSetOrdering",
]
