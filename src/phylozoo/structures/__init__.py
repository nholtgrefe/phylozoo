"""
Structures module for PhyloZoo.

This submodule contains fundamental data structures used throughout the package:
- Splits: 2-partitions of sets
- Partitions: General partitions of sets
- Circular orderings: Circular arrangements of elements

These are fundamental concepts that other submodules build upon.
"""

# Import fundamental classes for convenient access
from .splits import Split, SplitSet, SplitSystem, QuartetSplitSet
from .partition import Partition
from .circular import CircularOrdering, CircularSetOrdering

__all__ = [
    # Splits
    "Split",
    "SplitSet",
    "SplitSystem",
    "QuartetSplitSet",
    # Partitions
    "Partition",
    # Circular orderings
    "CircularOrdering",
    "CircularSetOrdering",
]

