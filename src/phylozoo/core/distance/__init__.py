"""
Distance module.

This module provides classes and functions for working with distance matrices. A
distance matrix represents pairwise distances between a set of labeled items, where
distances satisfy properties such as symmetry and non-negativity. The public API
(DistanceMatrix and the classifications, io submodules) is re-exported here; the
implementation is split across the base, classifications, and io submodules.
"""

from .base import DistanceMatrix
from .decomposition import isolation_index, split_decomposition
from . import classifications, decomposition, io

__all__ = [
    "DistanceMatrix",
    "isolation_index",
    "split_decomposition",
    "classifications",
    "decomposition",
    "io",
]
