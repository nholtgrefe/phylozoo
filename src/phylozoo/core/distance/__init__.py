"""
Distance module.

This module provides classes and functions for working with distance matrices. A
distance matrix represents pairwise distances between a set of labeled items, where
distances satisfy properties such as symmetry and non-negativity. The public API
(DistanceMatrix and the classifications, operations, io submodules) is re-exported
here; the implementation is split across the base, classifications, operations,
and io submodules.
"""

from .base import DistanceMatrix
from . import classifications, operations, io

__all__ = [
    "DistanceMatrix",
    "classifications",
    "operations",
    "io",
]
