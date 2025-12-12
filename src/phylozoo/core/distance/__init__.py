"""
Distance module.

This module provides classes and functions for working with distance matrices.
"""

from .base import DistanceMatrix
from . import classifications, operations, io

__all__ = [
    "DistanceMatrix",
    "classifications",
    "operations",
    "io",
]
