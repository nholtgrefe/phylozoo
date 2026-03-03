"""
Directed multi-graph module.

This module provides the DirectedMultiGraph class and related functions for working
with directed multi-graphs. A directed multi-graph allows multiple directed edges
between the same pair of vertices. It is the underlying graph structure for
directed phylogenetic networks, where multiple edges can represent parallel arcs.
"""

from .base import DirectedMultiGraph
from .isomorphism import is_isomorphic
from . import io

__all__ = [
    "DirectedMultiGraph",
    "is_isomorphic",
    "io",
]
