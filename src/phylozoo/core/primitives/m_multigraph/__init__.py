"""
Mixed multi-graph module.

This module provides the MixedMultiGraph class and related functions for working
with mixed multi-graphs. A mixed multi-graph contains both directed and undirected
edges and allows multiple edges between vertices. It is the underlying graph
structure for semi-directed phylogenetic networks, where tree edges are undirected
and reticulation arcs are directed.
"""

from .base import MixedMultiGraph
from .isomorphism import is_isomorphic
from . import io

__all__ = [
    "MixedMultiGraph",
    "is_isomorphic",
    "io",
]
