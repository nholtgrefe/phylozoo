"""
Directed network module.

This module provides the DirectedPhyNetwork class and related functions for working
with directed phylogenetic networks. Directed phylogenetic networks are rooted
networks where all edges are directed. Internal nodes are either tree nodes
(in-degree 1, out-degree >= 2) or hybrid nodes (in-degree >= 2, out-degree 1).
The public API is re-exported here; the implementation is split across the base,
features, classifications, transformations, derivations, conversions, isomorphism,
and io submodules.
"""

from .base import DirectedPhyNetwork
from . import classifications, features, transformations, derivations, io, conversions, isomorphism

__all__ = [
    "DirectedPhyNetwork",
    "classifications",
    "features",
    "transformations",
    "derivations",
    "io",
    "conversions",
    "isomorphism",
]

