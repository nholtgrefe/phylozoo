"""
Level-k generators for directed phylogenetic networks.

This module provides classes and functions for working with level-k generators
of directed phylogenetic networks. Generators are minimal biconnected components
that represent the core structure of level-k directed phylogenetic networks.
It also provides utilities for attaching leaves to generators to obtain full
directed phylogenetic networks.
"""

from .attachment import attach_leaves_to_generator
from .base import DirectedGenerator, generators_from_network
from .construction import all_level_k_generators
from .side import DirEdgeSide, EdgeSide, HybridSide, IsolatedNodeSide, NodeSide, Side

__all__ = [
    "DirectedGenerator",
    "generators_from_network",
    "Side",
    "NodeSide",
    "IsolatedNodeSide",
    "HybridSide",
    "EdgeSide",
    "DirEdgeSide",
    "attach_leaves_to_generator",
    "all_level_k_generators",
]
