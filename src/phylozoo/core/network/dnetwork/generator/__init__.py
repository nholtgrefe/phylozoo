"""
Generator module for level-k phylogenetic network generators.

This module provides classes and functions for working with level-k generators
of directed phylogenetic networks, as well as utilities for attaching leaves
to generators to obtain full directed phylogenetic networks.
"""

from .attachment import attach_leaves_to_generator
from .base import DirectedGenerator, generators_from_network
from .construction import all_level_k_generators
from .side import DirEdgeSide, EdgeSide, HybridSide, Level0NodeSide, NodeSide, Side

__all__ = [
    "DirectedGenerator",
    "generators_from_network",
    "Side",
    "NodeSide",
    "Level0NodeSide",
    "HybridSide",
    "EdgeSide",
    "DirEdgeSide",
    "attach_leaves_to_generator",
    "all_level_k_generators",
]
