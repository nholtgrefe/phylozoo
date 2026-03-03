"""
Level-k generators for semi-directed phylogenetic networks.

This module provides classes for representing level-k generators of
semi-directed phylogenetic networks. Generators are minimal biconnected
components that represent the core structure of level-k semi-directed
phylogenetic networks.
"""

from .attachment import attach_leaves_to_generator
from .base import SemiDirectedGenerator
from .construction import all_level_k_generators, dgenerator_to_sdgenerator
from .side import (
    BidirectedEdgeSide,
    IsolatedNodeSide,
    UndirEdgeSide,
)

__all__ = [
    "SemiDirectedGenerator",
    "UndirEdgeSide",
    "BidirectedEdgeSide",
    "IsolatedNodeSide",
    "dgenerator_to_sdgenerator",
    "all_level_k_generators",
    "attach_leaves_to_generator",
]

