"""
Mixed multi-graph module.

This module provides the MixedMultiGraph class and related operations.
"""

from .mm_graph import MixedMultiGraph

# Import operations for convenience
from .mm_operations import (
    number_of_connected_components,
    is_connected,
    connected_components,
    identify_two_nodes,
    identify_node_set,
    source_components,
)
# Import conversion functions for convenience
from .mm_conversions import (
    graph_to_mixedmultigraph,
    multigraph_to_mixedmultigraph,
    multidigraph_to_mixedmultigraph,
    directedmultigraph_to_mixedmultigraph,
)

__all__ = [
    "MixedMultiGraph",
    "number_of_connected_components",
    "is_connected",
    "connected_components",
    "identify_two_nodes",
    "identify_node_set",
    "source_components",
    "graph_to_mixedmultigraph",
    "multigraph_to_mixedmultigraph",
    "multidigraph_to_mixedmultigraph",
    "directedmultigraph_to_mixedmultigraph",
]

