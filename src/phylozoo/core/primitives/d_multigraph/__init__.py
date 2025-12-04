"""
Directed multi-graph module.

This module provides the DirectedMultiGraph class and related operations.
"""

from .dm_graph import DirectedMultiGraph

# Import operations for convenience
from .dm_operations import (
    number_of_connected_components,
    is_connected,
    connected_components,
    identify_two_nodes,
    identify_node_set,
)
# Import conversion functions for convenience
from .dm_conversion import (
    digraph_to_directedmultigraph,
    multidigraph_to_directedmultigraph,
)

__all__ = [
    "DirectedMultiGraph",
    "number_of_connected_components",
    "is_connected",
    "connected_components",
    "identify_two_nodes",
    "identify_node_set",
    "digraph_to_directedmultigraph",
    "multidigraph_to_directedmultigraph",
]

