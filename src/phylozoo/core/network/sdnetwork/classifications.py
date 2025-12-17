"""
Classification functions for semi-directed and mixed phylogenetic networks.

This module provides functions to classify and check properties of
semi-directed and mixed phylogenetic networks (e.g., is_tree, is_binary, level, etc.).
"""

from typing import TYPE_CHECKING

from ...primitives.m_multigraph.features import has_parallel_edges as graph_has_parallel_edges

if TYPE_CHECKING:
    from .base import SemiDirectedPhyNetwork
    from .base import MixedPhyNetwork


def has_parallel_edges(network: 'MixedPhyNetwork') -> bool:
    """
    Check if the network has any parallel edges.
    
    Parallel edges are multiple edges between the same pair of nodes,
    either as multiple directed edges in the same direction or multiple undirected edges.
    
    Parameters
    ----------
    network : MixedPhyNetwork
        The mixed phylogenetic network to check.
    
    Returns
    -------
    bool
        True if the network has at least one pair of parallel edges, False otherwise.
    
    Examples
    --------
    >>> net = MixedPhyNetwork(undirected_edges=[(3, 1), (3, 2), (3, 4)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})])
    >>> has_parallel_edges(net)
    False
    """
    return graph_has_parallel_edges(network._graph)


# TODO: Implement additional classification functions here.
# These could include functions like:
# - is_tree(network) -> bool
# - is_binary(network) -> bool
# - level(network) -> int
# - is_time_consistent(network) -> bool
# - etc.
#
# Follow the NetworkX-style function-based API, where functions take a network
# instance as the first argument.

