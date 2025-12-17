"""
Classification functions for directed phylogenetic networks.

This module provides functions to classify and check properties of
directed phylogenetic networks (e.g., is_tree, is_binary, level, etc.).
"""

from typing import TYPE_CHECKING

from ...primitives.d_multigraph.features import has_parallel_edges as graph_has_parallel_edges

if TYPE_CHECKING:
    from . import DirectedPhyNetwork

def is_lsa_network(network: 'DirectedPhyNetwork') -> bool:
    """
    Check whether a directed phylogenetic network is an LSA network.
    
    An LSA (Least Stable Ancestor) network is one where the root node is the
    LSA node, i.e., the lowest node through which all root-to-leaf paths pass.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to check.
    
    Returns
    -------
    bool
        True if the network's root node equals its LSA node, False otherwise.
    
    Notes
    -----
    For empty networks this function returns True.
    
    Examples
    --------
    >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
    >>> is_lsa_network(net)
    True
    """
    if network.number_of_nodes() == 0:
        return True
    return network.root_node == network.LSA_node


def has_parallel_edges(network: 'DirectedPhyNetwork') -> bool:
    """
    Check if the network has any parallel edges.
    
    Parallel edges are multiple edges between the same pair of nodes in the same direction.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to check.
    
    Returns
    -------
    bool
        True if the network has at least one pair of parallel edges, False otherwise.
    
    Examples
    --------
    >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
    >>> has_parallel_edges(net)
    False
    """
    return graph_has_parallel_edges(network._graph)
