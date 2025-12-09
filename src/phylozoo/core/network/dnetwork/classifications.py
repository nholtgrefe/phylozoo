"""
Classification functions for directed phylogenetic networks.

This module provides functions to classify and check properties of
directed phylogenetic networks (e.g., is_tree, is_binary, level, etc.).
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from . import DirectedPhyNetwork

def is_LSA_network(network: 'DirectedPhyNetwork') -> bool:
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
    >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"})
    >>> is_LSA_network(net)
    True
    """
    if network.number_of_nodes() == 0:
        return True
    return network.root_node == network.LSA_node
