"""
Classification functions for semi-directed and mixed phylogenetic networks.

This module provides functions to classify and check properties of
semi-directed and mixed phylogenetic networks (e.g., is_tree, is_binary, level, etc.).
"""

from functools import lru_cache
from typing import TYPE_CHECKING

from ...primitives.m_multigraph.features import has_parallel_edges as graph_has_parallel_edges
from .features import blobs
from ....utils.exceptions import PhyloZooNotImplementedError

if TYPE_CHECKING:
    from .base import SemiDirectedPhyNetwork
    from .base import MixedPhyNetwork


@lru_cache(maxsize=128)
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


def level(network: 'SemiDirectedPhyNetwork') -> int:
    """
    Return the level of the network.
    
    The level is the maximum over all blobs of (number of hybrid edges minus
    number of hybrid nodes) in that blob.
    
    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network.
    
    Returns
    -------
    int
        The level of the network.
    
    Notes
    -----
    For empty networks, this function returns 0.
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> net = SemiDirectedPhyNetwork(undirected_edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
    >>> level(net)
    0
    """
    if network.number_of_nodes() == 0:
        return 0
    
    hybrid_edges = network.hybrid_edges
    hybrid_nodes = network.hybrid_nodes
    
    max_level = 0
    for blob_set in blobs(network, trivial=False, leaves=False):
        # Count hybrid edges in this blob (both endpoints must be in blob)
        hybrid_edges_in_blob = sum(
            1 for u, v, _ in hybrid_edges
            if u in blob_set and v in blob_set
        )
        # Count hybrid nodes in this blob
        hybrid_nodes_in_blob = sum(1 for node in hybrid_nodes if node in blob_set)
        
        blob_level = hybrid_edges_in_blob - hybrid_nodes_in_blob
        max_level = max(max_level, blob_level)
    
    return max_level


def vertex_level(network: 'SemiDirectedPhyNetwork') -> int:
    """
    Return the vertex level of the network.
    
    The vertex level is the maximum over all blobs of the number of hybrid nodes
    in that blob.
    
    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network.
    
    Returns
    -------
    int
        The vertex level of the network.
    
    Notes
    -----
    For empty networks, this function returns 0.
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> net = SemiDirectedPhyNetwork(undirected_edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
    >>> vertex_level(net)
    0
    """
    if network.number_of_nodes() == 0:
        return 0
    
    hybrid_nodes = network.hybrid_nodes
    
    max_vertex_level = 0
    for blob_set in blobs(network, trivial=False, leaves=False):
        hybrid_nodes_in_blob = sum(1 for node in hybrid_nodes if node in blob_set)
        max_vertex_level = max(max_vertex_level, hybrid_nodes_in_blob)
    
    return max_vertex_level


def reticulation_number(network: 'SemiDirectedPhyNetwork') -> int:
    """
    Return the reticulation number of the network.
    
    The reticulation number is the total number of hybrid edges minus the total
    number of hybrid nodes.
    
    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network.
    
    Returns
    -------
    int
        The reticulation number of the network.
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> net = SemiDirectedPhyNetwork(undirected_edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
    >>> reticulation_number(net)
    0
    """
    return len(network.hybrid_edges) - len(network.hybrid_nodes)


@lru_cache(maxsize=128)
def is_binary(network: 'MixedPhyNetwork') -> bool:
    """
    Check if the network is binary.
    
    A network is binary if every internal node has degree exactly 3.
    
    Parameters
    ----------
    network : MixedPhyNetwork
        The mixed phylogenetic network to check.
    
    Returns
    -------
    bool
        True if the network is binary, False otherwise.
    
    Notes
    -----
    For empty networks or single-node networks, this function returns True.
    Unlike directed networks, semi-directed networks have no root node, so all
    internal nodes must have degree 3.
    
    Examples
    --------
    >>> net = MixedPhyNetwork(undirected_edges=[(3, 1), (3, 2), (3, 4)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})])
    >>> is_binary(net)
    True
    """
    if network.number_of_nodes() == 0:
        return True
    
    internal_nodes = network.internal_nodes
    
    for node in internal_nodes:
        node_degree = network.degree(node)
        # All internal nodes must have degree 3
        if node_degree != 3:
            return False
    
    return True


@lru_cache(maxsize=128)
def is_tree(network: 'MixedPhyNetwork') -> bool:
    """
    Check if the network is a tree.
    
    A network is a tree if it has no hybrid edges.
    
    Parameters
    ----------
    network : MixedPhyNetwork
        The mixed phylogenetic network to check.
    
    Returns
    -------
    bool
        True if the network is a tree, False otherwise.
    
    Examples
    --------
    >>> net = MixedPhyNetwork(undirected_edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
    >>> is_tree(net)
    True
    """
    return len(network.hybrid_edges) == 0


@lru_cache(maxsize=128)
def is_simple(network: 'SemiDirectedPhyNetwork') -> bool:
    """
    Check if the network is simple.
    
    A network is simple if it has at most one non-leaf blob.
    
    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network to check.
    
    Returns
    -------
    bool
        True if the network has at most one non-leaf blob, False otherwise.
    
    Notes
    -----
    For empty networks, this function returns True.
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> net = SemiDirectedPhyNetwork(undirected_edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
    >>> is_simple(net)
    True
    """
    if network.number_of_nodes() == 0:
        return True
    
    non_leaf_blobs = list(blobs(network, leaves=False))
    return len(non_leaf_blobs) <= 1


def is_galled(network: 'SemiDirectedPhyNetwork') -> bool:
    """Stub for is_galled function."""
    raise PhyloZooNotImplementedError("is_galled function is not implemented.")

def is_OLP(network: 'SemiDirectedPhyNetwork') -> bool:
    """Stub for is_OLP function."""
    raise PhyloZooNotImplementedError("is_OLP function is not implemented.")

@lru_cache(maxsize=128)
def is_stackfree(network: 'SemiDirectedPhyNetwork') -> bool:
    """
    Check if the network is stack-free.
    
    A network is stack-free if no hybrid node has another hybrid node as its child.
    In a semi-directed network, hybrid nodes have exactly one outgoing directed edge,
    so we check if that child is also a hybrid node.
    
    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network to check.
    
    Returns
    -------
    bool
        True if the network is stack-free (no hybrid has a hybrid child), False otherwise.
    
    Examples
    --------
    >>> # Network with no hybrids (stack-free)
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 4)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
    ... )
    >>> is_stackfree(net)
    True
    
    >>> # Network with hybrid that has tree node child (stack-free)
    >>> net = SemiDirectedPhyNetwork(
    ...     directed_edges=[
    ...         (5, 4), (6, 4)  # Both lead to hybrid 4
    ...     ],
    ...     undirected_edges=[
    ...         (7, 5), (7, 6),  # Root to tree nodes
    ...         (4, 8),  # Hybrid to tree node
    ...         (8, 1), (8, 2)  # Tree node to leaves
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> is_stackfree(net)
    True
    
    >>> # Network with stacked hybrids (not stack-free)
    >>> net = SemiDirectedPhyNetwork(
    ...     directed_edges=[
    ...         (5, 4), (6, 4),  # Both lead to hybrid 4
    ...         (4, 7), (8, 7)  # Hybrid 4 and tree node 8 lead to hybrid 7
    ...     ],
    ...     undirected_edges=[
    ...         (9, 5), (9, 6),  # Root to tree nodes
    ...         (7, 1)  # Hybrid 7 to leaf
    ...     ],
    ...     nodes=[(1, {'label': 'A'})]
    ... )
    >>> is_stackfree(net)
    False
    """
    hybrid_nodes = network.hybrid_nodes
    
    # If no hybrid nodes, network is stack-free
    if not hybrid_nodes:
        return True
    
    # Check each hybrid node: if its child (via directed edge) is also a hybrid, it's stacked
    for hybrid in hybrid_nodes:
        # Hybrid nodes have exactly one outgoing directed edge
        # Get the child via the directed edge using incident_child_edges
        out_edges = list(network.incident_child_edges(hybrid, keys=False, data=False))
        if out_edges:
            # Get the target of the first (and only) outgoing directed edge
            child = out_edges[0][1]  # (u, v) tuple, get v
            # If the child is also a hybrid node, we have a stack
            if child in hybrid_nodes:
                return False
    
    return True


def is_strongly_treechild(network: 'SemiDirectedPhyNetwork') -> bool:
    """Stub for is_strongly_treechild function."""
    raise PhyloZooNotImplementedError("is_strongly_treechild function is not implemented.")

def is_weakly_treechild(network: 'SemiDirectedPhyNetwork') -> bool:
    """Stub for is_weakly_treechild function."""
    raise PhyloZooNotImplementedError("is_weakly_treechild function is not implemented.")

def is_strongly_treebased(network: 'SemiDirectedPhyNetwork') -> bool:
    """Stub for is_strongly_treebased function."""
    raise PhyloZooNotImplementedError("is_strongly_treebased function is not implemented.")

def is_weakly_treebased(network: 'SemiDirectedPhyNetwork') -> bool:
    """Stub for is_weakly_treebased function."""
    raise PhyloZooNotImplementedError("is_weakly_treebased function is not implemented.")