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
    Some authors call a simple network a "bloblet (network)".
    
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


@lru_cache(maxsize=128)
def is_galled(network: 'SemiDirectedPhyNetwork') -> bool:
    """
    Check if the network is galled.
    
    A network is galled if no hybrid node is ancestral to another hybrid node in the same blob.
    A hybrid node is ancestral to another one if there exists an up-down path from one to the other.
    
    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network to check.
    
    Returns
    -------
    bool
        True if the network is galled, False otherwise.
    
    Notes
    -----
    For empty networks or networks with no hybrid nodes, this function returns True.
    All trees are galled networks.
    
    Examples
    --------
    >>> # Network with no hybrid nodes (galled)
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> is_galled(net)
    True
    
    >>> # Network with single hybrid in its own blob (galled)
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
    >>> is_galled(net)
    True
    
    >>> # Network with hybrid ancestral to another hybrid in same blob (not galled)
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
    >>> is_galled(net)
    False
    """
    from ...primitives.m_multigraph.features import updown_path_vertices
    
    if network.number_of_nodes() == 0:
        return True
    
    hybrid_nodes = network.hybrid_nodes
    
    # If no hybrid nodes, network is galled (it's a tree)
    if not hybrid_nodes:
        return True
    
    # Get all blobs
    blob_list = blobs(network, trivial=False, leaves=False)
    
    # Check each blob
    for blob_set in blob_list:
        # Get hybrid nodes in this blob
        hybrids_in_blob = [h for h in hybrid_nodes if h in blob_set]
        
        # If there's only one hybrid in the blob, it can't be ancestral to another
        if len(hybrids_in_blob) <= 1:
            continue
        
        # Check if any hybrid in this blob is ancestral to another hybrid in the same blob
        # (i.e., there exists an up-down path from one to the other)
        for h1 in hybrids_in_blob:
            for h2 in hybrids_in_blob:
                if h1 != h2:
                    # Check if there's an up-down path from h1 to h2
                    # If h2 is in the vertices on up-down paths from h1, then h1 is ancestral to h2
                    path_vertices = updown_path_vertices(network._graph, h1, h2)
                    if h2 in path_vertices:
                        return False
    
    return True

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


def _edge_root_locations(network: 'SemiDirectedPhyNetwork') -> list[tuple]:
    """
    Get all valid edge root locations (undirected + directed edges).

    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network.

    Returns
    -------
    list
        List of edge root locations as (u, v, key) tuples.
    """
    from .features import root_locations

    node_locs, undir_locs, dir_locs = root_locations(network)
    return list(undir_locs) + list(dir_locs)


@lru_cache(maxsize=128)
def is_strongly_treechild(network: 'SemiDirectedPhyNetwork') -> bool:
    """
    Check if the network is strongly tree-child.

    A semi-directed network is strongly tree-child if every rooting (on an edge)
    yields a directed network that is tree-child.

    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network.

    Returns
    -------
    bool
        True if every edge rooting is tree-child, False otherwise.

    Raises
    ------
    PhyloZooNotImplementedError
        If the network is non-binary or has parallel edges.

    Notes
    -----
    Only rootings on edges are considered; node rootings are ignored.
    For empty or single-node networks, returns True (vacuously).

    Examples
    --------
    >>> from phylozoo.core.network import SemiDirectedPhyNetwork
    >>> from phylozoo.core.network.sdnetwork.classifications import (
    ...     is_strongly_treechild,
    ... )
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 4)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})],
    ... )
    >>> is_strongly_treechild(net)
    True
    """
    from .derivations import to_d_network
    from ...network.dnetwork.classifications import is_treechild

    if network.number_of_nodes() <= 1:
        return True
    if not is_binary(network):
        raise PhyloZooNotImplementedError(
            "is_strongly_treechild is not implemented for non-binary networks."
        )
    if has_parallel_edges(network):
        raise PhyloZooNotImplementedError(
            "is_strongly_treechild is not implemented for networks with parallel edges."
        )

    edge_locs = _edge_root_locations(network)
    if not edge_locs:
        return True  # vacuously: no edge rootings to check

    for loc in edge_locs:
        d_net = to_d_network(network, root_location=loc)
        if not is_treechild(d_net):
            return False
    return True


@lru_cache(maxsize=128)
def is_weakly_treechild(network: 'SemiDirectedPhyNetwork') -> bool:
    """
    Check if the network is weakly tree-child.

    A semi-directed network is weakly tree-child if there exists at least one
    rooting (on an edge) that yields a directed network that is tree-child.

    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network.

    Returns
    -------
    bool
        True if some edge rooting is tree-child, False otherwise.

    Raises
    ------
    PhyloZooNotImplementedError
        If the network is non-binary or has parallel edges.

    Notes
    -----
    Only rootings on edges are considered; node rootings are ignored.
    For empty or single-node networks, returns True (vacuously).

    Examples
    --------
    >>> from phylozoo.core.network import SemiDirectedPhyNetwork
    >>> from phylozoo.core.network.sdnetwork.classifications import (
    ...     is_weakly_treechild,
    ... )
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 4)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})],
    ... )
    >>> is_weakly_treechild(net)
    True
    """
    from .derivations import to_d_network
    from ...network.dnetwork.classifications import is_treechild

    if network.number_of_nodes() <= 1:
        return True
    if not is_binary(network):
        raise PhyloZooNotImplementedError(
            "is_weakly_treechild is not implemented for non-binary networks."
        )
    if has_parallel_edges(network):
        raise PhyloZooNotImplementedError(
            "is_weakly_treechild is not implemented for networks with parallel edges."
        )

    edge_locs = _edge_root_locations(network)
    if not edge_locs:
        return False

    for loc in edge_locs:
        d_net = to_d_network(network, root_location=loc)
        if is_treechild(d_net):
            return True
    return False


@lru_cache(maxsize=128)
def is_strongly_treebased(network: 'SemiDirectedPhyNetwork') -> bool:
    """
    Check if the network is strongly tree-based.

    A semi-directed network is strongly tree-based if every rooting (on an edge)
    yields a directed network that is tree-based.

    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network.

    Returns
    -------
    bool
        True if every edge rooting is tree-based, False otherwise.

    Raises
    ------
    PhyloZooNotImplementedError
        If the network is non-binary or has parallel edges.

    Notes
    -----
    Only rootings on edges are considered; node rootings are ignored.
    For empty or single-node networks, returns True (vacuously).

    Examples
    --------
    >>> from phylozoo.core.network import SemiDirectedPhyNetwork
    >>> from phylozoo.core.network.sdnetwork.classifications import (
    ...     is_strongly_treebased,
    ... )
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 4)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})],
    ... )
    >>> is_strongly_treebased(net)
    True
    """
    from .derivations import to_d_network
    from ...network.dnetwork.classifications import is_treebased

    if network.number_of_nodes() <= 1:
        return True
    if not is_binary(network):
        raise PhyloZooNotImplementedError(
            "is_strongly_treebased is not implemented for non-binary networks."
        )
    if has_parallel_edges(network):
        raise PhyloZooNotImplementedError(
            "is_strongly_treebased is not implemented for networks with parallel edges."
        )

    edge_locs = _edge_root_locations(network)
    if not edge_locs:
        return True  # vacuously

    for loc in edge_locs:
        d_net = to_d_network(network, root_location=loc)
        if not is_treebased(d_net):
            return False
    return True


@lru_cache(maxsize=128)
def is_weakly_treebased(network: 'SemiDirectedPhyNetwork') -> bool:
    """
    Check if the network is weakly tree-based.

    A semi-directed network is weakly tree-based if there exists at least one
    rooting (on an edge) that yields a directed network that is tree-based.

    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network.

    Returns
    -------
    bool
        True if some edge rooting is tree-based, False otherwise.

    Raises
    ------
    PhyloZooNotImplementedError
        If the network is non-binary or has parallel edges.

    Notes
    -----
    Only rootings on edges are considered; node rootings are ignored.
    For empty or single-node networks, returns True (vacuously).

    Examples
    --------
    >>> from phylozoo.core.network import SemiDirectedPhyNetwork
    >>> from phylozoo.core.network.sdnetwork.classifications import (
    ...     is_weakly_treebased,
    ... )
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 4)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})],
    ... )
    >>> is_weakly_treebased(net)
    True
    """
    from .derivations import to_d_network
    from ...network.dnetwork.classifications import is_treebased

    if network.number_of_nodes() <= 1:
        return True
    if not is_binary(network):
        raise PhyloZooNotImplementedError(
            "is_weakly_treebased is not implemented for non-binary networks."
        )
    if has_parallel_edges(network):
        raise PhyloZooNotImplementedError(
            "is_weakly_treebased is not implemented for networks with parallel edges."
        )

    edge_locs = _edge_root_locations(network)
    if not edge_locs:
        return False

    for loc in edge_locs:
        d_net = to_d_network(network, root_location=loc)
        if is_treebased(d_net):
            return True
    return False