"""
Classification functions for directed phylogenetic networks.

This module provides functions to classify and check properties of
directed phylogenetic networks (e.g., is_tree, is_binary, level, etc.).

"""

from functools import lru_cache
from typing import TYPE_CHECKING, Any
import itertools

from ...primitives.d_multigraph.features import has_parallel_edges as graph_has_parallel_edges
from .features import blobs, omnians
from ....utils.exceptions import PhyloZooNotImplementedError, PhyloZooValueError, PhyloZooAlgorithmError


if TYPE_CHECKING:
    from . import DirectedPhyNetwork

@lru_cache(maxsize=128)
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


@lru_cache(maxsize=128)
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


def level(network: 'DirectedPhyNetwork') -> int:
    """
    Return the (strict) level of the network.
    
    The level is the maximum over all blobs of (number of hybrid edges minus
    number of hybrid nodes) in that blob.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network.
    
    Returns
    -------
    int
        The level of the network.
    
    Notes
    -----
    For empty networks, this function returns 0.
    
    Examples
    --------
    >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
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


def vertex_level(network: 'DirectedPhyNetwork') -> int:
    """
    Return the vertex level of the network.
    
    The vertex level is the maximum over all blobs of the number of hybrid nodes
    in that blob.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network.
    
    Returns
    -------
    int
        The vertex level of the network.
    
    Notes
    -----
    For empty networks, this function returns 0.
    
    Examples
    --------
    >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
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


def reticulation_number(network: 'DirectedPhyNetwork') -> int:
    """
    Return the reticulation number of the network.
    
    The reticulation number is the total number of hybrid edges minus the total
    number of hybrid nodes.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network.
    
    Returns
    -------
    int
        The reticulation number of the network.
    
    Examples
    --------
    >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
    >>> reticulation_number(net)
    0
    """
    return len(network.hybrid_edges) - len(network.hybrid_nodes)


@lru_cache(maxsize=128)
def is_binary(network: 'DirectedPhyNetwork') -> bool:
    """
    Check if the network is binary.
    
    A network is binary if every internal node has degree exactly 3, except
    for the root node which must have degree exactly 2.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to check.
    
    Returns
    -------
    bool
        True if the network is binary, False otherwise.
    
    Notes
    -----
    For empty networks or single-node networks, this function returns True.
    
    Examples
    --------
    >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
    >>> is_binary(net)
    True
    """
    if network.number_of_nodes() == 0:
        return True
    
    internal_nodes = network.internal_nodes
    root = network.root_node
    
    # Single-node network: root and leaf are the same node (degree 0)
    if network.number_of_nodes() == 1:
        return True
    
    # Check root node separately (root is not in internal_nodes)
    root_degree = network.degree(root)
    if root_degree != 2:
        return False
    
    # Check all other internal nodes
    for node in internal_nodes:
        node_degree = network.degree(node)
        # All internal nodes (non-root) must have degree 3
        if node_degree != 3:
            return False
    
    return True


@lru_cache(maxsize=128)
def is_tree(network: 'DirectedPhyNetwork') -> bool:
    """
    Check if the network is a tree.
    
    A network is a tree if it has no hybrid edges.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to check.
    
    Returns
    -------
    bool
        True if the network is a tree, False otherwise.
    
    Examples
    --------
    >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
    >>> is_tree(net)
    True
    """
    return len(network.hybrid_edges) == 0


@lru_cache(maxsize=128)
def is_simple(network: 'DirectedPhyNetwork') -> bool:
    """
    Check if the network is simple.
    
    A network is simple if it has at most one non-leaf blob.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to check.
    
    Returns
    -------
    bool
        True if the network has at most one non-leaf blob, False otherwise.
    
    Notes
    -----
    For empty networks, this function returns True.
    
    Examples
    --------
    >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
    >>> is_simple(net)
    True
    """
    if network.number_of_nodes() == 0:
        return True
    
    non_leaf_blobs = list(blobs(network, leaves=False))
    return len(non_leaf_blobs) <= 1


@lru_cache(maxsize=128)
def is_galled(network: 'DirectedPhyNetwork') -> bool:
    """
    Check if the network is galled.
    
    A network is galled if no hybrid node is ancestral to another hybrid node in the same blob.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to check.
    
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
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> is_galled(net)
    True
    
    >>> # Network with single hybrid in its own blob (galled)
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         (7, 5), (7, 6),  # Root to tree nodes
    ...         (5, 4), (6, 4),  # Both lead to hybrid 4
    ...         (4, 8),  # Hybrid to tree node
    ...         (8, 1), (8, 2)  # Tree node to leaves
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> is_galled(net)
    True
    
    >>> # Network with hybrid ancestral to another hybrid in same blob (not galled)
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         (9, 5), (9, 6),  # Root to tree nodes
    ...         (5, 4), (6, 4),  # Both lead to hybrid 4
    ...         (4, 7), (8, 7),  # Hybrid 4 and tree node 8 lead to hybrid 7
    ...         (7, 1)  # Hybrid 7 to leaf
    ...     ],
    ...     nodes=[(1, {'label': 'A'})]
    ... )
    >>> is_galled(net)
    False
    """
    import networkx as nx
    
    if network.number_of_nodes() == 0:
        return True
    
    hybrid_nodes = network.hybrid_nodes
    
    # If no hybrid nodes, network is galled (it's a tree)
    if not hybrid_nodes:
        return True
    
    # Get all blobs
    blob_list = blobs(network, trivial=False, leaves=False)
    
    # Get the underlying NetworkX graph for path checking
    nx_graph = network._graph._graph
    
    # Check each blob
    for blob_set in blob_list:
        # Get hybrid nodes in this blob
        hybrids_in_blob = [h for h in hybrid_nodes if h in blob_set]
        
        # If there's only one hybrid in the blob, it can't be ancestral to another
        if len(hybrids_in_blob) <= 1:
            continue
        
        # Check if any hybrid in this blob is ancestral to another hybrid in the same blob
        for h1 in hybrids_in_blob:
            for h2 in hybrids_in_blob:
                if h1 != h2:
                    # Check if h1 is ancestral to h2 (there's a directed path from h1 to h2)
                    if nx.has_path(nx_graph, h1, h2):
                        return False
    
    return True

@lru_cache(maxsize=128)
def is_stackfree(network: 'DirectedPhyNetwork') -> bool:
    """
    Check if the network is stack-free.
    
    A network is stack-free if no hybrid node has another hybrid node as its child.
    In other words, there are no "stacked" hybrid nodes.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to check.
    
    Returns
    -------
    bool
        True if the network is stack-free (no hybrid has a hybrid child), False otherwise.
    
    Examples
    --------
    >>> # Network with no hybrids (stack-free)
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> is_stackfree(net)
    True
    
    >>> # Network with hybrid that has tree node child (stack-free)
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         (7, 5), (7, 6),  # Root to tree nodes
    ...         (5, 4), (6, 4),  # Both lead to hybrid 4
    ...         (4, 8),  # Hybrid to tree node
    ...         (8, 1), (8, 2)  # Tree node to leaves
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> is_stackfree(net)
    True
    
    >>> # Network with stacked hybrids (not stack-free)
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         (9, 5), (9, 6),  # Root to tree nodes
    ...         (5, 4), (6, 4),  # Both lead to hybrid 4
    ...         (4, 7), (8, 7),  # Hybrid 4 and tree node 8 lead to hybrid 7
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
    
    # Check each hybrid node: if its child is also a hybrid, it's stacked
    for hybrid in hybrid_nodes:
        # Hybrid nodes have out-degree 1, so get the single child
        children = list(network.children(hybrid))
        if children:
            child = children[0]
            # If the child is also a hybrid node, we have a stack
            if child in hybrid_nodes:
                return False
    
    return True


@lru_cache(maxsize=128)
def is_treechild(network: 'DirectedPhyNetwork') -> bool:
    """
    Check if the network is tree-child.
    
    A phylogenetic network is tree-child if every internal vertex has at least
    one child that is not a hybrid node.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to check.
    
    Returns
    -------
    bool
        True if the network is tree-child, False otherwise.
    
    Notes
    -----
    For empty networks or single-node networks, this function returns True.
    All trees are tree-child networks.
    
    Examples
    --------
    >>> # Tree (tree-child)
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> is_treechild(net)
    True
    
    >>> # Network with hybrid that has tree node child (tree-child)
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         (7, 5), (7, 6),  # Root to tree nodes
    ...         (5, 4), (6, 4),  # Both lead to hybrid 4
    ...         (4, 8),  # Hybrid to tree node
    ...         (8, 1), (8, 2)  # Tree node to leaves
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> is_treechild(net)
    True
    
    >>> # Network where internal node has only hybrid children (not tree-child)
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         (9, 5), (9, 6),  # Root to tree nodes
    ...         (5, 4), (6, 4),  # Both lead to hybrid 4
    ...         (5, 7), (6, 7),  # Both lead to hybrid 7
    ...         (4, 1), (7, 2)  # Hybrids to leaves
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> is_treechild(net)
    False
    """
    if network.number_of_nodes() == 0:
        return True
    
    hybrid_nodes = network.hybrid_nodes
    leaves = network.leaves
    
    # Check all non-leaf vertices (root + internal nodes)
    # Each must have at least one child that is not a hybrid node
    for node in network._graph.nodes:
        # Skip leaves (they have no children)
        if node in leaves:
            continue
        
        children_list = list(network.children(node))
        
        # If node has no children, skip (shouldn't happen, but be safe)
        if not children_list:
            continue
        
        # Check if all children are hybrid nodes
        all_children_are_hybrids = all(child in hybrid_nodes for child in children_list)
        
        if all_children_are_hybrids:
            return False
    
    return True

@lru_cache(maxsize=128)
def is_treebased(network: 'DirectedPhyNetwork') -> bool:
    """
    Check if the network is tree-based.

    A network is tree-based if it is a base-tree with additional arcs.

    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to check.
    
    Returns
    -------
    bool
        True if the network is tree-based, False otherwise.
    
    Raises
    ------
    PhyloZooNotImplementedError
        If the network is non-binary or has parallel edges.
    
    Notes
    -----
    For empty networks or single-node networks, this function returns True.
    The implementation uses the omnian characterization of tree-based networks :cite:`Jetten2016`.

    Examples
    --------
    >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
    >>> is_treebased(net)
    True
    """

    if not is_binary(network):
        raise PhyloZooNotImplementedError("is_treebased function is not implemented for non-binary networks.")
    if has_parallel_edges(network):
        raise PhyloZooNotImplementedError("is_treebased function is not implemented for networks with parallel edges.")

    # Let U be the set of all omnians of N. Then N is tree-based
    # if and only if for every subset S ⊆ U the number of different
    # children of the vertices in S is greater than or equal to |S|.

    omnian_set = omnians(network)
    if not omnian_set:
        return True

    U = sorted(omnian_set)
    children = {u: set(network.children(u)) for u in U}

    for r in range(1, len(U) + 1):
        for subset in itertools.combinations(U, r):
            neigh: set[object] = set()
            for u in subset:
                neigh |= children[u]
            if len(neigh) < len(subset):
                return False

    return True


@lru_cache(maxsize=128)
def is_ultrametric(network: 'DirectedPhyNetwork') -> bool:
    """
    Check if the network is ultrametric.

    A network is ultrametric if all root-to-leaf distances are equal.
    All paths from root to each leaf must have the same distance.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to check.
    
    Returns
    -------
    bool
        True if the network is ultrametric, False otherwise.
    
    Raises
    ------
    PhyloZooValueError
        If any edge lacks a branch length.
    PhyloZooAlgorithmError
        If no path exists from root to a leaf.
    
    Notes
    -----
    For empty networks or single-node networks, this function returns True.
    Every edge must have a branch length specified.
    If parallel edges have different branch lengths, the network is not ultrametric (returns False).
    The distance from root to a leaf is the sum of branch lengths along a directed path.
    
    Examples
    --------
    >>> # Tree with equal distances (ultrametric)
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         {'u': 3, 'v': 1, 'branch_length': 1.0},
    ...         {'u': 3, 'v': 2, 'branch_length': 1.0}
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> is_ultrametric(net)
    True
    
    >>> # Tree with different distances (not ultrametric)
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         {'u': 3, 'v': 1, 'branch_length': 1.0},
    ...         {'u': 3, 'v': 2, 'branch_length': 2.0}
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> is_ultrametric(net)
    False
    """
    import networkx as nx
    
    if network.number_of_nodes() == 0:
        return True
    
    leaves = network.leaves
    if not leaves:
        return True
    
    root = network.root_node
    nx_graph = network._graph._graph
    
    # Step 1: Check that all edges have branch lengths
    for u, v, key, data in network._graph.edges(keys=True, data=True):
        bl = data.get('branch_length')
        if bl is None:
            raise PhyloZooValueError(
                f"Edge ({u}, {v}, {key}) lacks a branch_length. "
                "All edges must have branch_length for is_ultrametric."
            )
    
    # Step 2: Check that all parallel edges have the same branch length
    # If parallel edges have different branch lengths, the network is not ultrametric
    # Group edges by (u, v) pair
    edge_groups: dict[tuple[Any, Any], list[float]] = {}
    for u, v, key, data in network._graph.edges(keys=True, data=True):
        edge_key = (u, v)
        bl = data.get('branch_length')
        if edge_key not in edge_groups:
            edge_groups[edge_key] = []
        edge_groups[edge_key].append(bl)
    
    # Check parallel edges have same branch length
    for (u, v), branch_lengths in edge_groups.items():
        if len(branch_lengths) > 1:  # Parallel edges
            first_bl = branch_lengths[0]
            if not all(abs(bl - first_bl) < 1e-10 for bl in branch_lengths):
                return False  # Not ultrametric if parallel edges have different branch lengths
    
    # Step 3: Compute distances for all paths from root to each leaf
    distances: list[float] = []
    
    for leaf in leaves:
        # Find all paths from root to leaf
        paths = list(nx.all_simple_paths(nx_graph, root, leaf))
        if not paths:
            raise PhyloZooAlgorithmError(f"No path from root {root} to leaf {leaf}")
        
        for path in paths:
            # Sum branch lengths along the path
            total_distance = 0.0
            for i in range(len(path) - 1):
                u, v = path[i], path[i + 1]
                # Get branch length (if parallel edges exist, they all have the same length)
                # Use the first edge's branch length
                edges_data = nx_graph[u].get(v, {})
                if not edges_data:
                    raise PhyloZooAlgorithmError(f"Edge ({u}, {v}) not found in graph")
                first_key = next(iter(edges_data))
                bl = edges_data[first_key].get('branch_length')
                if bl is None:
                    raise PhyloZooValueError(f"Edge ({u}, {v}) lacks a branch_length")
                total_distance += bl
            distances.append(total_distance)
    
    # Check if all distances are equal
    if not distances:
        return True
    
    first_distance = distances[0]
    return all(abs(d - first_distance) < 1e-10 for d in distances)

@lru_cache(maxsize=128)
def is_normal(network: 'DirectedPhyNetwork') -> bool:
    """
    Check if the network is normal.
    
    A reticulation arc (u, v) is a shortcut if there is a directed path from u to v
    that does not traverse (u, v). A normal network is a tree-child network without shortcuts.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to check.
    
    Returns
    -------
    bool
        True if the network is normal, False otherwise.
    
    Notes
    -----
    For empty networks or single-node networks, this function returns True.
    All trees are normal networks.
    Networks with parallel edges are never normal (parallel edges are shortcuts).
    
    Examples
    --------
    >>> # Tree (normal)
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> is_normal(net)
    True
    
    >>> # Tree-child network without shortcuts (normal)
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         (7, 5), (7, 6),  # Root to tree nodes
    ...         (5, 4), (6, 4),  # Both lead to hybrid 4
    ...         (4, 8),  # Hybrid to tree node
    ...         (8, 1), (8, 2)  # Tree node to leaves
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> is_normal(net)
    True
    
    >>> # Network with shortcut (not normal)
    >>> # Edge (5, 4) is a shortcut because path 5 -> 9 -> 4 exists
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         (10, 5), (10, 6),  # Root to tree nodes
    ...         (5, 4), (5, 9), (6, 4),  # Tree nodes lead to hybrid 4
    ...         (9, 4),  # This creates a shortcut: path 5 -> 9 -> 4 bypasses edge (5, 4)
    ...         (4, 1)  # Hybrid to leaf
    ...     ],
    ...     nodes=[(1, {'label': 'A'})]
    ... )
    >>> is_normal(net)
    False
    """
    import networkx as nx
    
    if network.number_of_nodes() == 0:
        return True
    
    # Step 1: Check if tree-child (normal networks must be tree-child)
    if not is_treechild(network):
        return False
    
    # Step 2: Check for parallel edges (they are always shortcuts)
    if has_parallel_edges(network):
        return False
    
    # Step 3: Check each hybrid edge for shortcuts
    # A hybrid edge (u, v) is a shortcut if there's a directed path from u to v
    # that doesn't use the edge (u, v)
    hybrid_edges = network.hybrid_edges
    
    # If no hybrid edges, network is a tree (normal)
    if not hybrid_edges:
        return True
    
    # Get the underlying NetworkX graph
    nx_graph = network._graph._graph
    
    # For each hybrid edge, check if there's an alternative path
    for u, v, key in hybrid_edges:
        # Create a copy of the graph without this specific edge
        graph_without_edge = nx_graph.copy()
        graph_without_edge.remove_edge(u, v, key)
        
        # Check if there's still a path from u to v
        if nx.has_path(graph_without_edge, u, v):
            return False
    
    return True
