"""
Network features module.

This module provides functions to extract and identify features of directed
phylogenetic networks (e.g., LSA node, blobs, omnians, etc.).
"""

import warnings
from collections import deque
from functools import lru_cache
from typing import Any, TypeVar

import networkx as nx

from ...primitives.d_multigraph.features import (
    bi_edge_connected_components,
    cut_edges as graph_cut_edges,
    cut_vertices as graph_cut_vertices,
)
from .base import DirectedPhyNetwork
from ....utils.exceptions import PhyloZooValueError, PhyloZooAlgorithmError, PhyloZooWarning

T = TypeVar('T')


def lsa_node(network: DirectedPhyNetwork) -> T:
    """
    Find the Least Stable Ancestor (LSA) node of a directed phylogenetic network.
    
    The LSA is the lowest node through which all paths from the root to the leaves pass.
    In other words, it is the unique node that is an ancestor of all leaves and is
    the lowest such node (has maximum depth from the root).
    
    Parameters
    ----------
    network : DirectedPhyNetwork[T]
        The directed phylogenetic network.
    
    Returns
    -------
    T
        The LSA node identifier.
    
    Raises
    ------
    PhyloZooValueError
        If the network is empty or has no leaves.
    PhyloZooAlgorithmError
        If there is no path from root to leaf.
    
    Examples
    --------
    >>> # LSA below the root (hybrid node 4 is on all root-to-leaf paths)
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         (7, 5), (7, 6),                 # root to tree nodes
    ...         (5, 4, 0), (5, 4, 1),           # parallel edges keep tree node 5 out-degree >= 2
    ...         (6, 4, 0), (6, 4, 1),           # parallel edges keep tree node 6 out-degree >= 2
    ...         (4, 10),                        # hybrid 4 (in-degree 4, out-degree 1) to tree node 10
    ...         (10, 1), (10, 2)                # tree node 10 splits to leaves
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> lsa_node(net)
    4
    >>> # In a simple tree, the LSA is just the root
    >>> net2 = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
    >>> lsa_node(net2)
    3
    """
    if network.number_of_nodes() == 0:
        raise PhyloZooValueError("Cannot find LSA node in empty network")
    
    leaves = network.leaves
    
    # If there's only one leaf, the LSA is that leaf's parent (or the leaf itself if it's the root)
    if len(leaves) == 1:
        leaf = next(iter(leaves))
        # If the leaf is the root (shouldn't happen in valid networks, but handle it)
        if network.indegree(leaf) == 0:
            return leaf
        # Otherwise, return the parent (there should be exactly one parent for a leaf)
        parents = list(network.parents(leaf))
        return parents[0] if parents else leaf
    
    # Find the LSA: the lowest node through which ALL paths from root to leaves pass
    # This means for each leaf, the LSA must be on ALL simple paths from root to that leaf
    root = network.root_node
    dag = network._graph._graph
    
    # For each leaf, find all nodes that appear on ALL simple paths from root to that leaf
    # The LSA must be in the intersection of all these sets
    nodes_on_all_paths_to_each_leaf = []
    
    for leaf in leaves:
        # Find all simple paths from root to this leaf
        try:
            all_paths = list(nx.all_simple_paths(dag, root, leaf))
        except nx.NetworkXNoPath:
            raise PhyloZooAlgorithmError(f"No path from root {root} to leaf {leaf}")
        
        if not all_paths:
            raise PhyloZooAlgorithmError(f"No path from root {root} to leaf {leaf}")
        
        # Find nodes that appear on ALL paths to this leaf
        # Start with nodes from the first path
        nodes_on_all_paths = set(all_paths[0])
        # Intersect with nodes from each subsequent path
        for path in all_paths[1:]:
            nodes_on_all_paths &= set(path)
        
        nodes_on_all_paths_to_each_leaf.append(nodes_on_all_paths)
    
    # Find intersection: nodes that are on ALL paths to ALL leaves
    lsa_candidates = nodes_on_all_paths_to_each_leaf[0]
    for node_set in nodes_on_all_paths_to_each_leaf[1:]:
        lsa_candidates &= node_set
    
    if not lsa_candidates:
        raise PhyloZooAlgorithmError("No node found on all paths from root to all leaves")
    
    # Find the deepest node (maximum depth from root) among candidates
    # Compute depths using BFS from root (use deque for O(1) popleft)
    depths: dict[T, int] = {}
    queue = deque([root])
    depths[root] = 0
    
    while queue:
        current = queue.popleft()
        for child in network.children(current):
            if child not in depths:
                depths[child] = depths[current] + 1
                queue.append(child)
    
    # Return the deepest node among LSA candidates
    return max(lsa_candidates, key=lambda node: depths.get(node, 0))


@lru_cache(maxsize=128)
def blobs(
    network: DirectedPhyNetwork,
    trivial: bool = True,
    leaves: bool = True,
) -> list[set[T]]:
    """
    Get blobs of the network.
    
    A blob is a maximal subgraph without any cut-edges. This function provides
    filtering options to control which blobs are returned.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network.
    trivial : bool, optional
        Whether to include trivial (single-node) blobs. By default True.
    leaves : bool, optional
        Whether to include blobs that contain only leaves. By default True.
    
    Returns
    -------
    list[set[T]]
        List of sets of nodes forming each blob.
    
    Raises
    ------
    PhyloZooValueError
        If `trivial=False` and `leaves=True` (this combination is not possible
        since leaves are single-node components).
    
    Notes
    -----
    Blobs are computed as bi-edge connected components (2-edge-connected components).
    A bi-edge connected component is a maximal subgraph that remains connected
    after removing any single edge (i.e., has no cut-edges/bridges).
    
    Results are cached using LRU cache with maxsize=128.
    
    Examples
    --------
    >>> # Network with hybrid node creating a non-trivial blob (cycle)
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         (8, 5), (8, 6),  # Root to tree nodes
    ...         (5, 1), (5, 2),  # Tree node 5: in=1, out=2
    ...         (6, 3), (6, 9),  # Tree node 6: in=1, out=2
    ...         (5, 4), (6, 4),  # Both lead to hybrid node 4
    ...         (4, 7),  # Hybrid 4: in=2, out=1
    ...         (7, 10), (7, 11),  # Tree node 7: in=1, out=2
    ...     ],
    ...     nodes=[
    ...         (1, {'label': 'A'}),
    ...         (2, {'label': 'B'}),
    ...         (3, {'label': 'C'}),
    ...         (9, {'label': 'D'}),
    ...         (10, {'label': 'E'}),
    ...         (11, {'label': 'F'}),
    ...     ]
    ... )
    >>> sorted([sorted(b) for b in blobs(net)])
    [[1], [2], [3], [4, 5, 6, 8], [7], [9], [10], [11]]
    >>> # Filtering: exclude trivial (single-node) blobs
    >>> len([b for b in blobs(net, trivial=False, leaves=False)])
    1
    >>> # Filtering: exclude blobs containing only leaves
    >>> len([b for b in blobs(net, leaves=False)])
    2
    """
    # Check for invalid parameter combination
    if not trivial and leaves:
        raise PhyloZooValueError(
            "Cannot have trivial=False and leaves=True: leaves are single-node "
            "components, so excluding trivial components would exclude all leaves."
        )
    
    leaves_set = network.leaves
    result: list[set[T]] = []
    
    # Process bi-edge connected components directly
    for blob in bi_edge_connected_components(network._graph):
        blob_set = set(blob)
        
        # Filter single-node components based on parameters
        if len(blob_set) == 1:
            node = next(iter(blob_set))
            # Skip if trivial=False and it's not a leaf
            if not trivial and node not in leaves_set:
                continue
            # Skip if leaves=False and it is a leaf
            if not leaves and node in leaves_set:
                continue
            result.append(blob_set)
        else:
            # Multi-node components: always include (cannot consist entirely of leaves)
            result.append(blob_set)
    
    return result


@lru_cache(maxsize=128)
def k_blobs(
    network: DirectedPhyNetwork,
    k: int,
    trivial: bool = True,
    leaves: bool = True,
) -> list[set[T]]:
    """
    Get k-blobs of the network.
    
    A k-blob is a blob with exactly k edges incident to it. An incident edge
    is any edge that connects a node inside the blob to a node outside the blob.
    Parallel edges are counted separately.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network.
    k : int
        The number of edges that should be incident to each returned blob.
    trivial : bool, optional
        Whether to include trivial (single-node) blobs. By default True.
    leaves : bool, optional
        Whether to include blobs that contain only leaves. By default True.
    
    Returns
    -------
    list[set[T]]
        List of sets of nodes forming each k-blob.
    
    Raises
    ------
    PhyloZooValueError
        If `trivial=False` and `leaves=True` (this combination is not possible
        since leaves are single-node components).
    
    Notes
    -----
    This function identifies blobs and then filters
    them based on the number of incident edges. Parallel edges are counted
    separately, so if there are two parallel edges crossing the blob boundary,
    they count as two incident edges.
    
    Results are cached using LRU cache with maxsize=128.
    
    Examples
    --------
    >>> # Tree network: leaves are 1-blobs, internal nodes are 2-blobs or more
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> sorted([sorted(b) for b in k_blobs(net, k=1)])
    [[1], [2]]
    >>> sorted([sorted(b) for b in k_blobs(net, k=2)])
    [[3]]
    """
    # Check for invalid parameter combination
    if not trivial and leaves:
        raise PhyloZooValueError(
            "Cannot have trivial=False and leaves=True: leaves are single-node "
            "components, so excluding trivial components would exclude all leaves."
        )
    
    # Get the underlying graph for edge access
    graph = network._graph
    result: list[set[T]] = []
    
    # Iterate through blobs
    for blob in blobs(network, trivial=trivial, leaves=leaves):
        blob_set = set(blob)
        
        # Count edges incident to this blob
        # An incident edge has exactly one endpoint in the blob
        incident_edge_count = 0
        
        # Check all edges incident to nodes in the blob
        for node in blob_set:
            # Check incoming edges (parent edges)
            for u, v, key in graph.incident_parent_edges(node, keys=True):
                if u not in blob_set:  # Edge crosses blob boundary
                    incident_edge_count += 1
            
            # Check outgoing edges (child edges)
            for u, v, key in graph.incident_child_edges(node, keys=True):
                if v not in blob_set:  # Edge crosses blob boundary
                    incident_edge_count += 1
        
        # Add blob if it has exactly k incident edges
        if incident_edge_count == k:
            result.append(blob_set)
    
    return result


@lru_cache(maxsize=128)
def cut_edges(network: DirectedPhyNetwork) -> set[tuple[T, T, int]]:
    """
    Find all cut-edges (bridges) in the network.
    
    A cut-edge is an edge whose removal increases the number of
    weakly connected components. Results are cached per network instance.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network.
    
    Returns
    -------
    set[tuple[T, T, int]]
        Set of cut-edges as 3-tuples (u, v, key).
    
    Examples
    --------
    >>> net = DirectedPhyNetwork(edges=[(1, 2), (2, 3)], nodes=[(3, {'label': 'A'})])
    >>> edges = cut_edges(net)
    >>> (1, 2, 0) in edges and (2, 3, 0) in edges
    True
    
    Notes
    -----
    Results are cached using LRU cache with maxsize=128.
    """
    return graph_cut_edges(network._graph, keys=True, data=False)


@lru_cache(maxsize=128)
def cut_vertices(network: DirectedPhyNetwork) -> set[T]:
    """
    Find all cut-vertices (articulation points) in the network.
    
    A cut-vertex is a vertex whose removal increases the number of
    weakly connected components. Results are cached per network instance.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network.
    
    Returns
    -------
    set[T]
        Set of cut-vertices.
    
    Examples
    --------
    >>> net = DirectedPhyNetwork(edges=[(1, 2), (2, 3), (2, 4)], nodes=[(3, {'label': 'A'}), (4, {'label': 'B'})])
    >>> vertices = cut_vertices(net)
    >>> 2 in vertices
    True
    >>> 1 in vertices
    False
    
    Notes
    -----
    Results are cached using LRU cache with maxsize=128.
    """
    return graph_cut_vertices(network._graph, data=False)


@lru_cache(maxsize=128)
def omnians(network: DirectedPhyNetwork) -> set[T]:
    """
    Find all omnian nodes in a directed phylogenetic network.
    
    An omnian is an internal node (non-leaf) where all of its children are hybrid nodes.
    
    Parameters
    ----------
    network : DirectedPhyNetwork[T]
        The directed phylogenetic network.
    
    Returns
    -------
    set[T]
        Set of omnian node identifiers.
    
    Warns
    -----
    PhyloZooWarning
        If the network contains parallel edges, as omnians are not defined for
        networks with parallel edges in the original paper. Behavior may be unexpected.
    
    Notes
    -----
    This function is based on the definition from:
    
    Jetten, Laura, and Leo van Iersel. "Nonbinary tree-based phylogenetic networks."
    IEEE/ACM transactions on computational biology and bioinformatics 15.1 (2016): 205-217.
    
    Examples
    --------
    >>> # Network with omnians (nodes 5 and 8 both have all children as hybrids)
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         (7, 5), (7, 8), (7, 9),  # Root to tree nodes
    ...         (5, 4), (5, 6),  # Node 5 to hybrid nodes 4 and 6
    ...         (8, 4), (8, 6),  # Node 8 to hybrid nodes 4 and 6
    ...         (9, 4), (9, 6),  # Node 9 to hybrid nodes 4 and 6
    ...         (4, 1), (6, 2)   # Hybrids to leaves
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> omnians(net)
    {5, 8, 9}
    
    >>> # Network with no omnians
    >>> net2 = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> omnians(net2)
    set()
    """
    # Check for parallel edges and warn (import here to avoid circular import)
    from .classifications import has_parallel_edges
    
    # Single-node network has no omnians
    if network.number_of_nodes() == 1:
        return set()
    
    # Parallel edges warn
    if has_parallel_edges(network):
        warnings.warn(
            "Network contains parallel edges. Omnians are not defined for networks "
            "with parallel edges in the original paper (Jetten & van Iersel, 2016). "
            "Behavior may be unexpected. Proceed with care.",
            PhyloZooWarning,
            stacklevel=2
        )
    
    # Get sets of leaves and hybrid nodes
    leaves = network.leaves
    hybrid_nodes = network.hybrid_nodes
    
    # Find omnians: internal nodes (not leaves) where all children are hybrid nodes
    omnian_set: set[T] = set()
    
    for node in network.internal_nodes:      
        # Get all children of this node
        children = list(network.children(node))

        # Check if all children are hybrid nodes
        if all(child in hybrid_nodes for child in children):
            omnian_set.add(node)
    
    return omnian_set

