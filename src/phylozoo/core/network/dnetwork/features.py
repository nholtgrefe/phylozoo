"""
Network features module.

This module provides functions to extract and identify features of directed
phylogenetic networks (e.g., LSA node, blobs, omnians, etc.).
"""

from collections import deque
from typing import Iterator, TypeVar

import networkx as nx

from ...primitives.d_multigraph.features import biconnected_components
from .base import DirectedPhyNetwork

T = TypeVar('T')


def find_lsa_node(network: DirectedPhyNetwork) -> T:
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
    ValueError
        If the network is empty or has no leaves.
    
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
    >>> find_lsa_node(net)
    4
    >>> # In a simple tree, the LSA is just the root
    >>> net2 = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
    >>> find_lsa_node(net2)
    3
    """
    if network.number_of_nodes() == 0:
        raise ValueError("Cannot find LSA node in empty network")
    
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
    
    # Find all ancestors of each leaf
    # Use NetworkX's ancestors function on the underlying graph
    dag = network._graph._graph
    
    # For each leaf, find all its ancestors
    ancestor_sets = []
    for leaf in leaves:
        ancestors = set(nx.ancestors(dag, leaf))
        ancestors.add(leaf)  # Include the leaf itself
        ancestor_sets.append(ancestors)
    
    # Find intersection: nodes that are ancestors of all leaves
    common_ancestors = ancestor_sets[0]
    for ancestor_set in ancestor_sets[1:]:
        common_ancestors &= ancestor_set
    
    if not common_ancestors:
        # This shouldn't happen in a valid network, but handle it
        raise ValueError("No common ancestor found for all leaves")
    
    # Find the node with maximum depth (distance from root)
    root = network.root_node
    
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
    
    # Find the deepest node among common ancestors using max() with key function
    return max(common_ancestors, key=lambda node: depths.get(node, 0))


def blobs(
    network: DirectedPhyNetwork,
    trivial: bool = True,
    leaves: bool = True,
) -> Iterator[set[T]]:
    """
    Yield blobs of the network.
    
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
    
    Yields
    ------
    set[T]
        Sets of nodes forming each blob.
    
    Raises
    ------
    ValueError
        If `trivial=False` and `leaves=True` (this combination is not possible
        since leaves are single-node components).
    
    Notes
    -----
    Blobs are computed as:
    - All biconnected components with at least 3 nodes
    - Biconnected components with exactly 2 nodes that have parallel edges
    - All remaining nodes (trivial blobs, including leaves)
    
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
    >>> len(list(blobs(net, trivial=False, leaves=False)))
    1
    >>> # Filtering: exclude blobs containing only leaves
    >>> len(list(blobs(net, leaves=False)))
    2
    """
    # Check for invalid parameter combination
    if not trivial and leaves:
        raise ValueError(
            "Cannot have trivial=False and leaves=True: leaves are single-node "
            "components, so excluding trivial components would exclude all leaves."
        )
    
    leaves_set = network.leaves
    graph = network._graph
    
    # Collect nodes that are in non-trivial blobs
    visited: set[T] = set()
    
    # Process biconnected components
    for blob in biconnected_components(graph):
        blob_set = set(blob)
        
        # Keep blobs with >= 3 nodes
        if len(blob_set) >= 3:
            visited.update(blob_set)
            # Filter blobs containing only leaves
            if not leaves and blob_set.issubset(leaves_set):
                continue
            yield blob_set
            continue
        
        # Keep 2-node blobs only if they have parallel edges
        if len(blob_set) == 2:
            u, v = list(blob_set)
            # Check for parallel edges in the underlying graph
            # Count edges in both directions (u->v and v->u)
            num_edges = 0
            if u in graph._graph and v in graph._graph[u]:
                num_edges += len(graph._graph[u][v])
            if v in graph._graph and u in graph._graph[v]:
                num_edges += len(graph._graph[v][u])
            
            if num_edges > 1:  # Parallel edges exist
                visited.update(blob_set)
                # Filter blobs containing only leaves
                yield blob_set
    
    # Add remaining nodes as trivial blobs
    for node in graph.nodes():
        if node not in visited:
            if node not in leaves_set and trivial:
                yield {node}
            elif node in leaves_set and leaves:
                yield {node}

