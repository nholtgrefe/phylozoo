"""
Network features module.

This module provides functions to extract and identify features of directed
phylogenetic networks (e.g., LSA node, blobs, omnians, etc.).
"""

from collections import deque
from typing import TypeVar

import networkx as nx

from . import DirectedPhyNetwork

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

