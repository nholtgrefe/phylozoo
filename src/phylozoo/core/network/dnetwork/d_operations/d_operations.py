"""
Operations for directed phylogenetic networks.

This module provides operations for DirectedPhyNetwork instances.
"""

from typing import TypeVar

from ..d_phynetwork import DirectedPhyNetwork

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
    ...     taxa={1: "A", 2: "B"}
    ... )
    >>> find_lsa_node(net)
    4
    >>> # In a simple tree, the LSA is just the root
    >>> net2 = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"})
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
        if len(parents) == 1:
            return parents[0]
        # This shouldn't happen in valid networks, but handle it
        return parents[0]
    
    # Find all ancestors of each leaf
    # Use NetworkX's ancestors function on the underlying graph
    import networkx as nx
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
    
    # Compute depths using BFS from root
    depths: dict[T, int] = {}
    queue = [root]
    depths[root] = 0
    
    while queue:
        current = queue.pop(0)
        for child in network.children(current):
            if child not in depths:
                depths[child] = depths[current] + 1
                queue.append(child)
    
    # Find the deepest node among common ancestors
    lsa_node = root
    max_depth = depths.get(root, 0)
    
    for node in common_ancestors:
        node_depth = depths.get(node, 0)
        if node_depth > max_depth:
            max_depth = node_depth
            lsa_node = node
    
    return lsa_node


def to_LSA_network(network: DirectedPhyNetwork) -> DirectedPhyNetwork:
    """
    Create a new LSA-network by removing everything above the LSA node.
    
    This function finds the LSA (Least Stable Ancestor) node and creates a new network
    that contains only the LSA node and all nodes/edges below it. The LSA becomes the
    new root of the resulting network.
    
    Note: All branch lengths, bootstrap values, gamma values, and other edge attributes
    from edges above the LSA are removed (as those edges are removed). Edge attributes
    for edges below the LSA are preserved.
    
    Parameters
    ----------
    network : DirectedPhyNetwork[T]
        The directed phylogenetic network.
    
    Returns
    -------
    DirectedPhyNetwork[T]
        A new network with the LSA as the root.
    
    Raises
    ------
    ValueError
        If the network is empty or has no leaves.
    
    Examples
    --------
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         (7, 5), (7, 6),                 # root to tree nodes
    ...         (5, 4, 0), (5, 4, 1),           # parallel edges keep tree node 5 out-degree >= 2
    ...         (6, 4, 0), (6, 4, 1),           # parallel edges keep tree node 6 out-degree >= 2
    ...         (4, 10),                        # hybrid 4 (in-degree 4, out-degree 1) to tree node 10
    ...         (10, 1), (10, 2)                # tree node 10 splits to leaves
    ...     ],
    ...     taxa={1: "A", 2: "B"}
    ... )
    >>> lsa_net = to_LSA_network(net)
    >>> lsa_net.root_node
    4
    >>> sorted(lsa_net.leaves)
    [1, 2]
    """
    if network.number_of_nodes() == 0:
        raise ValueError("Cannot create LSA network from empty network")
    
    lsa_node = network.LSA_node
    
    # If LSA is the root, return a copy
    if lsa_node == network.root_node:
        return network.copy()
    
    # Find all nodes that are descendants of the LSA (including LSA itself)
    import networkx as nx
    dag = network._graph._graph
    descendants = set(nx.descendants(dag, lsa_node))
    descendants.add(lsa_node)
    
    # Collect all edges that are between nodes in the descendant set
    new_edges = []
    for u, v, key, data in network._graph.edges(keys=True, data=True):
        if u in descendants and v in descendants:
            # Reconstruct edge in the same format
            edge_dict = {'u': u, 'v': v}
            if key != 0:  # Only include key if it's not the default
                edge_dict['key'] = key
            # Copy all edge attributes
            edge_dict.update(data)
            new_edges.append(edge_dict)
    
    # LSA networks keep the same leaves and taxa as the original network
    new_taxa = {leaf: network.get_label(leaf) for leaf in network.leaves if network.get_label(leaf) is not None}
    
    # Collect internal node labels for nodes that are still in the network
    new_internal_labels = {}
    for node in descendants:
        if node not in network.leaves and node != lsa_node:
            label = network.get_label(node)
            if label is not None:
                new_internal_labels[node] = label
    
    # Create new network
    return DirectedPhyNetwork(
        edges=new_edges,
        taxa=new_taxa,
        internal_node_labels=new_internal_labels if new_internal_labels else None
    )
