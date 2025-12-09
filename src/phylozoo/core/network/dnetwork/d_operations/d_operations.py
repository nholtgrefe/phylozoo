"""
Operations for directed phylogenetic networks.

This module provides operations for DirectedPhyNetwork instances.
"""

from collections import deque
from typing import Any, Dict, List, Optional, Set, Tuple, TypeVar

import networkx as nx

from ..d_phynetwork import DirectedPhyNetwork
from ..d_classifications import is_LSA_network
from ....primitives.m_multigraph.mm_graph import MixedMultiGraph
from ....primitives.m_multigraph.mm_operations import suppress_degree2_node
from ...sdnetwork.sd_phynetwork import SemiDirectedPhyNetwork

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
    depths: Dict[T, int] = {}
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
    # Use dict comprehension but avoid calling get_label twice
    new_taxa: Dict[Any, str] = {}
    for leaf in network.leaves:
        label = network.get_label(leaf)
        if label is not None:
            new_taxa[leaf] = label
    
    # Collect internal node labels for nodes that are still in the network (dict comprehension)
    new_internal_labels: Dict[Any, str] = {
        node: network.get_label(node)
        for node in descendants
        if node not in network.leaves
        and node != lsa_node
        and network.get_label(node) is not None
    }
    
    # Create new network
    return DirectedPhyNetwork(
        edges=new_edges,
        taxa=new_taxa,
        internal_node_labels=new_internal_labels if new_internal_labels else None
    )


def _merge_edge_attributes_for_suppression(
    edge1_data: Dict[str, Any],
    edge2_data: Dict[str, Any],
    edge1_is_directed: bool,
    edge2_is_directed: bool,
) -> Dict[str, Any]:
    """
    Merge edge attributes for suppressing a degree-2 node.
    
    Special handling for 'gamma' and 'branch_length':
    - Gamma: If one edge is directed and has gamma, use that gamma. If both have gamma,
      multiply them.
    - Branch length: If both have branch_length, sum them.
    - Other attributes: Use edge1_data values, with edge2_data overriding.
    
    Parameters
    ----------
    edge1_data : Dict[str, Any]
        Attributes of the first edge.
    edge2_data : Dict[str, Any]
        Attributes of the second edge.
    edge1_is_directed : bool
        Whether the first edge is directed.
    edge2_is_directed : bool
        Whether the second edge is directed.
    
    Returns
    -------
    Dict[str, Any]
        Merged attributes with special handling for gamma and branch_length.
    """
    merged = {}
    
    # Start with edge1_data
    if edge1_data:
        merged.update(edge1_data)
    
    # Handle gamma specially
    gamma1 = edge1_data.get('gamma') if edge1_data else None
    gamma2 = edge2_data.get('gamma') if edge2_data else None
    
    if gamma1 is not None and gamma2 is not None:
        # Both have gamma: multiply them
        merged['gamma'] = gamma1 * gamma2
    elif gamma1 is not None and edge1_is_directed:
        # Only edge1 has gamma and it's directed: use it
        merged['gamma'] = gamma1
    elif gamma2 is not None and edge2_is_directed:
        # Only edge2 has gamma and it's directed: use it
        merged['gamma'] = gamma2
    elif gamma1 is not None:
        # Only edge1 has gamma (but not directed): use it
        merged['gamma'] = gamma1
    elif gamma2 is not None:
        # Only edge2 has gamma (but not directed): use it
        merged['gamma'] = gamma2
    
    # Handle branch_length specially
    bl1 = edge1_data.get('branch_length') if edge1_data else None
    bl2 = edge2_data.get('branch_length') if edge2_data else None
    
    if bl1 is not None and bl2 is not None:
        # Both have branch_length: sum them
        merged['branch_length'] = bl1 + bl2
    elif bl1 is not None:
        merged['branch_length'] = bl1
    elif bl2 is not None:
        merged['branch_length'] = bl2
    
    # Override with other edge2_data attributes (except gamma and branch_length which we handled)
    if edge2_data:
        for key, value in edge2_data.items():
            if key not in ('gamma', 'branch_length'):
                merged[key] = value
    
    return merged


def to_sd_network(d_network: DirectedPhyNetwork) -> SemiDirectedPhyNetwork:
    """
    Convert a DirectedPhyNetwork to a SemiDirectedPhyNetwork.
    
    Steps:
    1. If the directed network is not an LSA network, replace it by its LSA-network.
    2. Undirect all non-hybrid edges; hybrid edges remain directed.
    3. Suppress any degree-2 node (this stems from a degree-2 root). Suppression may
       create parallel edges. Suppression connects the two neighbors directly:
       - undirected+undirected -> undirected
       - directed+directed (u->x, x->v) -> directed (u->v)
       - directed into x and undirected out (u->x, x—v) -> directed (u->v)
       - undirected into x and directed out (u—x, x->v) -> directed (u->v)
    
    Parameters
    ----------
    d_network : DirectedPhyNetwork
        The directed phylogenetic network to convert.
    
    Returns
    -------
    SemiDirectedPhyNetwork
        The corresponding semi-directed phylogenetic network.
    """
    # 1) Ensure LSA network
    working = d_network if is_LSA_network(d_network) else to_LSA_network(d_network)
    
    # Empty network shortcut
    if working.number_of_nodes() == 0:
        return SemiDirectedPhyNetwork(directed_edges=[], undirected_edges=[], taxa={}, internal_node_labels=None)
    
    # 2) Separate hybrid (directed) vs tree (to be undirected) edges
    hybrid_edge_set = set(working.hybrid_edges)
    directed_edges: List[Dict[str, Any]] = []
    undirected_edges: List[Dict[str, Any]] = []
    
    for u, v, key, data in working._graph.edges(keys=True, data=True):
        edge_dict: Dict[str, Any] = {"u": u, "v": v}
        if key != 0:
            edge_dict["key"] = key
        edge_dict.update(data)
        if (u, v) in hybrid_edge_set:
            directed_edges.append(edge_dict)
        else:
            undirected_edges.append(edge_dict)
    
    # 3) Build a mixed graph to allow suppression
    mixed = MixedMultiGraph(directed_edges=directed_edges, undirected_edges=undirected_edges)
    
    # Track which nodes are suppressed (leaves are never suppressed, only internal nodes)
    suppressed_nodes: Set[Any] = set()
    
    # Maintain a set of degree-2 nodes for efficient lookup and updates
    # Initialize by scanning all nodes once
    degree2_nodes: Set[Any] = {node for node in mixed.nodes() if mixed.degree(node) == 2}
    
    # Suppress all degree-2 nodes (iteratively, as suppression may create new ones)
    # Before suppression, merge attributes with special handling for gamma and branch_length
    while degree2_nodes:
        # Process each degree-2 node (convert to list to avoid modification during iteration)
        nodes_to_process = list(degree2_nodes)
        degree2_nodes.clear()
        
        for node in nodes_to_process:
            # Skip if node was already suppressed (can happen if multiple nodes point to same neighbor)
            if node not in mixed.nodes():
                continue
            
            # We already know degree is 2, so collect the two incident edges directly
            # Collect incident edges to merge attributes using public API
            undirected_edges = list(mixed.incident_undirected_edges(node, keys=True, data=True))
            directed_in = list(mixed.incident_parent_edges(node, keys=True, data=True))
            directed_out = list(mixed.incident_child_edges(node, keys=True, data=True))
            
            # Since we know degree is 2, we expect exactly 2 edges total
            # Build edges_info directly (we know there are exactly 2)
            edges_info: List[Tuple[Dict[str, Any], bool]] = []
            for u, v, key, data in directed_in:
                edges_info.append((data, True))  # (data, is_directed)
            for u, v, key, data in directed_out:
                edges_info.append((data, True))  # (data, is_directed)
            for u, v, key, data in undirected_edges:
                edges_info.append((data, False))  # (data, is_directed)
            
            # Defensive check: if degree changed, skip (shouldn't happen, but be safe)
            if len(edges_info) != 2:
                continue
            
            (d1, dir1), (d2, dir2) = edges_info
            
            # Merge attributes with special handling for gamma and branch_length
            merged_attrs = _merge_edge_attributes_for_suppression(d1, d2, dir1, dir2)
            
            # Track that this node will be suppressed
            suppressed_nodes.add(node)
            
            # Get neighbors before suppression to update degree2_nodes set
            neighbors_before = set()
            for u, v, key, data in directed_in:
                neighbors_before.add(u)
            for u, v, key, data in directed_out:
                neighbors_before.add(v)
            for u, v, key, data in undirected_edges:
                neighbors_before.add(v if u == node else u)
            
            # Use suppress_degree2_node with merged attributes
            suppress_degree2_node(mixed, node, merged_attrs=merged_attrs)
            
            # After suppression, check if any neighbors became degree-2
            for neighbor in neighbors_before:
                if neighbor in mixed.nodes() and mixed.degree(neighbor) == 2:
                    degree2_nodes.add(neighbor)
    
    # Taxa stay the same (leaves are never suppressed, only internal nodes)
    # Reuse the taxa mapping from the working network (dict comprehension for efficiency)
    new_taxa: Dict[Any, str] = {
        leaf: working.get_label(leaf)
        for leaf in working.leaves
        if working.get_label(leaf) is not None
    }
    
    # Internal labels: remove any nodes that were suppressed (dict comprehension for efficiency)
    new_internal_labels: Dict[Any, str] = {
        node: label
        for node, label in working._node_to_label.items()
        if node not in working.leaves and node not in suppressed_nodes
    }
    
    # Rebuild directed/undirected edges from mixed graph to pass to SD network
    # Use public API methods
    final_directed = []
    for u, v, k, d in mixed.directed_edges_iter(keys=True, data=True):
        edge_dict = {"u": u, "v": v}
        if k != 0:
            edge_dict["key"] = k
        edge_dict.update(d)
        final_directed.append(edge_dict)
    final_undirected = []
    for u, v, k, d in mixed.undirected_edges_iter(keys=True, data=True):
        edge_dict = {"u": u, "v": v}
        if k != 0:
            edge_dict["key"] = k
        edge_dict.update(d)
        final_undirected.append(edge_dict)
    
    return SemiDirectedPhyNetwork(
        directed_edges=final_directed,
        undirected_edges=final_undirected,
        taxa=new_taxa,
        internal_node_labels=new_internal_labels if new_internal_labels else None,
    )
