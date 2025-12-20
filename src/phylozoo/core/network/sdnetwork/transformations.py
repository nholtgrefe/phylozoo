"""
Network transformations module.

This module provides functions to transform semi-directed and mixed phylogenetic networks
(e.g., suppress degree-2 nodes/blobs, identify parallel edges, etc.).
"""

from typing import Any

from . import SemiDirectedPhyNetwork
from ...primitives.m_multigraph.transformations import (
    identify_parallel_edge as mm_identify_parallel_edge,
    suppress_degree2_node as mm_suppress_degree2_node,
)


def _merge_attrs_for_degree2_suppression_mixed(
    edge1_data: dict[str, Any],
    edge2_data: dict[str, Any],
) -> dict[str, Any]:
    """
    Merge edge attributes for suppressing a degree-2 node in a mixed network.
    
    Special handling:
    - branch_length: If both edges have branch_length, sum them. Otherwise, use the one that has it.
    - gamma: If edge2 has gamma, use it. Otherwise, don't include gamma.
    - All other attributes (bootstrap, etc.) are removed.
    
    Parameters
    ----------
    edge1_data : dict[str, Any]
        Attributes of the first edge.
    edge2_data : dict[str, Any]
        Attributes of the second edge.
    
    Returns
    -------
    dict[str, Any]
        Merged attributes containing branch_length (if present) and gamma (if edge2 has it).
        All other attributes are removed.
    
    Notes
    -----
    This is an internal helper function for identify_parallel_edges.
    """
    merged: dict[str, Any] = {}
    
    # Handle branch_length: sum if both present, otherwise use the one that has it
    bl1 = edge1_data.get('branch_length') if edge1_data else None
    bl2 = edge2_data.get('branch_length') if edge2_data else None
    
    if bl1 is not None and bl2 is not None:
        merged['branch_length'] = bl1 + bl2
    elif bl1 is not None:
        merged['branch_length'] = bl1
    elif bl2 is not None:
        merged['branch_length'] = bl2
    
    # Handle gamma: use edge2's gamma if present
    gamma2 = edge2_data.get('gamma') if edge2_data else None
    if gamma2 is not None:
        merged['gamma'] = gamma2
    
    # All other attributes (bootstrap, etc.) are removed
    return merged


def _merge_attrs_for_parallel_identification_mixed(
    edges_data: list[dict[str, Any]],
) -> dict[str, Any]:
    """
    Merge edge attributes for identifying parallel edges in a mixed network.
    
    Special handling:
    - branch_length: Extract from the first edge (all should be same by validation).
    - gamma: Sum all gamma values from parallel edges.
    - All other attributes (bootstrap, etc.) are removed.
    
    Parameters
    ----------
    edges_data : list[dict[str, Any]]
        List of edge data dictionaries for parallel edges between the same pair of nodes.
    
    Returns
    -------
    dict[str, Any]
        Merged attributes containing branch_length (if present) and gamma (sum of all gammas).
        All other attributes are removed.
    
    Notes
    -----
    This is an internal helper function for identify_parallel_edges.
    All parallel edges should have the same branch_length by validation, so we
    just take it from the first edge.
    """
    merged: dict[str, Any] = {}
    
    if not edges_data:
        return merged
    
    # Extract branch_length from first edge (all should be same by validation)
    first_data = edges_data[0] or {}
    bl = first_data.get('branch_length')
    if bl is not None:
        merged['branch_length'] = bl
    
    # Sum gamma values from all parallel edges
    gamma_sum = 0.0
    has_gamma = False
    for edge_data in edges_data:
        if edge_data:
            gamma = edge_data.get('gamma')
            if gamma is not None:
                gamma_sum += gamma
                has_gamma = True
    
    if has_gamma:
        merged['gamma'] = gamma_sum
    
    # All other attributes (bootstrap, etc.) are removed
    return merged


def identify_parallel_edges(network: SemiDirectedPhyNetwork) -> SemiDirectedPhyNetwork:
    """
    Identify all parallel edges and suppress all degree-2 nodes exhaustively.
    
    This function iteratively:
    1. Identifies all parallel edges
    2. Suppresses all degree-2 nodes
    
    The process continues until no more changes occur, as suppression may create
    new parallel edges, and identification may create new degree-2 nodes.
    
    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network to transform.
    
    Returns
    -------
    SemiDirectedPhyNetwork
        A new network with all parallel edges identified and all degree-2 nodes suppressed.
    
    Notes
    -----
    - Branch lengths are preserved: summed when suppressing degree-2 nodes, kept from
      first edge when identifying parallel edges (all should be same by validation).
    - Gamma values: summed when identifying parallel edges, preserved from edge2 when
      suppressing degree-2 nodes (if edge2 has gamma, otherwise no gamma is included).
    - All other edge attributes (bootstrap, etc.) are removed.
    - Node labels and other node attributes are preserved.
    - Leaves are never suppressed.
    - Handles both directed and undirected parallel edges separately.
    
    Examples
    --------
    >>> net = SemiDirectedPhyNetwork(
    ...     directed_edges=[
    ...         {'u': 1, 'v': 2, 'branch_length': 0.5},
    ...         {'u': 1, 'v': 2, 'branch_length': 0.5}   # Parallel directed edges: node 2 is hybrid (in-degree 2)
    ...     ],
    ...     undirected_edges=[
    ...         {'u': 2, 'v': 3, 'branch_length': 0.3},  # Node 2 to tree node 3
    ...         {'u': 3, 'v': 4, 'branch_length': 0.2},  # Node 3 to leaf
    ...         {'u': 3, 'v': 5, 'branch_length': 0.1},  # Node 3 to leaf (keeps node 3 valid)
    ...         {'u': 1, 'v': 6, 'branch_length': 0.1}   # Root to another leaf
    ...     ],
    ...     nodes=[(4, {'label': 'A'}), (5, {'label': 'B'}), (6, {'label': 'C'})]
    ... )
    >>> result = identify_parallel_edges(net)
    >>> # Step 1: Parallel directed edges 1->2 identified -> node 2 becomes degree-2 (in-degree 1, undirected out)
    >>> # Step 2: Degree-2 node 2 suppressed (directed_in + undirected -> undirected) -> undirected edge 1-3 with branch_length=0.8 (0.5+0.3)
    >>> # Result: undirected edges 1-3 (0.8), 3-4 (0.2), 3-5 (0.1), 1-6 (0.1)
    """
    # Handle empty and single-node networks
    if network.number_of_nodes() == 0:
        return network.copy()
    
    if network.number_of_nodes() == 1:
        return network.copy()
    
    # Create a copy of the internal graph to modify
    from ...primitives.m_multigraph import MixedMultiGraph
    
    # Extract all edges from the network
    directed_edges_list: list[dict[str, Any]] = []
    undirected_edges_list: list[dict[str, Any]] = []
    
    for u, v, key, data in network._graph.directed_edges_iter(keys=True, data=True):
        edge_dict: dict[str, Any] = {'u': u, 'v': v}
        if key != 0:
            edge_dict['key'] = key
        edge_dict.update(data or {})
        directed_edges_list.append(edge_dict)
    
    for u, v, key, data in network._graph.undirected_edges_iter(keys=True, data=True):
        edge_dict: dict[str, Any] = {'u': u, 'v': v}
        if key != 0:
            edge_dict['key'] = key
        edge_dict.update(data or {})
        undirected_edges_list.append(edge_dict)
    
    # Create a working graph copy
    working_graph = MixedMultiGraph(
        directed_edges=directed_edges_list,
        undirected_edges=undirected_edges_list
    )
    
    # Get original leaves (these should never be suppressed)
    original_leaves = network.leaves
    
    # Iterate until no more changes occur
    max_iterations = 1000  # Safety limit
    iteration = 0
    
    while iteration < max_iterations:
        iteration += 1
        changes_made = False
        
        # Step 1: Find and identify all parallel directed edges
        parallel_directed: list[tuple[Any, Any]] = []
        for u, v in working_graph._directed.edges():
            if working_graph._directed.number_of_edges(u, v) > 1:
                parallel_directed.append((u, v))
        
        for u, v in parallel_directed:
            # Skip if nodes no longer exist
            if u not in working_graph.nodes() or v not in working_graph.nodes():
                continue
            
            # Skip if no longer parallel
            if working_graph._directed.number_of_edges(u, v) <= 1:
                continue
            
            # Collect edge data for all parallel edges
            edges_dict = working_graph._directed[u].get(v, {})
            if len(edges_dict) <= 1:
                continue
            
            edges_data = [edges_dict[key] for key in sorted(edges_dict.keys())]
            
            # Merge attributes using helper function
            merged_attrs = _merge_attrs_for_parallel_identification_mixed(edges_data)
            
            # Identify parallel edges
            mm_identify_parallel_edge(working_graph, u, v, merged_attrs=merged_attrs)
            changes_made = True
        
        # Step 2: Find and identify all parallel undirected edges
        parallel_undirected: list[tuple[Any, Any]] = []
        for u, v in working_graph._undirected.edges():
            if working_graph._undirected.number_of_edges(u, v) > 1:
                parallel_undirected.append((u, v))
        
        for u, v in parallel_undirected:
            # Skip if nodes no longer exist
            if u not in working_graph.nodes() or v not in working_graph.nodes():
                continue
            
            # Skip if no longer parallel
            if working_graph._undirected.number_of_edges(u, v) <= 1:
                continue
            
            # Collect edge data for all parallel edges (check both directions)
            edges_dict = working_graph._undirected[u].get(v, {})
            if not edges_dict:
                edges_dict = working_graph._undirected[v].get(u, {})
            
            if len(edges_dict) <= 1:
                continue
            
            edges_data = [edges_dict[key] for key in sorted(edges_dict.keys())]
            
            # Merge attributes using helper function
            merged_attrs = _merge_attrs_for_parallel_identification_mixed(edges_data)
            
            # Identify parallel edges
            mm_identify_parallel_edge(working_graph, u, v, merged_attrs=merged_attrs)
            changes_made = True
        
        # Step 3: Find and suppress all degree-2 nodes (excluding leaves)
        degree2_nodes: list[Any] = []
        for node in working_graph.nodes():
            if node not in original_leaves and working_graph.degree(node) == 2:
                degree2_nodes.append(node)
        
        for node in degree2_nodes:
            # Skip if node no longer exists or degree changed
            if node not in working_graph.nodes():
                continue
            
            if working_graph.degree(node) != 2:
                continue
            
            # Collect incident edges
            undirected_edges = list(working_graph.incident_undirected_edges(node, keys=True, data=True))
            directed_in = list(working_graph.incident_parent_edges(node, keys=True, data=True))
            directed_out = list(working_graph.incident_child_edges(node, keys=True, data=True))
            
            # Should have exactly 2 edges total for degree-2 node
            if len(undirected_edges) + len(directed_in) + len(directed_out) != 2:
                continue
            
            # Build edges_info
            edges_info: list[tuple[dict[str, Any], bool]] = []
            for u, v, key, data in directed_in:
                edges_info.append((data or {}, True))
            for u, v, key, data in directed_out:
                edges_info.append((data or {}, True))
            for u, v, key, data in undirected_edges:
                edges_info.append((data or {}, False))
            
            if len(edges_info) != 2:
                continue
            
            (d1, dir1), (d2, dir2) = edges_info
            
            # Merge attributes using helper function
            merged_attrs = _merge_attrs_for_degree2_suppression_mixed(d1, d2)
            
            # Suppress degree-2 node
            mm_suppress_degree2_node(working_graph, node, merged_attrs=merged_attrs)
            changes_made = True
        
        # If no changes were made, we're done
        if not changes_made:
            break
    
    if iteration >= max_iterations:
        raise RuntimeError(
            "identify_parallel_edges exceeded maximum iterations. "
            "This may indicate an infinite loop or a bug."
        )
    
    # Extract edges and nodes from the modified graph
    new_directed_edges: list[dict[str, Any]] = []
    for u, v, key, data in working_graph.directed_edges_iter(keys=True, data=True):
        edge_dict: dict[str, Any] = {'u': u, 'v': v}
        if key != 0:
            edge_dict['key'] = key
        edge_dict.update(data or {})
        new_directed_edges.append(edge_dict)
    
    new_undirected_edges: list[dict[str, Any]] = []
    for u, v, key, data in working_graph.undirected_edges_iter(keys=True, data=True):
        edge_dict: dict[str, Any] = {'u': u, 'v': v}
        if key != 0:
            edge_dict['key'] = key
        edge_dict.update(data or {})
        new_undirected_edges.append(edge_dict)
    
    # Preserve node labels
    new_nodes: list[Any | tuple[Any, dict[str, Any]]] = []
    for node in working_graph.nodes():
        label = network.get_label(node)
        if label is not None:
            new_nodes.append((node, {'label': label}))
    
    # Create and return new network
    return SemiDirectedPhyNetwork(
        directed_edges=new_directed_edges,
        undirected_edges=new_undirected_edges,
        nodes=new_nodes if new_nodes else None
    )

