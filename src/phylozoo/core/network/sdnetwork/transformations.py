"""
Network transformations module.

This module provides functions to transform semi-directed and mixed phylogenetic networks
(e.g., suppress degree-2 nodes/blobs, identify parallel edges, etc.).
"""

from typing import Any, TypeVar

from phylozoo.utils.exceptions import PhyloZooAlgorithmError

from . import SemiDirectedPhyNetwork
from .base import MixedPhyNetwork
from .features import k_blobs
from .conversions import sdnetwork_from_graph
from ._utils import (
    _merge_attrs_for_degree2_suppression_mixed,
    _merge_attrs_for_parallel_identification_mixed,
    _suppress_deg2_nodes,
)
from ...primitives.m_multigraph.transformations import (
    identify_parallel_edge as mm_identify_parallel_edge,
    identify_vertices as mm_identify_vertices,
)

T = TypeVar('T')


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
    
    Raises
    ------
    PhyloZooAlgorithmError
        If the algorithm exceeds the maximum number of iterations.
        
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
    
    # Create a working graph copy (preserves all node and edge attributes)
    working_graph = network._graph.copy()
    
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
        
        # Step 3: Suppress all degree-2 nodes (excluding leaves)
        nodes_before = set(working_graph.nodes())
        _suppress_deg2_nodes(working_graph)
        if set(working_graph.nodes()) != nodes_before:
            changes_made = True
        
        # If no changes were made, we're done
        if not changes_made:
            break
    
    if iteration >= max_iterations:
        raise PhyloZooAlgorithmError(
            "identify_parallel_edges exceeded maximum iterations. "
            "This may indicate an infinite loop or a bug."
        )
    
    # Create and return new network from the modified graph
    return sdnetwork_from_graph(working_graph, network_type='semi-directed')


def suppress_2_blobs(network: MixedPhyNetwork) -> MixedPhyNetwork:
    """
    Suppress all 2-blobs in the network.
    
    A 2-blob is a blob with exactly 2 incident edges. This function:
    1. Finds all 2-blobs using k_blobs
    2. For each 2-blob, identifies all vertices in the blob with the first vertex (creating a degree-2 node), then suppresses the degree-2 node using proper attribute merging
    3. Returns a new validated network
    
    Parameters
    ----------
    network : MixedPhyNetwork
        The mixed phylogenetic network to transform.
    
    Returns
    -------
    MixedPhyNetwork
        A new network with all 2-blobs suppressed.
        The network is validated before being returned.

    Notes
    -----
    - After identifying vertices in a 2-blob, the kept vertex becomes degree-2 and is then suppressed
    - Edge attributes (branch_length, gamma) are properly merged during suppression
    - The function works on a copy of the graph, so the original network is not modified
    - For mixed networks, incident edges can be directed or undirected, and the suppression handles both types correctly
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(1, 2), (2, 3), (3, 4), (4, 5), (1, 6)],
    ...     nodes=[(5, {'label': 'A'}), (6, {'label': 'B'})]
    ... )
    >>> result = suppress_2_blobs(net)
    >>> result.validate()  # Should not raise
    """
    # Work on a copy of the graph to avoid modifying the original
    working_graph = network._graph.copy()
    
    # Find all 2-blobs
    two_blobs = k_blobs(network, k=2, trivial=False, leaves=False)
    
    # Process each 2-blob
    for blob in two_blobs:
        # Convert blob set to sorted list, keep first vertex
        blob_list = sorted(blob)
        if len(blob_list) < 2:
            continue  # Skip trivial blobs (shouldn't happen with trivial=False)
        
        first_vertex = blob_list[0]
        
        # Identify all vertices in blob with the first vertex
        # This modifies working_graph in place
        mm_identify_vertices(working_graph, blob_list)
        
        # After identification, the first_vertex should be degree-2
        # (since 2-blobs have exactly 2 incident edges)
        if working_graph.degree(first_vertex) != 2:
            # This shouldn't happen, but skip if it does
            continue
        
        # Collect incident edges to merge attributes
        # Order: directed_in first, then undirected, then directed_out
        edges_info: list[dict[str, Any]] = []
        
        # Collect in order: directed_in, undirected, directed_out
        for u, v, key, data in working_graph.incident_parent_edges(first_vertex, keys=True, data=True):
            edges_info.append(data or {})
        for u, v, key, data in working_graph.incident_undirected_edges(first_vertex, keys=True, data=True):
            edges_info.append(data or {})
        for u, v, key, data in working_graph.incident_child_edges(first_vertex, keys=True, data=True):
            edges_info.append(data or {})
        
        # We should have exactly 2 incident edges
        if len(edges_info) != 2:
            # Invalid configuration, skip
            continue
        
        # Merge attributes using the helper function
        merged_attrs = _merge_attrs_for_degree2_suppression_mixed(edges_info[0], edges_info[1])
        
        # Suppress the degree-2 node with merged attributes
        # Note: We suppress a single specific node here, not all degree-2 nodes
        from ...primitives.m_multigraph.transformations import suppress_degree2_node as mm_suppress_degree2_node
        mm_suppress_degree2_node(working_graph, first_vertex, merged_attrs=merged_attrs)
    
    # Create and return new network from the modified graph (will be validated)
    # Return same type as input
    network_type = 'semi-directed' if isinstance(network, SemiDirectedPhyNetwork) else 'mixed'
    return sdnetwork_from_graph(working_graph, network_type=network_type)

