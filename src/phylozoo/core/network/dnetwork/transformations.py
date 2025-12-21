"""
Network transformations module.

This module provides functions to transform directed phylogenetic networks
(e.g., suppress degree-2 nodes, convert to LSA network, convert to semi-directed network, etc.).
"""

from typing import Any

from . import DirectedPhyNetwork
from .classifications import is_lsa_network
from .features import k_blobs
from .io import dnetwork_from_dmgraph
from ...primitives.d_multigraph.transformations import (
    identify_parallel_edge as dm_identify_parallel_edge,
    identify_vertices as dm_identify_vertices,
    suppress_degree2_node as dm_suppress_degree2_node,
)
from ...primitives.m_multigraph import MixedMultiGraph
from ...primitives.m_multigraph.transformations import suppress_degree2_node
from ..sdnetwork import SemiDirectedPhyNetwork


def to_lsa_network(network: DirectedPhyNetwork) -> DirectedPhyNetwork:
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
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> lsa_net = to_lsa_network(net)
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
    
    # LSA networks keep the same leaves and labels as the original network
    # Build nodes list with labels in NetworkX-style format
    new_nodes: list[Any | tuple[Any | dict[str | str]]] = []
    for leaf in network.leaves:
        label = network.get_label(leaf)
        if label is not None:
            new_nodes.append((leaf, {'label': label}))
    
    # Collect internal node labels for nodes that are still in the network
    for node in descendants:
        if node not in network.leaves and node != lsa_node:
            label = network.get_label(node)
            if label is not None:
                new_nodes.append((node, {'label': label}))
    
    # Create new network
    return DirectedPhyNetwork(
        edges=new_edges,
        nodes=new_nodes if new_nodes else None
    )


def _merge_edge_attributes_for_suppression(
    edge1_data: dict[str, Any],
    edge2_data: dict[str, Any],
    edge1_is_directed: bool,
    edge2_is_directed: bool,
) -> dict[str, Any]:
    """
    Merge edge attributes for suppressing a degree-2 node.
    
    Special handling for 'gamma' and 'branch_length':
    - Gamma: If one edge is directed and has gamma, use that gamma. If both have gamma,
      multiply them.
    - Branch length: If both have branch_length, sum them.
    - Other attributes: Use edge1_data values, with edge2_data overriding.
    
    Parameters
    ----------
    edge1_data : dict[str, Any]
        Attributes of the first edge.
    edge2_data : dict[str, Any]
        Attributes of the second edge.
    edge1_is_directed : bool
        Whether the first edge is directed.
    edge2_is_directed : bool
        Whether the second edge is directed.
    
    Returns
    -------
    dict[str, Any]
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
       - directed into x and undirected out (u->x, x—v) -> undirected (u-v)
       - undirected into x and directed out (u—x, x->v) -> directed (u->v)
    
    Parameters
    ----------
    d_network : DirectedPhyNetwork
        The directed phylogenetic network to convert.
    
    Returns
    -------
    SemiDirectedPhyNetwork
        The corresponding semi-directed phylogenetic network.

    Examples
    --------
    >>> # Simple tree (no hybrids) - all edges become undirected
    >>> dnet = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> sdnet = to_sd_network(dnet)
    >>> sdnet.number_of_directed_edges()
    0
    >>> sdnet.number_of_undirected_edges()
    2

    >>> # Network with hybrids - hybrid edges remain directed
    >>> dnet = DirectedPhyNetwork(
    ...     edges=[
    ...         (4, 1), (4, 2),  # Tree edges from root
    ...         {'u': 1, 'v': 3, 'gamma': 0.6},  # Hybrid edge
    ...         {'u': 2, 'v': 3, 'gamma': 0.4}   # Hybrid edge
    ...     ],
    ...     nodes=[(3, {'label': 'C'})]
    ... )
    >>> sdnet = to_sd_network(dnet)
    >>> sdnet.number_of_directed_edges()  # Hybrid edges
    2
    >>> sdnet.number_of_undirected_edges()  # Tree edges
    2
    """
    # Single-leaf shortcut: the LSA of a single-leaf network is that leaf,
    # so the semi-directed network collapses to a single node with no edges.
    if len(d_network.leaves) == 1:
        leaf = next(iter(d_network.leaves))
        label = d_network.get_label(leaf)
        nodes_list = [(leaf, {'label': label})] if label is not None else None
        return SemiDirectedPhyNetwork(
            directed_edges=[],
            undirected_edges=[],
            nodes=nodes_list,
        )
    
    # Empty network shortcut
    if d_network.number_of_nodes() == 0:
        return SemiDirectedPhyNetwork(directed_edges=[], undirected_edges=[], nodes=None)
    
    # 1) Ensure LSA network
    working = d_network if is_lsa_network(d_network) else to_lsa_network(d_network)
    
    # Single-node network shortcut
    if working.number_of_nodes() == 1:
        node = next(iter(working._graph.nodes))
        label = working.get_label(node)
        nodes_list = [(node, {'label': label})] if label else None
        return SemiDirectedPhyNetwork(
            directed_edges=[],
            undirected_edges=[],
            nodes=nodes_list
        )
    
    # 2) Separate hybrid (directed) vs tree (to be undirected) edges
    hybrid_edge_set = working.hybrid_edges  # Now contains (u, v, key) tuples
    directed_edges: list[dict[str, Any]] = []
    undirected_edges: list[dict[str, Any]] = []
    
    for u, v, key, data in working._graph.edges(keys=True, data=True):
        edge_dict: dict[str, Any] = {"u": u, "v": v}
        if key != 0:
            edge_dict["key"] = key
        edge_dict.update(data)
        if (u, v, key) in hybrid_edge_set:
            directed_edges.append(edge_dict)
        else:
            undirected_edges.append(edge_dict)
    
    # 3) Build a mixed graph to allow suppression
    mixed = MixedMultiGraph(directed_edges=directed_edges, undirected_edges=undirected_edges)
    
    # Track which nodes are suppressed (leaves are never suppressed, only internal nodes)
    suppressed_nodes: set[Any] = set()
    
    # Maintain a set of degree-2 nodes for efficient lookup and updates
    # Initialize by scanning all nodes once
    degree2_nodes: set[Any] = {node for node in mixed.nodes() if mixed.degree(node) == 2}
    
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
            edges_info: list[tuple[dict[str, Any], bool]] = []
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
    
    # Build nodes list with labels (leaves are never suppressed, only internal nodes)
    # Leaves: all leaves from working network
    # Internal nodes: only those not suppressed
    nodes_list: list[Any | tuple[Any | dict[str | str]]] = []
    
    # Add leaves with labels
    for leaf in working.leaves:
        label = working.get_label(leaf)
        if label is not None:
            nodes_list.append((leaf, {'label': label}))
    
    # Add internal nodes with labels (excluding suppressed nodes)
    for node, label in working._node_to_label.items():
        if node not in working.leaves and node not in suppressed_nodes:
            nodes_list.append((node, {'label': label}))
    
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
        nodes=nodes_list if nodes_list else None,
    )


def _merge_attrs_for_degree2_suppression_directed(
    edge1_data: dict[str, Any],
    edge2_data: dict[str, Any],
) -> dict[str, Any]:
    """
    Merge edge attributes for suppressing a degree-2 node in a directed network.
    
    Special handling:
    - branch_length: If both edges have branch_length, sum them. Otherwise, use the one that has it.
    - gamma: If edge2 has gamma, use it. Otherwise, don't include gamma.
    - All other attributes (bootstrap, etc.) are removed.
    
    Parameters
    ----------
    edge1_data : dict[str, Any]
        Attributes of the first edge (uv).
    edge2_data : dict[str, Any]
        Attributes of the second edge (vw).
    
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


def _merge_attrs_for_parallel_identification_directed(
    edges_data: list[dict[str, Any]],
) -> dict[str, Any]:
    """
    Merge edge attributes for identifying parallel edges in a directed network.
    
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


def identify_parallel_edges(network: DirectedPhyNetwork) -> DirectedPhyNetwork:
    """
    Identify all parallel edges and suppress all degree-2 nodes exhaustively.
    
    This function iteratively:
    1. Identifies all parallel edges (removes all but one, keeping branch_length)
    2. Suppresses all degree-2 nodes (sums branch_lengths, removes other attributes)
    
    The process continues until no more changes occur, as suppression may create
    new parallel edges, and identification may create new degree-2 nodes.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to transform.
    
    Returns
    -------
    DirectedPhyNetwork
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
    
    Examples
    --------
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         {'u': 1, 'v': 2, 'branch_length': 0.5},
    ...         {'u': 1, 'v': 2, 'branch_length': 0.5},  # Parallel edges: node 2 is hybrid (in-degree 2)
    ...         {'u': 2, 'v': 3, 'branch_length': 0.3},  # Node 2 to tree node 3
    ...         {'u': 3, 'v': 4, 'branch_length': 0.2},  # Node 3 to leaf
    ...         {'u': 3, 'v': 5, 'branch_length': 0.1},  # Node 3 to leaf (keeps node 3 valid)
    ...         {'u': 1, 'v': 6, 'branch_length': 0.1}
    ...     ],
    ...     nodes=[(4, {'label': 'A'}), (5, {'label': 'B'}), (6, {'label': 'C'})]
    ... )
    >>> result = identify_parallel_edges(net)
    >>> # Step 1: Parallel edges 1->2 identified -> node 2 becomes degree-2 (in-degree 1, out-degree 1)
    >>> # Step 2: Degree-2 node 2 suppressed -> edge 1->3 with branch_length=0.8 (0.5+0.3)
    >>> # Result: edges 1->3 (0.8), 3->4 (0.2), 3->5 (0.1), 1->6 (0.1)
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
        
        # Step 1: Find and identify all parallel edges
        parallel_pairs: list[tuple[Any, Any]] = []
        for u, v in working_graph._graph.edges():
            if working_graph._graph.number_of_edges(u, v) > 1:
                parallel_pairs.append((u, v))
        
        for u, v in parallel_pairs:
            # Skip if nodes no longer exist (may have been removed)
            if u not in working_graph.nodes() or v not in working_graph.nodes():
                continue
            
            # Skip if no longer parallel (may have been identified by previous iteration)
            if working_graph._graph.number_of_edges(u, v) <= 1:
                continue
            
            # Collect edge data for all parallel edges
            edges_dict = working_graph._graph[u].get(v, {})
            if len(edges_dict) <= 1:
                continue
            
            edges_data = [edges_dict[key] for key in sorted(edges_dict.keys())]
            
            # Merge attributes using helper function
            merged_attrs = _merge_attrs_for_parallel_identification_directed(edges_data)
            
            # Identify parallel edges
            dm_identify_parallel_edge(working_graph, u, v, merged_attrs=merged_attrs)
            changes_made = True
        
        # Step 2: Find and suppress all degree-2 nodes (excluding leaves)
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
            incoming = list(working_graph.incident_parent_edges(node, keys=True, data=True))
            outgoing = list(working_graph.incident_child_edges(node, keys=True, data=True))
            
            # Should have exactly one incoming and one outgoing for degree-2 node
            if len(incoming) != 1 or len(outgoing) != 1:
                continue
            
            # Get edge data
            (u, _, k1, d1) = incoming[0]
            (_, v, k2, d2) = outgoing[0]
            
            # Merge attributes using helper function
            merged_attrs = _merge_attrs_for_degree2_suppression_directed(
                d1 or {}, d2 or {}
            )
            
            # Suppress degree-2 node
            dm_suppress_degree2_node(working_graph, node, merged_attrs=merged_attrs)
            changes_made = True
        
        # If no changes were made, we're done
        if not changes_made:
            break
    
    if iteration >= max_iterations:
        raise RuntimeError(
            "identify_parallel_edges exceeded maximum iterations. "
            "This may indicate an infinite loop or a bug."
        )
    
    # Create and return new network from the modified graph
    return dnetwork_from_dmgraph(working_graph)


def suppress_2_blobs(network: DirectedPhyNetwork) -> DirectedPhyNetwork:
    """
    Suppress all 2-blobs in the network.
    
    A 2-blob is a blob with exactly 2 incident edges. This function:
    1. Finds all 2-blobs using k_blobs
    2. For each 2-blob (except those containing the root):
       - Identifies all vertices in the blob with the first vertex
       - This creates a degree-2 node
       - Suppresses the degree-2 node using proper attribute merging
    3. Returns a new validated network
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to transform.
    
    Returns
    -------
    DirectedPhyNetwork
        A new network with all 2-blobs suppressed (except root-containing ones).
        The network is validated before being returned.
    
    Raises
    ------
    ValueError
        If the transformation creates an invalid network structure.
    
    Notes
    -----
    - 2-blobs containing the root node are skipped (special case for directed networks)
    - After identifying vertices in a 2-blob, the kept vertex becomes degree-2
      and is then suppressed
    - Edge attributes (branch_length, gamma) are properly merged during suppression
    - The function works on a copy of the graph, so the original network is not modified
    
    Examples
    --------
    >>> net = DirectedPhyNetwork(
    ...     edges=[(1, 2), (2, 3), (3, 4), (4, 5), (1, 6)],
    ...     nodes=[(5, {'label': 'A'}), (6, {'label': 'B'})]
    ... )
    >>> result = suppress_2_blobs(net)
    >>> result.validate()  # Should not raise
    """
    # Work on a copy of the graph to avoid modifying the original
    working_graph = network._graph.copy()
    
    # Find all 2-blobs
    two_blobs = k_blobs(network, k=2, trivial=False, leaves=False)
    
    # Get root node for checking
    root = network.root_node
    
    # Process each 2-blob
    for blob in two_blobs:
        # Skip 2-blobs containing the root (special case for directed networks)
        if root in blob:
            continue
        
        # Convert blob set to sorted list, keep first vertex
        blob_list = sorted(blob)
        if len(blob_list) < 2:
            continue  # Skip trivial blobs (shouldn't happen with trivial=False)
        
        first_vertex = blob_list[0]
        
        # Identify all vertices in blob with the first vertex
        # This modifies working_graph in place
        dm_identify_vertices(working_graph, blob_list)
        
        # After identification, the first_vertex should be degree-2
        # (since 2-blobs have exactly 2 incident edges)
        if working_graph.degree(first_vertex) != 2:
            # This shouldn't happen, but skip if it does
            continue
        
        # Collect incident edges to merge attributes
        directed_in = list(working_graph.incident_parent_edges(first_vertex, keys=True, data=True))
        directed_out = list(working_graph.incident_child_edges(first_vertex, keys=True, data=True))
        
        if len(directed_in) != 1 or len(directed_out) != 1:
            # Invalid configuration, skip
            continue
        
        # Get edge data
        (u, _, k1, d1) = directed_in[0]
        (_, v, k2, d2) = directed_out[0]
        
        # Merge attributes using the helper function
        merged_attrs = _merge_attrs_for_degree2_suppression_directed(d1 or {}, d2 or {})
        
        # Suppress the degree-2 node with merged attributes
        dm_suppress_degree2_node(working_graph, first_vertex, merged_attrs=merged_attrs)
    
    # Create and return new network from the modified graph (will be validated)
    return dnetwork_from_dmgraph(working_graph)

