"""
Network transformations module.

This module provides functions to transform directed phylogenetic networks
(e.g., suppress degree-2 nodes, convert to LSA network, etc.).
"""

from typing import TYPE_CHECKING, Any

from . import DirectedPhyNetwork
from .features import k_blobs
from .conversions import dnetwork_from_graph
from ._utils import (
    _merge_attrs_for_degree2_suppression_directed,
    _merge_attrs_for_parallel_identification_directed,
    _suppress_deg2_nodes,
)
from ...primitives.d_multigraph.transformations import (
    identify_parallel_edge as dm_identify_parallel_edge,
    identify_vertices as dm_identify_vertices,
)
from ....utils.exceptions import PhyloZooValueError, PhyloZooNetworkError, PhyloZooAlgorithmError

if TYPE_CHECKING:
    from ...primitives.d_multigraph import DirectedMultiGraph

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
    PhyloZooValueError
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
        raise PhyloZooValueError("Cannot create LSA network from empty network")
    
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
    
    # Create new network, preserving graph attributes
    return DirectedPhyNetwork(
        edges=new_edges,
        nodes=new_nodes if new_nodes else None,
        attributes=network.get_network_attribute() if network.get_network_attribute() else None
    )


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
    
    Raises
    ------
    PhyloZooAlgorithmError
        If the algorithm exceeds the maximum number of iterations.
    PhyloZooValueError
        If the network is empty or has no leaves.
    PhyloZooNetworkError
        If the network is invalid.
        
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
        
        # Step 2: Suppress all degree-2 nodes (excluding leaves)
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
    return dnetwork_from_graph(working_graph)


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
    PhyloZooNetworkError
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
        
        # Get edge data for merging
        (u, _, k1, d1) = directed_in[0]
        (_, v, k2, d2) = directed_out[0]
        
        # Merge attributes using the helper function
        merged_attrs = _merge_attrs_for_degree2_suppression_directed(d1 or {}, d2 or {})
        
        # Suppress the degree-2 node with merged attributes
        # Note: We suppress a single specific node here, not all degree-2 nodes
        from ...primitives.d_multigraph.transformations import suppress_degree2_node as dm_suppress_degree2_node
        dm_suppress_degree2_node(working_graph, first_vertex, merged_attrs=merged_attrs)
    
    # Create and return new network from the modified graph (will be validated)
    return dnetwork_from_graph(working_graph)


def _compute_caterpillar_gammas(
    original_gammas: list[float],
) -> list[tuple[float, float]]:
    """
    Compute gamma values for hybrid edges in a caterpillar structure.
    
    For a hybrid node with in-degree k > 2, we create a caterpillar where each
    caterpillar node becomes a hybrid node with 2 incoming edges. The gammas
    are computed to preserve ratios as we go down the caterpillar.
    
    Parameters
    ----------
    original_gammas : list[float]
        List of original gamma values for edges entering the hybrid node,
        in order. All edges go into the caterpillar (the node is replaced).
        Must sum to 1.0.
    
    Returns
    -------
    list[tuple[float, float]]
        List of (gamma1, gamma2) pairs for each caterpillar hybrid node.
        Each pair sums to 1.0. The first pair preserves the ratio of the
        original top two edges. Subsequent pairs use (1.0 - next, next)
        where gamma1 is for the tree edge from previous cat node and gamma2
        is for the hybrid edge from the next parent.
    
    Examples
    --------
    >>> # Original: [0.1, 0.2, 0.3, 0.4] (sums to 1.0)
    >>> # First pair preserves 0.1:0.2 = 1:2 ratio = (0.3333, 0.6666)
    >>> gammas = _compute_caterpillar_gammas([0.1, 0.2, 0.3, 0.4])
    >>> len(gammas) == 3
    True
    >>> # First pair: (0.1/(0.1+0.2), 0.2/(0.1+0.2)) = (0.3333, 0.6666)
    >>> abs(gammas[0][0] - 0.3333) < 0.01 and abs(gammas[0][1] - 0.6666) < 0.01
    True
    >>> # Second pair: (1.0 - 0.3, 0.3) = (0.7, 0.3)
    >>> abs(gammas[1][0] - 0.7) < 0.01 and abs(gammas[1][1] - 0.3) < 0.01
    True
    >>> # Third pair: (1.0 - 0.4, 0.4) = (0.6, 0.4)
    >>> abs(gammas[2][0] - 0.6) < 0.01 and abs(gammas[2][1] - 0.4) < 0.01
    True
    """
    if len(original_gammas) < 2:
        return []
    
    g1 = original_gammas[0]  # First edge is kept with its original gamma
    remaining = original_gammas[1:]  # Rest go into caterpillar
    
    if len(remaining) == 1:
        # Only one remaining edge: (1.0 - next, next)
        g2 = remaining[0]
        return [(1.0 - g2, g2)]
    
    # First caterpillar node: preserves ratio of top two edges
    g2 = remaining[0]
    top_two_total = g1 + g2
    if top_two_total == 0:
        first_pair = (0.5, 0.5)
    else:
        first_pair = (g1 / top_two_total, g2 / top_two_total)
    
    caterpillar_pairs = [first_pair]
    
    # Process remaining edges: each pair is (1.0 - next, next)
    for i in range(1, len(remaining)):
        g_next = remaining[i]
        pair = (1.0 - g_next, g_next)
        caterpillar_pairs.append(pair)
    
    return caterpillar_pairs


def _binary_resolve_tree_nodes(
    graph: 'DirectedMultiGraph',  # type: ignore[name-defined]
    has_branch_lengths: bool,
) -> None:
    """
    Resolve high out-degree nodes in a directed multigraph using caterpillar structures.
    
    This function modifies the graph in place by replacing nodes with out-degree > 2
    with caterpillar structures. The node is completely removed and replaced by a
    caterpillar chain. All incoming edges to the original node are redirected to the
    first caterpillar node, and outgoing edges are distributed through the caterpillar.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The directed multigraph to modify in place.
    has_branch_lengths : bool
        Whether the graph has branch lengths. If True, new edges get branch_length=0.0.
    
    Notes
    -----
    - The function modifies the graph in place.
    - Only branch_length attributes are preserved on edges.
    - New caterpillar edges get branch_length=0.0 if has_branch_lengths is True.
    - The high out-degree node is completely removed and replaced by caterpillar nodes.
    - Structure: For node u with children [v1, v2, v3, v4], creates:
      cat1 -> v1, cat1 -> cat2, cat2 -> v2, cat2 -> cat3, cat3 -> v3, cat3 -> v4.
      All incoming edges to u are redirected to cat1.
    """
    # Helper function to filter edge attributes (keep only branch_length)
    def _filter_attrs(data: dict) -> dict:
        """Filter edge attributes to keep only branch_length."""
        filtered = {}
        if 'branch_length' in data:
            filtered['branch_length'] = data['branch_length']
        return filtered
    
    # Find high out-degree nodes
    high_outdegree_nodes = [
        n for n in graph.nodes() 
        if graph.outdegree(n) > 2
    ]
    
    # Process each high out-degree node
    for u in high_outdegree_nodes:
        out_neighbors = list(graph.successors(u))
        if len(out_neighbors) <= 2:
            continue  # May have been resolved already
        
        # Collect all outgoing edges and their data
        out_edges_data = {}
        for v in out_neighbors:
            edge_data = next(iter(graph._graph[u][v].values()))
            out_edges_data[v] = _filter_attrs(edge_data)
            graph.remove_edge(u, v)
        
        # Collect all incoming edges (to reconnect to first caterpillar node)
        in_predecessors = list(graph.predecessors(u))
        in_edges_data = {}
        for pred in in_predecessors:
            edge_data = next(iter(graph._graph[pred][u].values()))
            in_edges_data[pred] = _filter_attrs(edge_data)
            graph.remove_edge(pred, u)
        
        # Remove the node
        graph.remove_node(u)
        
        # Build caterpillar: cat1 -> v1, cat1 -> cat2, cat2 -> v2, cat2 -> cat3, cat3 -> v3, cat3 -> v4, ...
        # Structure: first cat node connects to first child, then chain of cats connecting to remaining children
        first_v = out_neighbors[0]
        rest_v = out_neighbors[1:]
        
        # Create first caterpillar node
        cat1 = next(graph.generate_node_ids(1))
        graph.add_node(cat1)
        
        # Connect first child to first cat node
        graph.add_edge(cat1, first_v, **out_edges_data[first_v])
        
        # Build caterpillar chain for remaining children
        last_cat = cat1
        for i, v in enumerate(rest_v):
            if i < len(rest_v) - 1:
                # Create next caterpillar node
                next_cat = next(graph.generate_node_ids(1))
                graph.add_node(next_cat)
                # Connect: last_cat -> next_cat
                tree_attrs = {'branch_length': 0.0} if has_branch_lengths else {}
                graph.add_edge(last_cat, next_cat, **tree_attrs)
                # Connect: next_cat -> v
                graph.add_edge(next_cat, v, **out_edges_data[v])
                last_cat = next_cat
            else:
                # Last child: connect directly to last caterpillar node
                graph.add_edge(last_cat, v, **out_edges_data[v])
        
        # Reconnect all incoming edges to first caterpillar node
        for pred in in_predecessors:
            graph.add_edge(pred, cat1, **in_edges_data[pred])


def _binary_resolve_hybrid_nodes(
    graph: 'DirectedMultiGraph',  # type: ignore[name-defined]
    has_branch_lengths: bool,
) -> None:
    """
    Resolve high in-degree nodes (hybrid nodes) in a directed multigraph using caterpillar structures.
    
    This function modifies the graph in place by replacing nodes with in-degree > 2
    with caterpillar structures. The node is completely removed and replaced by a
    caterpillar chain, with all outgoing edges redirected to the last caterpillar node.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The directed multigraph to modify in place.
    has_branch_lengths : bool
        Whether the graph has branch lengths. If True, new edges get branch_length=0.0.
    
    Notes
    -----
    - The function modifies the graph in place.
    - Branch_length and gamma attributes are preserved on edges.
    - Gamma values are computed using _compute_caterpillar_gammas to preserve ratios.
      The top two hybrid edges (u1 -> cat1, u2 -> cat1) preserve the ratio of the
      original top two edges. Subsequent caterpillar nodes use (1.0 - gamma_next, gamma_next)
      for the tree edge and hybrid edge respectively.
    - New caterpillar edges get branch_length=0.0 if has_branch_lengths is True.
    - The high in-degree node is completely removed and replaced by caterpillar nodes.
    - Structure: For node v with parents [u1, u2, u3, u4], creates:
      u1 -> cat1, u2 -> cat1, cat1 -> cat2, u3 -> cat2, cat2 -> cat3, u4 -> cat3.
      All outgoing edges from v are redirected to cat3 (the last caterpillar node).
    """
    # Helper function to filter edge attributes (keep only branch_length and gamma)
    def _filter_attrs(data: dict) -> dict:
        """Filter edge attributes to keep only branch_length and gamma."""
        filtered = {}
        if 'branch_length' in data:
            filtered['branch_length'] = data['branch_length']
        if 'gamma' in data:
            filtered['gamma'] = data['gamma']
        return filtered
    
    # Find high in-degree nodes
    high_indegree_nodes = [
        n for n in graph.nodes() 
        if graph.indegree(n) > 2
    ]
    
    # Process each high in-degree node
    for v in high_indegree_nodes:
        in_neighbors = list(graph.predecessors(v))
        if len(in_neighbors) <= 2:
            continue  # May have been resolved already
        
        # Collect all incoming edges and their data (with gammas)
        in_edges_data = {}
        original_gammas = []
        for u in in_neighbors:
            edge_data = next(iter(graph._graph[u][v].values()))
            original_gammas.append(edge_data.get('gamma', 0.0))
            in_edges_data[u] = _filter_attrs(edge_data)
            graph.remove_edge(u, v)
        
        # Collect all outgoing edges (to reconnect to last caterpillar node)
        out_successors = list(graph.successors(v))
        out_edges_data = {}
        for succ in out_successors:
            edge_data = next(iter(graph._graph[v][succ].values()))
            out_edges_data[succ] = _filter_attrs(edge_data)
            graph.remove_edge(v, succ)
        
        # Remove the node
        graph.remove_node(v)
        
        # Compute caterpillar gammas
        caterpillar_pairs = _compute_caterpillar_gammas(original_gammas)
        
        # Build caterpillar: u1 -> cat1, u2 -> cat1, cat1 -> cat2, u3 -> cat2, cat2 -> cat3, u4 -> cat3
        first_u = in_neighbors[0]
        rest_u = in_neighbors[1:]
        
        if len(rest_u) == 1:
            # Special case: only one additional edge
            # Create single caterpillar node: u1 -> cat1, u2 -> cat1
            cat_node = next(graph.generate_node_ids(1))
            graph.add_node(cat_node)
            
            gamma_cumulative, gamma_next = caterpillar_pairs[0]
            
            # Edge from u1 to cat_node (with gamma)
            cum_attrs = in_edges_data[first_u].copy()
            cum_attrs['gamma'] = gamma_cumulative
            graph.add_edge(first_u, cat_node, **cum_attrs)
            
            # Hybrid edge: u2 -> cat_node (with gamma)
            hybrid_attrs = in_edges_data[rest_u[0]].copy()
            hybrid_attrs['gamma'] = gamma_next
            graph.add_edge(rest_u[0], cat_node, **hybrid_attrs)
            
            last_cat = cat_node
        else:
            # Multiple additional edges: build caterpillar chain
            # Structure: u1 -> cat1, u2 -> cat1, cat1 -> cat2, u3 -> cat2, cat2 -> cat3, u4 -> cat3
            last_cat = None
            for i, u in enumerate(rest_u):
                # Create caterpillar node for this neighbor
                cat_node = next(graph.generate_node_ids(1))
                graph.add_node(cat_node)
                
                # Get gamma pair for this caterpillar node
                gamma_cumulative, gamma_next = caterpillar_pairs[i]
                
                if i == 0:
                    # First caterpillar node: receives from u1 and u2
                    # Top two hybrid edges preserve the ratio of original top two edges
                    # Edge from u1 to cat_node (with gamma preserving ratio)
                    cum_attrs = in_edges_data[first_u].copy()
                    cum_attrs['gamma'] = gamma_cumulative
                    graph.add_edge(first_u, cat_node, **cum_attrs)
                    
                    # Hybrid edge: u2 -> cat_node (with gamma preserving ratio)
                    hybrid_attrs = in_edges_data[u].copy()
                    hybrid_attrs['gamma'] = gamma_next
                    graph.add_edge(u, cat_node, **hybrid_attrs)
                else:
                    # Subsequent caterpillar nodes: receive from previous cat and next u
                    # Tree edge from previous hybrid node: (1.0 - gamma_next, gamma_next)
                    tree_attrs = {}
                    if has_branch_lengths:
                        tree_attrs['branch_length'] = 0.0
                    tree_attrs['gamma'] = gamma_cumulative  # This is (1.0 - gamma_next)
                    graph.add_edge(last_cat, cat_node, **tree_attrs)
                    
                    # Hybrid edge: u -> cat_node (with gamma_next)
                    hybrid_attrs = in_edges_data[u].copy()
                    hybrid_attrs['gamma'] = gamma_next
                    graph.add_edge(u, cat_node, **hybrid_attrs)
                
                last_cat = cat_node
        
        # Reconnect all outgoing edges from v to the last caterpillar node
        for succ in out_successors:
            graph.add_edge(last_cat, succ, **out_edges_data[succ])


def binary_resolution(network: DirectedPhyNetwork) -> DirectedPhyNetwork:
    """
    Convert a non-binary network to a binary network by resolving high-degree nodes.
    
    This function creates a binary resolution of the network by replacing nodes with
    in-degree > 2 or out-degree > 2 with caterpillar structures. The first neighbor
    (incoming or outgoing) is kept, and additional neighbors are connected through
    a caterpillar chain.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to convert.
    
    Returns
    -------
    DirectedPhyNetwork
        A new binary network. If the input network is already binary, returns a copy.
    
    Raises
    ------
    PhyloZooValueError
        If the network has parallel edges (binary resolution is not defined for
        networks with parallel edges).
    
    Notes
    -----
    - Attribute handling:
      - All attributes are removed except branch_length and gamma
      - Branch length handling:
        - If the original network had branch lengths on any edges, new edges in the
          caterpillar structures are assigned branch_length=0.0
        - Original edges keep their branch_length values
      - Gamma handling (for high in-degree nodes only):
        - The top two hybrid edges in the caterpillar maintain the same ratio as
          the original top two hybrid edges
        - As we go down the caterpillar, each new hybrid edge's gamma accounts for
          the cumulative probability from previous edges
        - Gamma values are computed using _compute_caterpillar_gammas()
    - Node labels are preserved
    - The network must not have parallel edges
    
    Examples
    --------
    >>> # Network with high out-degree node
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         (5, 1), (5, 2), (5, 3), (5, 4)  # Node 5 has out-degree 4
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (4, {'label': 'D'})]
    ... )
    >>> binary_net = binary_resolution(net)
    >>> binary_net.is_binary()
    True
    """
    from .classifications import has_parallel_edges, is_binary
    from .conversions import dnetwork_from_graph
    
    # Check for parallel edges
    if has_parallel_edges(network):
        raise PhyloZooValueError(
            "binary_resolution cannot be applied to networks with parallel edges"
        )
    
    # Check if network is already binary
    if is_binary(network):
        return network.copy()
    
    # Check if network has any branch lengths
    has_branch_lengths = False
    for u, v, key, data in network._graph.edges(keys=True, data=True):
        if data.get('branch_length') is not None:
            has_branch_lengths = True
            break
    
    # Work on a copy of the graph
    working_graph = network._graph.copy()
    
    # Resolve high out-degree nodes (tree nodes)
    _binary_resolve_tree_nodes(working_graph, has_branch_lengths)
    
    # Resolve high in-degree nodes (hybrid nodes)
    _binary_resolve_hybrid_nodes(working_graph, has_branch_lengths)
    
    # Create and return new network from the modified graph
    return dnetwork_from_graph(working_graph)

