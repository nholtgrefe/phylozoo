"""
Internal utilities for directed network operations.

This module provides internal helper functions for network transformations,
including edge attribute merging and degree-2 node suppression.
"""

from typing import TYPE_CHECKING, Any

from ...primitives.d_multigraph.transformations import suppress_degree2_node as dm_suppress_degree2_node

if TYPE_CHECKING:
    from ...primitives.d_multigraph.base import DirectedMultiGraph


def _merge_attrs_for_degree2_suppression_directed(
    edge1_data: dict[str, Any],
    edge2_data: dict[str, Any],
) -> dict[str, Any]:
    """
    Merge edge attributes for suppressing a degree-2 node in a directed network.
    
    This function handles the directed-only case where both edges are directed.
    It uses simpler semantics than the mixed-graph version: only branch_length
    and gamma (from edge2) are kept; all other attributes are removed.
    
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
    This is an internal helper function for suppressing degree-2 nodes in directed networks.
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
    
    # Handle gamma: use edge2's gamma if present (directed networks use simpler semantics)
    gamma2 = edge2_data.get('gamma') if edge2_data else None
    if gamma2 is not None:
        merged['gamma'] = gamma2
    
    # All other attributes (bootstrap, etc.) are removed
    return merged


def _suppress_deg2_nodes(
    graph: 'DirectedMultiGraph',
    exclude_nodes: set[Any] | None = None,
) -> None:
    """
    Suppress all degree-2 nodes in a directed multigraph.
    
    This function is intended for DirectedPhyNetwork operations. It iteratively suppresses
    all nodes with indegree=1 and outdegree=1 (the only valid degree-2 configuration for
    directed graphs). Suppression continues until no more degree-2 nodes remain, as
    suppression may create new degree-2 nodes.
    
    Edge attributes are properly merged using `_merge_attrs_for_degree2_suppression_directed`,
    which handles branch_length (summed) and gamma (from edge2) correctly.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The directed multigraph to modify. **This function modifies the graph in place.**
    exclude_nodes : set[Any] | None, optional
        Set of nodes to exclude from suppression. If None, no nodes are excluded.
        Default is None.
    
    Notes
    -----
    - The function modifies the graph in place.
    - Suppression may create parallel edges.
    - Nodes in `exclude_nodes` are never suppressed, even if they are degree-2.
    - The function handles the iterative nature of suppression (new degree-2 nodes may
      be created as others are suppressed).
    """
    if exclude_nodes is None:
        exclude_nodes = set()
    
    # Iteratively suppress degree-2 nodes until no more remain
    while True:
        # Find all degree-2 nodes (indegree=1 and outdegree=1)
        degree2_nodes = [
            node for node in graph.nodes()
            if node not in exclude_nodes
            and graph._graph.in_degree(node) == 1
            and graph._graph.out_degree(node) == 1
        ]
        
        if not degree2_nodes:
            break
        
        # Process each degree-2 node
        for node in degree2_nodes:
            # Defensive check: node may have been removed by previous suppression
            if node not in graph.nodes():
                continue
            
            # Verify degree is still 2
            if graph._graph.in_degree(node) != 1 or graph._graph.out_degree(node) != 1:
                continue
            
            # Collect incident edges
            incoming = list(graph.incident_parent_edges(node, keys=True, data=True))
            outgoing = list(graph.incident_child_edges(node, keys=True, data=True))
            
            # Should have exactly one incoming and one outgoing for degree-2 node
            if len(incoming) != 1 or len(outgoing) != 1:
                continue
            
            # Get edge data
            (u, _, k1, d1) = incoming[0]
            (_, v, k2, d2) = outgoing[0]
            
            # Merge attributes using helper function
            merged_attrs = _merge_attrs_for_degree2_suppression_directed(
                dict(d1) if d1 else {},
                dict(d2) if d2 else {}
            )

            # Suppress degree-2 node
            dm_suppress_degree2_node(graph, node, merged_attrs=merged_attrs)


def _merge_attrs_for_parallel_identification_directed(
    edges_data: list[dict[str, Any]],
) -> dict[str, Any]:
    """
    Merge edge attributes for identifying parallel edges in a directed network.
    
    Special handling:
    - branch_length: Weighted average by gamma values if gammas are present,
      otherwise simple average. If no branch_lengths are present, it is omitted.
    - gamma: Sum all gamma values from parallel edges.
    - All other attributes (bootstrap, etc.) are removed.
    
    Parameters
    ----------
    edges_data : list[dict[str, Any]]
        List of edge data dictionaries for parallel edges between the same pair of nodes.
    
    Returns
    -------
    dict[str, Any]
        Merged attributes containing branch_length (weighted average if present) and
        gamma (sum of all gammas). All other attributes are removed.
    
    Notes
    -----
    This is an internal helper function for identify_parallel_edges.
    Branch length is computed as a weighted average: sum(branch_length_i * gamma_i) / sum(gamma_i)
    when gammas are present, or as a simple average when no gammas are present.
    """
    merged: dict[str, Any] = {}
    
    if not edges_data:
        return merged
    
    # Collect branch_lengths and gammas, paired by edge index
    edge_pairs: list[tuple[float | None, float | None]] = []
    has_gamma = False
    
    for edge_data in edges_data:
        if edge_data:
            bl = edge_data.get('branch_length')
            gamma = edge_data.get('gamma')
            edge_pairs.append((bl, gamma))
            if gamma is not None:
                has_gamma = True
    
    # Compute branch_length: weighted average by gamma if gammas exist, otherwise simple average
    branch_lengths = [bl for bl, _ in edge_pairs if bl is not None]
    if branch_lengths:
        if has_gamma:
            # Weighted average: sum(branch_length_i * gamma_i) / sum(gamma_i)
            # Only use edges that have both branch_length and gamma
            weighted_sum = 0.0
            gamma_sum = 0.0
            for bl, gamma in edge_pairs:
                if bl is not None and gamma is not None:
                    weighted_sum += bl * gamma
                    gamma_sum += gamma
            if gamma_sum > 0:
                merged['branch_length'] = weighted_sum / gamma_sum
            else:
                # Fallback to simple average if no valid pairs
                merged['branch_length'] = sum(branch_lengths) / len(branch_lengths)
        else:
            # Simple average when no gammas
            merged['branch_length'] = sum(branch_lengths) / len(branch_lengths)
    
    # Sum gamma values from all parallel edges
    if has_gamma:
        gammas = [gamma for _, gamma in edge_pairs if gamma is not None]
        merged['gamma'] = sum(gammas)
    
    # All other attributes (bootstrap, etc.) are removed
    return merged
