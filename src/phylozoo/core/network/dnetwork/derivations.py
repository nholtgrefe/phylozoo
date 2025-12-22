"""
Network derivations module.

This module provides functions to derive other data structures from directed
phylogenetic networks (e.g., splits, quartets, distances, blobtrees, subnetworks, etc.).
"""

# TODO: Implement derivation functions here.
# These could include:
# - Split extraction
# - Quartet extraction
# - Distance calculations
# - Blobtree construction
# - Subnetwork extraction
# - Displayed tree extraction
from typing import Any

from . import DirectedPhyNetwork
from .classifications import is_lsa_network
from .features import blobs
from .io import dnetwork_from_dmgraph
from .transformations import (
    to_lsa_network,
    suppress_2_blobs,
    _merge_edge_attributes_for_suppression,
)
from ...primitives.d_multigraph.transformations import identify_vertices as dm_identify_vertices
from ...primitives.m_multigraph import MixedMultiGraph
from ...primitives.m_multigraph.transformations import suppress_degree2_node
from ..sdnetwork import SemiDirectedPhyNetwork
from ..sdnetwork.io import sdnetwork_from_mmgraph


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

    # Ensure node labels are preserved in the mixed graph
    for node in mixed.nodes():
        if node in working._node_to_label:
            label = working._node_to_label[node]
            mixed._undirected.nodes[node]['label'] = label
            mixed._directed.nodes[node]['label'] = label
            mixed._combined.nodes[node]['label'] = label

    # Convert the mixed graph to a semi-directed network
    return sdnetwork_from_mmgraph(mixed, network_type='semi-directed')


def tree_of_blobs(network: DirectedPhyNetwork) -> DirectedPhyNetwork:
    """
    Create a tree of blobs by suppressing all 2-blobs and collapsing internal blobs.

    This function:
    1. Suppresses all 2-blobs using suppress_2_blobs
    2. Finds all internal blobs (blobs with more than 1 node, excluding leaves)
    3. For each internal blob, identifies all vertices with a single vertex
    4. Returns a new network representing the tree of blobs

    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to transform.

    Returns
    -------
    DirectedPhyNetwork
        A new network where each blob has been collapsed to a single vertex,
        forming a tree structure.

    Examples
    --------
    >>> # Create a directed network with a hybrid
    >>> from phylozoo.core.network.dnetwork.classifications import is_tree
    >>> dnet = DirectedPhyNetwork(
    ...     edges=[
    ...         (10, 5), (10, 6),  # Root to tree nodes
    ...         (5, 4), (6, 4),    # Both lead to hybrid 4 (in-degree 2)
    ...         (4, 8),            # Hybrid to tree node
    ...         (8, 1), (8, 2),    # Tree node to leaves
    ...         (5, 3), (6, 7)     # Additional leaves to satisfy degree constraints
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
    ... )
    >>> tree_net = tree_of_blobs(dnet)
    >>> is_tree(tree_net)  # Should be True
    True
    """
    # First suppress all 2-blobs using the existing function
    blob_network = suppress_2_blobs(network)

    # Find all internal blobs (more than 1 node, not containing only leaves)
    all_blobs = blobs(blob_network, trivial=False, leaves=False)

    # Work on a copy of the internal graph and collapse blobs
    working_graph = blob_network._graph.copy()
    for blob in all_blobs:
        if len(blob) > 1:
            # Sort to ensure deterministic behavior
            blob_sorted = sorted(blob)
            # Identify all vertices in the blob with the first vertex
            dm_identify_vertices(working_graph, blob_sorted)

    # Convert back to DirectedPhyNetwork
    return dnetwork_from_dmgraph(working_graph)
