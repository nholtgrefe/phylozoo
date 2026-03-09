"""
Network derivations module.

This module provides functions to derive other data structures from directed
phylogenetic networks (e.g., splits, quartets, distances, blobtrees, subnetworks, etc.).
"""

import itertools
from typing import Any, Iterator, Literal

import numpy as np
import networkx as nx

from . import DirectedPhyNetwork
from .classifications import is_lsa_network
from .features import blobs
from .conversions import dnetwork_from_graph
from .transformations import (
    to_lsa_network,
    suppress_2_blobs as suppress_2_blobs_fn,
    identify_parallel_edges as identify_parallel_edges_dn,
)
from ...split import Split, SplitSystem, WeightedSplitSystem
from ...quartet import QuartetProfileSet
from ...primitives.partition import Partition
from ._utils import _suppress_deg2_nodes as dm_suppress_deg2_nodes
from ..sdnetwork._utils import _suppress_deg2_nodes as mm_suppress_deg2_nodes
from ...primitives.d_multigraph.transformations import identify_vertices as dm_identify_vertices
from ...primitives.d_multigraph.transformations import subgraph as dm_subgraph
from ...primitives.d_multigraph import DirectedMultiGraph
from ...primitives.m_multigraph import MixedMultiGraph
from ..sdnetwork import SemiDirectedPhyNetwork
from ..sdnetwork.conversions import sdnetwork_from_graph
from ....core.distance import DistanceMatrix
from ....utils.exceptions import PhyloZooValueError, PhyloZooAlgorithmError

def to_sd_network(d_network: DirectedPhyNetwork) -> SemiDirectedPhyNetwork:
    """
    Convert a DirectedPhyNetwork to a SemiDirectedPhyNetwork.

    Steps:

    1. If the directed network is not an LSA network, replace it by its LSA-network.
    2. Undirect all non-hybrid edges; hybrid edges remain directed.
    3. Suppress any degree-2 node (this stems from a degree-2 root). Suppression may
       create parallel edges. Suppression connects the two neighbors directly:
       undirected+undirected -> undirected; directed+directed (u->x, x->v) -> directed
       (u->v); directed into x and undirected out (u->x, x-v) -> undirected (u-v);
       undirected into x and directed out (u-x, x->v) -> directed (u->v).

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

    # Suppress all degree-2 nodes using the mixed graph utility function
    mm_suppress_deg2_nodes(mixed)

    # Ensure node labels are preserved in the mixed graph
    for node in mixed.nodes():
        if node in working._node_to_label:
            label = working._node_to_label[node]
            mixed._undirected.nodes[node]['label'] = label
            mixed._directed.nodes[node]['label'] = label
            mixed._combined.nodes[node]['label'] = label

    # Convert the mixed graph to a semi-directed network
    return sdnetwork_from_graph(mixed, network_type='semi-directed')


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
    blob_network = suppress_2_blobs_fn(network)

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
    return dnetwork_from_graph(working_graph)


def subnetwork(
    network: DirectedPhyNetwork,
    taxa: list[str],
    suppress_2_blobs: bool = False,
    identify_parallel_edges: bool = False,
    make_lsa: bool = False,
) -> DirectedPhyNetwork:
    """
    Extract the subnetwork induced by a subset of taxa (leaf labels).

    The subnetwork is defined as the union of all directed paths from the
    requested leaves up to the root (i.e., all their ancestors and the
    leaves themselves). The induced subgraph is taken on the underlying
    DirectedMultiGraph, then degree-2 internal nodes are suppressed. Optionally,
    the result can be post-processed by suppressing 2-blobs, identifying
    parallel edges, and/or converting to an LSA network. After any of these
    optional steps, degree-2 suppression is applied again to clean up artifacts.

    Parameters
    ----------
    network : DirectedPhyNetwork
        Source network.
    taxa : list[str]
        Subset of taxon labels (leaf labels) to induce the subnetwork on.
    suppress_2_blobs : bool, default False
        If True, suppress all 2-blobs in the resulting network.
    identify_parallel_edges : bool, default False
        If True, identify/merge parallel edges in the resulting network.
    make_lsa : bool, default False
        If True, convert the result to an LSA-network.

    Returns
    -------
    DirectedPhyNetwork
        The derived subnetwork. Returns an empty network if `taxa` is empty.

    Raises
    ------
    PhyloZooValueError
        If any of the provided taxa are not found in the network.
    """
    if not taxa:
        # Return an empty network if no taxa are specified
        return DirectedPhyNetwork(edges=[], nodes=[])

    # Work on the original network; LSA conversion (if requested) will be
    # applied after optional post-processing (user requested behavior).
    working_net = network

    # Map taxa labels to node ids using public API
    leaf_nodes: list[Any] = []
    for t in taxa:
        node_id = working_net.get_node_id(t)
        if node_id is None:
            raise PhyloZooValueError(f"Taxon label '{t}' not found in network")
        leaf_nodes.append(node_id)

    # Collect all ancestors (and the leaves themselves)
    import networkx as nx

    dag = working_net._graph._graph
    nodes_set: set[Any] = set()
    for leaf in leaf_nodes:
        nodes_set.add(leaf)
        nodes_set.update(nx.ancestors(dag, leaf))

    # Create induced DirectedMultiGraph using existing utility
    induced_dm = dm_subgraph(working_net._graph, nodes_set)

    # Work on a mutable copy for transformations
    working_dm = induced_dm.copy()

    # First pass: suppress all degree-2 nodes (directed suppression semantics)
    dm_suppress_deg2_nodes(working_dm, exclude_nodes=None)

    # Convert to DirectedPhyNetwork for higher-level transformations
    result_net = dnetwork_from_graph(working_dm)

    # Optional post-processing steps
    if suppress_2_blobs:
        result_net = suppress_2_blobs_fn(result_net)

    if identify_parallel_edges:
        result_net = identify_parallel_edges_dn(result_net)

    if make_lsa:
        result_net = to_lsa_network(result_net)

    return result_net


def k_taxon_subnetworks(
    network: DirectedPhyNetwork,
    k: int,
    suppress_2_blobs: bool = False,
    identify_parallel_edges: bool = False,
    make_lsa: bool = False,
) -> Iterator[DirectedPhyNetwork]:
    """
    Generate all subnetworks induced by exactly k taxa.

    This function yields all possible subnetworks of the network that are
    induced by exactly k taxon labels. For each combination of k taxa,
    the corresponding subnetwork is computed using the `subnetwork` function.

    Parameters
    ----------
    network : DirectedPhyNetwork
        Source network.
    k : int
        Number of taxa to include in each subnetwork. Must be between 0 and
        the number of taxa in the network (inclusive).
    suppress_2_blobs : bool, default False
        If True, suppress all 2-blobs in each resulting subnetwork.
    identify_parallel_edges : bool, default False
        If True, identify/merge parallel edges in each resulting subnetwork.
    make_lsa : bool, default False
        If True, convert each result to an LSA-network.

    Yields
    ------
    DirectedPhyNetwork
        Subnetworks induced by exactly k taxa. Each subnetwork is generated
        lazily as the iterator is consumed.

    Raises
    ------
    PhyloZooValueError
        If k < 0 or k > number of taxa in the network.

    Examples
    --------
    >>> net = DirectedPhyNetwork(
    ...     edges=[(5, 3), (5, 4), (3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
    ... )
    >>> # Generate all 2-taxon subnetworks
    >>> subnetworks = list(k_taxon_subnetworks(net, k=2))
    >>> len(subnetworks)
    3  # C(3,2) = 3 combinations
    >>> # Each subnetwork has exactly 2 leaves
    >>> all(len(subnet.leaves) == 2 for subnet in subnetworks)
    True
    >>> # Generate all 1-taxon subnetworks
    >>> single_taxon_subs = list(k_taxon_subnetworks(net, k=1))
    >>> len(single_taxon_subs)
    3  # C(3,1) = 3 combinations
    """
    all_taxa = list(network.taxa)
    num_taxa = len(all_taxa)

    # Validate k
    if k < 0:
        raise PhyloZooValueError(f"k must be non-negative, got {k}")
    if k > num_taxa:
        raise PhyloZooValueError(
            f"k ({k}) cannot exceed the number of taxa ({num_taxa}) in the network"
        )

    # Generate all combinations of k taxa
    for taxa_combination in itertools.combinations(all_taxa, k):
        yield subnetwork(
            network,
            list(taxa_combination),
            suppress_2_blobs=suppress_2_blobs,
            identify_parallel_edges=identify_parallel_edges,
            make_lsa=make_lsa,
        )


def _switchings(network: DirectedPhyNetwork, probability: bool = False) -> Iterator[DirectedMultiGraph]:
    """
    Generate all switchings of a directed phylogenetic network.
    
    A switching is obtained by deleting all but one incident parent edge for each
    hybrid node. This function generates all possible combinations of keeping
    exactly one parent edge per hybrid node. Each switching is a tree (not
    necessarily a phylogenetic tree).
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network.
    probability : bool, optional
        If True, store the probability of the switching in the graph's 'probability'
        attribute. The probability is the product of gamma values for the kept hybrid
        edges. If a hybrid edge has no gamma value, it is taken to be 1/k where k
        is the in-degree of the hybrid node. By default False.
    
    Yields
    ------
    DirectedMultiGraph
        A switching of the network (one parent edge kept per hybrid node). Each
        switching is a tree. If probability=True, the graph has a 'probability'
        attribute containing the switching probability.
    
    Examples
    --------
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         (10, 5), (10, 6),  # Root to tree nodes
    ...         (5, 4), (6, 4),    # Both lead to hybrid 4 (in-degree 2)
    ...         (4, 8),            # Hybrid to tree node
    ...         (8, 1), (8, 2),    # Tree node to leaves
    ...         (5, 3), (6, 7)     # Additional leaves to satisfy degree constraints
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
    ... )
    >>> switchings = list(_switchings(net))
    >>> len(switchings)
    2  # Two parent edges for hybrid node 4: (5,4) and (6,4)
    >>> # Each switching has exactly one parent edge for the hybrid node
    >>> hybrid = 4
    >>> for sw in switchings:
    ...     parent_edges = list(sw.incident_parent_edges(hybrid, keys=True))
    ...     assert len(parent_edges) == 1
    >>> # With probability=True, each switching has a probability attribute
    >>> net_with_gamma = DirectedPhyNetwork(
    ...     edges=[
    ...         (10, 5), (10, 6),
    ...         {'u': 5, 'v': 4, 'gamma': 0.6},
    ...         {'u': 6, 'v': 4, 'gamma': 0.4},
    ...         (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
    ... )
    >>> switchings_with_prob = list(_switchings(net_with_gamma, probability=True))
    >>> switchings_with_prob[0]._graph.graph.get('probability')
    0.6  # Probability of keeping edge (5,4)
    >>> switchings_with_prob[1]._graph.graph.get('probability')
    0.4  # Probability of keeping edge (6,4)
    """
    hybrid_nodes = network.hybrid_nodes
    
    # If no hybrid nodes, return the original graph
    if not hybrid_nodes:
        switching_graph = network._graph.copy()
        if probability:
            switching_graph.set_graph_attribute('probability', 1.0)
        yield switching_graph
        return
    
    # Collect parent edges for each hybrid node (hybrid nodes always have parent edges)
    hybrid_parent_edges: dict[Any, list[tuple[Any, Any, int]]] = {}
    hybrid_indegrees: dict[Any, int] = {}
    for hybrid in hybrid_nodes:
        parent_edges = list(network.incident_parent_edges(hybrid, keys=True))
        hybrid_parent_edges[hybrid] = parent_edges
        hybrid_indegrees[hybrid] = len(parent_edges)
    
    # Generate all combinations: for each hybrid, choose one parent edge to keep
    hybrid_list = list(hybrid_nodes)
    
    for edge_combination in itertools.product(*hybrid_parent_edges.values()):
        # Create a copy of the graph
        switching_graph = network._graph.copy()
        
        # Calculate probability and remove edges in a single loop
        prob = 1.0 if probability else None
        for hybrid, (keep_u, keep_v, keep_key) in zip(hybrid_list, edge_combination):
            # Calculate probability contribution for this hybrid if requested
            if probability:
                gamma = network.get_gamma(keep_u, keep_v, keep_key)
                if gamma is not None:
                    prob *= gamma
                else:
                    # No gamma specified, use uniform probability 1/k
                    prob *= 1.0 / hybrid_indegrees[hybrid]
            
            # Remove all parent edges except the one to keep
            for u, v, key in hybrid_parent_edges[hybrid]:
                if (u, v, key) != (keep_u, keep_v, keep_key):
                    switching_graph.remove_edge(u, v, key=key)
        
        # Store probability if requested (in all underlying graphs)
        if probability:
            switching_graph.set_graph_attribute('probability', prob)
        
        yield switching_graph


def displayed_trees(network: DirectedPhyNetwork, probability: bool = False) -> Iterator[DirectedPhyNetwork]:
    """
    Generate all displayed trees of a directed phylogenetic network.
    
    A displayed tree is obtained by:

    1. Taking a switching (deleting all but one parent edge per hybrid node)
    2. Exhaustively removing degree-1 nodes that are not leaves or the root
    3. Suppressing all degree-2 nodes
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network.
    probability : bool, optional
        If True, store the probability of the displayed tree in the network's
        'probability' attribute. The probability is inherited from the switching
        and equals the product of gamma values for the kept hybrid edges. If a
        hybrid edge has no gamma value, it is taken to be 1/k where k is the
        in-degree of the hybrid node. If there are no hybrid nodes, the
        probability is 1.0. By default False.
    
    Yields
    ------
    DirectedPhyNetwork
        A displayed tree of the network. If probability=True, the network has a
        'probability' attribute containing the tree probability.
    
    Examples
    --------
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         (10, 5), (10, 6),
    ...         (5, 4), (6, 4),
    ...         (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
    ... )
    >>> trees = list(displayed_trees(net))
    >>> len(trees)
    2  # Two switchings yield two displayed trees
    """
    # Get original network's leaves and root for reference
    original_leaves = network.leaves
    original_root = network.root_node
    
    # Iterate through all switchings
    for switching_graph in _switchings(network, probability=probability):
        # Work on a copy to avoid modifying the switching
        tree_graph = switching_graph.copy()
        
        # Exhaustively remove degree-1 nodes that are not leaves or root
        while True:
            degree1_nodes = [
                node for node in tree_graph.nodes()
                if tree_graph.degree(node) == 1
                and node not in original_leaves
                and node != original_root
            ]
            
            if not degree1_nodes:
                break
            
            # Remove all degree-1 nodes (excluding leaves and root)
            for node in degree1_nodes:
                # Double-check node still exists and is still degree-1
                if node in tree_graph.nodes() and tree_graph.degree(node) == 1:
                    tree_graph.remove_node(node)
        
        # Suppress all degree-2 nodes
        dm_suppress_deg2_nodes(tree_graph, exclude_nodes=None)
        
        # Convert back to DirectedPhyNetwork
        # Note: dnetwork_from_graph already copies graph attributes, so probability
        # is automatically preserved from the switching graph.
        displayed_tree = dnetwork_from_graph(tree_graph)
        
        yield displayed_tree


def _switching_distance_matrix(
    switching_graph: DirectedMultiGraph,
    taxa: list[str],
    original_network: DirectedPhyNetwork,
) -> np.ndarray:
    """
    Compute the full distance matrix for all pairwise distances between taxa in a switching.
    
    Parameters
    ----------
    switching_graph : DirectedMultiGraph
        The switching graph (a tree).
    taxa : list[str]
        List of taxon labels.
    original_network : DirectedPhyNetwork
        The original network (used to access branch lengths).
    
    Returns
    -------
    np.ndarray
        A symmetric numpy array of shape (len(taxa), len(taxa)) with pairwise distances.
        Diagonal is 0.0 (distance from a taxon to itself).
    
    Notes
    -----
    The switching is a tree, so there is exactly one path between any two leaves.
    Branch lengths are accessed from the original network. If an edge has no branch
    length, 1.0 is used as the default.
    """
    n = len(taxa)
    if n == 0:
        return np.array([])
    if n == 1:
        return np.array([[0.0]])
    
    # Get leaf nodes corresponding to taxon labels
    leaf_nodes: list[Any] = []
    for taxon in taxa:
        leaf_node = original_network._label_to_node.get(taxon)
        if leaf_node is None:
            raise PhyloZooValueError(f"Taxon '{taxon}' not found in network")
        leaf_nodes.append(leaf_node)
    
    # Initialize distance matrix
    distance_matrix = np.zeros((n, n), dtype=np.float64)
    
    # Use the combined graph view for path finding (treats all edges as undirected)
    combined_graph = switching_graph._combined
    
    # Compute pairwise distances
    for i, leaf1 in enumerate(leaf_nodes):
        for j, leaf2 in enumerate(leaf_nodes):
            if i == j:
                distance_matrix[i, j] = 0.0
                continue
            
            # Find the unique path between the two leaves
            path = nx.shortest_path(combined_graph, leaf1, leaf2)
            
            # Sum branch lengths along the path
            total_distance = 0.0
            for k in range(len(path) - 1):
                u, v = path[k], path[k + 1]
                
                # Get branch length from the original network
                bl = original_network.get_branch_length(u, v)
                if bl is None:
                    bl = original_network.get_branch_length(v, u)
                
                # Default to 1.0 if no branch length found
                if bl is None:
                    bl = 1.0
                
                total_distance += bl
            
            distance_matrix[i, j] = total_distance
            distance_matrix[j, i] = total_distance  # Symmetric
    
    return distance_matrix


def distances(
    network: DirectedPhyNetwork,
    mode: Literal['shortest', 'longest', 'average'] = 'average',
) -> DistanceMatrix:
    """
    Compute pairwise distances between taxa based on switchings.
    
    This function computes distances by considering all switchings of the network.
    For each pair of taxa, the distance is computed in each switching (sum of branch
    lengths along the unique path), and then aggregated according to the specified mode.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network.
    mode : Literal['shortest', 'longest', 'average'], optional
        Distance aggregation mode:

        - 'shortest': Take the minimum distance across all switchings
        - 'longest': Take the maximum distance across all switchings
        - 'average': Take the probability-weighted average across all switchings

        By default 'average'.

    Returns
    -------
    DistanceMatrix
        A distance matrix with pairwise distances between all taxa.
    
    Examples
    --------
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         (10, 5), (10, 6),
    ...         (5, 4), (6, 4),
    ...         (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
    ... )
    >>> dm = distances(net, mode='shortest')
    >>> len(dm)
    4
    >>> dm.get_distance('A', 'B')
    2.0  # Example distance
    """
    # Get all taxa
    all_taxa = list(network.taxa)
    all_taxa.sort()
    n = len(all_taxa)
    
    # Handle edge cases
    if n == 0:
        return DistanceMatrix(np.array([]), labels=[])
    if n == 1:
        return DistanceMatrix(np.array([[0.0]]), labels=all_taxa)
    
    # Initialize aggregation arrays based on mode
    if mode == 'shortest':
        result = np.full((n, n), np.inf, dtype=np.float64)
    elif mode == 'longest':
        result = np.zeros((n, n), dtype=np.float64)
    elif mode == 'average':
        weighted_sum = np.zeros((n, n), dtype=np.float64)
        prob_sum = 0.0
    else:
        raise PhyloZooValueError(f"Invalid mode: {mode}. Must be 'shortest', 'longest', or 'average'")
    
    # Iterate through all switchings
    for switching_graph in _switchings(network, probability=(mode == 'average')):
        # Get probability if mode is 'average'
        prob = 1.0
        if mode == 'average':
            prob = switching_graph._graph.graph.get('probability', 1.0)
            if prob is None:
                prob = 1.0
        
        # Compute full distance matrix for this switching
        switching_matrix = _switching_distance_matrix(switching_graph, all_taxa, network)
        
        # Update aggregation based on mode
        if mode == 'shortest':
            result = np.minimum(result, switching_matrix)
        elif mode == 'longest':
            result = np.maximum(result, switching_matrix)
        elif mode == 'average':
            weighted_sum += switching_matrix * prob
            prob_sum += prob
    
    # For 'average' mode, divide by sum of probabilities
    if mode == 'average':
        result = weighted_sum / prob_sum
    
    # Ensure diagonal is 0.0
    np.fill_diagonal(result, 0.0)
    
    # Create and return DistanceMatrix
    return DistanceMatrix(result, labels=all_taxa)


def induced_splits(network: DirectedPhyNetwork) -> SplitSystem:
    """
    Extract all splits induced by cut-edges of the network.
    
    This function:

    1. Converts the network to an LSA network
    2. Suppresses all 2-blobs (which don't influence splits)
    3. Finds all cut-edges
    4. For each cut-edge, computes the split it induces (2-partition of taxa)
    
    The split induced by a cut-edge is the 2-partition of taxa obtained when
    removing that edge from the network.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network.
    
    Returns
    -------
    SplitSystem
        A split system containing all splits induced by cut-edges.
    
    Examples
    --------
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> splits = induced_splits(net)
    >>> len(splits) >= 1
    True
    
    Notes
    -----
    The network is first converted to an LSA network and then the tree-of-blobs
    is computed. Since the tree-of-blobs has the same split system as the original
    network, we can efficiently compute splits using a single DFS traversal of
    the tree structure.
    """
    # Handle empty networks
    if network.number_of_nodes() == 0:
        return SplitSystem()
    
    # Optimization: If network is already a tree, use it directly
    if network.is_tree():
        blob_tree = network
    else:
        # Step 1: Convert to LSA network
        lsa_network = to_lsa_network(network) if not is_lsa_network(network) else network
        
        # Step 2: Get tree-of-blobs (this has the same splits as the original network)
        blob_tree = tree_of_blobs(lsa_network)
    
    # Get all taxa
    all_taxa = list(blob_tree.taxa)
    if len(all_taxa) < 2:
        # Need at least 2 taxa for splits
        return SplitSystem()
    
    # Step 3: Single DFS traversal to compute splits
    # In a tree, every edge is a cut-edge, so we can efficiently compute splits
    splits: set[Split] = set()
    
    # Use combined graph view (undirected) for traversal
    combined_graph = blob_tree._graph._combined
    
    # Pick an arbitrary leaf as starting point for DFS
    if len(blob_tree.leaves) == 0:
        return SplitSystem()
    
    start_leaf = next(iter(blob_tree.leaves))
    visited: set[Any] = set()
    
    def dfs(node: Any, parent: Any | None = None) -> set[str]:
        """DFS that returns the set of leaves in the subtree rooted at node."""
        visited.add(node)
        node_leaves: set[str] = set()
        
        # If this is a leaf, add its taxon
        if node in blob_tree._node_to_label and node in blob_tree.leaves:
            taxon = blob_tree._node_to_label[node]
            node_leaves.add(taxon)
        
        # Process neighbors (excluding parent)
        for neighbor in combined_graph.neighbors(node):
            if neighbor == parent:
                continue
            if neighbor not in visited:
                neighbor_leaves = dfs(neighbor, node)
                node_leaves.update(neighbor_leaves)
                
                # For edge (node, neighbor), create split
                # Split: (leaves in neighbor's subtree, all other leaves)
                if neighbor_leaves and len(neighbor_leaves) < len(all_taxa):
                    other_leaves = set(all_taxa) - neighbor_leaves
                    if other_leaves:  # Both sides must have at least one leaf
                        split = Split(neighbor_leaves, other_leaves)
                        splits.add(split)
        
        return node_leaves
    
    # Start DFS from the leaf
    dfs(start_leaf)
    
    return SplitSystem(splits)


def split_from_cutedge(
    network: DirectedPhyNetwork,
    u: Any,
    v: Any,
    key: int | None = None,
    return_node_taxa: bool = False,
) -> Split | tuple[Split, tuple[Any, frozenset[str]], tuple[Any, frozenset[str]]]:
    """
    Get the split induced by a cut-edge in the network.
    
    This function removes the specified edge from the network and finds the
    taxa on either side of the resulting partition. If the edge is not a
    cut-edge (i.e., removing it does not disconnect the graph), an error is raised.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network.
    u : Any
        First node of the edge.
    v : Any
        Second node of the edge.
    key : int | None, optional
        Edge key for parallel edges. If None and multiple parallel edges exist,
        raises PhyloZooValueError. If None and exactly one edge exists, that edge is used.
        By default None.
    return_node_taxa : bool, optional
        If True, also returns tuples (u, taxa1) and (v, taxa2) indicating which
        node is on which side of the split. By default False.
    
    Returns
    -------
    Split | tuple[Split, tuple[Any, frozenset[str]], tuple[Any, frozenset[str]]]
        If return_node_taxa is False: The split induced by the cut-edge.
        If return_node_taxa is True: A tuple (split, (u, taxa1), (v, taxa2)) where
        taxa1 are the taxa on the side of u and taxa2 are the taxa on the side of v.
    
    Raises
    ------
    PhyloZooValueError
        If the edge does not exist, if multiple parallel edges exist and key is None,
        or if the edge is not a cut-edge (removal does not disconnect the graph).
    
    Examples
    --------
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2), (3, 4)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
    ... )
    >>> split = split_from_cutedge(net, 3, 1)
    >>> 'A' in split.set1 or 'A' in split.set2
    True
    """
    # Create a copy of the underlying graph
    graph_copy = network._graph.copy()
    
    # Check if edge exists and handle keys
    if not graph_copy.has_edge(u, v):
        raise PhyloZooValueError(f"Edge ({u}, {v}) does not exist")
    
    # Get all edge keys between u and v
    edge_keys = list(graph_copy[u][v].keys())
    if len(edge_keys) == 0:
        raise PhyloZooValueError(f"Edge ({u}, {v}) does not exist")
    if len(edge_keys) > 1 and key is None:
        raise PhyloZooValueError(
            f"Multiple parallel edges exist between {u} and {v}. "
            "Must specify 'key' parameter."
        )
    if key is None:
        key = edge_keys[0]
    if key not in edge_keys:
        raise PhyloZooValueError(f"Edge ({u}, {v}, key={key}) does not exist")
    
    # Remove the edge
    graph_copy.remove_edge(u, v, key)
    
    # Check if removal disconnects the graph
    import networkx as nx
    components = list(nx.weakly_connected_components(graph_copy._graph))
    if len(components) != 2:
        raise PhyloZooValueError(
            f"Edge ({u}, {v}, key={key}) is not a cut-edge. "
            f"Removal creates {len(components)} components instead of 2."
        )
    
    # Determine which component contains u and v
    component_u = None
    component_v = None
    for comp in components:
        if u in comp:
            component_u = comp
        if v in comp:
            component_v = comp
    
    # Get taxa in each component (intersection of component with leaves)
    leaves_set = network.leaves
    taxa_u = frozenset(
        network._node_to_label[node]
        for node in component_u & leaves_set
        if node in network._node_to_label
    )
    taxa_v = frozenset(
        network._node_to_label[node]
        for node in component_v & leaves_set
        if node in network._node_to_label
    )
    
    split = Split(taxa_u, taxa_v)
    
    if return_node_taxa:
        return (split, (u, taxa_u), (v, taxa_v))
    return split


def displayed_splits(network: DirectedPhyNetwork) -> WeightedSplitSystem:
    """
    Compute weighted split system from all displayed trees of the network.
    
    This function iterates through all displayed trees of the network and collects
    their induced splits, weighted by the probability of each displayed tree. If a
    split appears in multiple displayed trees, their probabilities are summed.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network.
    
    Returns
    -------
    WeightedSplitSystem
        A weighted split system where each split's weight is the sum of probabilities
        of all displayed trees that contain that split.
    
    Examples
    --------
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         (10, 5), (10, 6),  # Root to tree nodes
    ...         (5, 4), (6, 4),    # Both lead to hybrid 4 (in-degree 2)
    ...         (4, 8),            # Hybrid to tree node
    ...         (8, 1), (8, 2),    # Tree node to leaves
    ...         (5, 3), (6, 7)     # Additional leaves to satisfy degree constraints
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
    ... )
    >>> splits = displayed_splits(net)
    >>> isinstance(splits, WeightedSplitSystem)
    True
    >>> len(splits) > 0
    True
    """
    # Handle empty networks
    if network.number_of_nodes() == 0:
        return WeightedSplitSystem()
    
    # Collect splits with their weights
    split_weights: dict[Split, float] = {}
    
    # Iterate through all displayed trees with probabilities
    for displayed_tree in displayed_trees(network, probability=True):
        # Get probability of this displayed tree
        prob = displayed_tree.get_network_attribute('probability')
        if prob is None:
            prob = 1.0
        
        # Get induced splits from this displayed tree
        tree_splits = induced_splits(displayed_tree)
        
        # Add each split with its weight (accumulate if split already exists)
        for split in tree_splits.splits:
            split_weights[split] = split_weights.get(split, 0.0) + prob
    
    # Create weighted split system
    if not split_weights:
        return WeightedSplitSystem()
    
    return WeightedSplitSystem(split_weights)


def displayed_quartets(network: DirectedPhyNetwork) -> QuartetProfileSet:
    """
    Compute quartet profile set from all displayed trees of the network.
    
    This function converts the directed network to a semi-directed network and then
    uses the semi-directed displayed_quartets function. This ensures that quartets
    are unrooted (not rooted quartets).
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network.
    
    Returns
    -------
    QuartetProfileSet
        A quartet profile set where each profile corresponds to a 4-taxon set,
        and contains quartets from displayed trees weighted by their probabilities.
    
    Examples
    --------
    >>> net = DirectedPhyNetwork(
    ...     edges=[(5, 4), (6, 4), (4, 1), (4, 2), (5, 3), (6, 7)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
    ... )
    >>> profileset = displayed_quartets(net)
    >>> isinstance(profileset, QuartetProfileSet)
    True
    >>> len(profileset) > 0
    True
    """
    # Convert to semi-directed network
    sd_network = to_sd_network(network)
    
    # Import here to avoid circular dependency
    from ..sdnetwork.derivations import displayed_quartets as sd_displayed_quartets
    
    # Use the semi-directed displayed_quartets function
    return sd_displayed_quartets(sd_network)


def partition_from_blob(
    network: DirectedPhyNetwork,
    blob: set[Any],
    return_edge_taxa: bool = False,
) -> Partition | tuple[Partition, list[tuple[Any, Any, frozenset[str]]]]:
    """
    Get the partition of taxa induced by removing a blob from the network.
    
    When all nodes in the blob are removed, the network splits into connected components.
    Each component's taxa form a part of the partition.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network.
    blob : set[Any]
        Set of nodes forming the blob to remove. All nodes in this set will be
        removed to compute the partition.
    return_edge_taxa : bool, optional
        If True, also return a list of tuples (u, v, taxa_set) where u is a node
        in the component, v is a node in the blob, and taxa_set is the frozenset
        of taxa in that component. By default False.
    
    Returns
    -------
    Partition | tuple[Partition, list[tuple[Any, Any, frozenset[str]]]]
        If return_edge_taxa is False: The partition of taxa induced by removing the blob.
        If return_edge_taxa is True: A tuple (partition, edge_taxa_list) where edge_taxa_list
        is a list of (u, v, taxa_set) tuples connecting each component to the blob.
    
    Raises
    ------
    PhyloZooValueError
        If blob is empty or if blob contains nodes not in the network.
        If blob is not a non-leaf blob (internal blob).
    PhyloZooAlgorithmError
        If could not find edge connecting component with taxa to blob.
        If removing the blob does not disconnect the network (blob is not a cut-blob).
    Examples
    --------
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2), (3, 4)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
    ... )
    >>> partition = partition_from_blob(net, {3})
    >>> len(partition)
    3
    """
    # Validate blob
    if not blob:
        raise PhyloZooValueError("Blob cannot be empty")
    
    # Check all nodes in blob are in network
    network_nodes = set(network._graph.nodes)
    missing_nodes = blob - network_nodes
    if missing_nodes:
        raise PhyloZooValueError(
            f"Blob contains nodes not in network: {missing_nodes}"
        )
    
    # Check that blob is a non-leaf blob (including trivial single-node blobs)
    non_leaf_blobs = blobs(network, trivial=True, leaves=False)
    blob_frozen = frozenset(blob)
    if blob_frozen not in {frozenset(b) for b in non_leaf_blobs}:
        raise PhyloZooValueError(
            f"Blob {blob} is not a non-leaf blob. "
            "Only non-leaf blobs (internal blobs) can be used for partition_from_blob."
        )
    
    # Create a copy of the graph and remove all blob nodes
    graph_copy = network._graph.copy()
    for node in blob:
        graph_copy.remove_node(node)
    
    # Find weakly connected components (for directed graphs)
    import networkx as nx
    components = list(nx.weakly_connected_components(graph_copy._graph))
    
    # Check that removing blob disconnects the network (at least 2 components)
    if len(components) < 2:
        raise PhyloZooAlgorithmError(
            f"Removing blob {blob} does not disconnect the network. "
            f"Result has {len(components)} component(s), expected at least 2."
        )
    
    # Get leaves set for taxon extraction
    leaves_set = network.leaves
    
    # Build partition parts and edge_taxa list
    partition_parts: list[set[str]] = []
    edge_taxa_list: list[tuple[Any, Any, frozenset[str]]] = []
    
    # For each component, find taxa and connecting edge
    for component in components:
        # Get taxa in this component
        component_taxa = frozenset(
            network._node_to_label[node]
            for node in component & leaves_set
            if node in network._node_to_label
        )
        
        # Skip components with no taxa
        if not component_taxa:
            continue
        
        partition_parts.append(set(component_taxa))
        
        # If return_edge_taxa is True, find an edge connecting this component to the blob
        if return_edge_taxa:
            # Find a node in the component that has an edge to a node in the blob
            u = None
            v = None
            for node_in_component in component:
                # Check for edges to blob nodes
                for blob_node in blob:
                    if network._graph.has_edge(node_in_component, blob_node):
                        u = node_in_component
                        v = blob_node
                        break
                    elif network._graph.has_edge(blob_node, node_in_component):
                        u = node_in_component
                        v = blob_node
                        break
                if u is not None:
                    break
            
            if u is not None and v is not None:
                edge_taxa_list.append((u, v, component_taxa))
            else:
                # Should not happen if blob is a cut-blob, but handle gracefully
                raise PhyloZooAlgorithmError(
                    f"Could not find edge connecting component with taxa {component_taxa} to blob"
                )
    
    # Ensure partition covers all taxa (add missing as singletons)
    partition_taxa = set()
    for part in partition_parts:
        partition_taxa.update(part)
    missing_taxa = network.taxa - partition_taxa
    for taxon in missing_taxa:
        partition_parts.append({taxon})
        # For missing taxa, the leaf node might be in the blob itself
        # Find the leaf node and an edge connecting it to outside the blob
        if return_edge_taxa:
            taxon_node = network.get_node_id(taxon)
            if taxon_node is not None and taxon_node in blob:
                # Find an edge from this leaf (in blob) to a node outside blob
                # or from outside blob to this leaf
                u = None
                v = None
                # Check outgoing edges
                for neighbor in network._graph.neighbors(taxon_node):
                    if neighbor not in blob:
                        u = taxon_node
                        v = neighbor
                        break
                # Check incoming edges
                if u is None:
                    for predecessor in network._graph.predecessors(taxon_node):
                        if predecessor not in blob:
                            u = predecessor
                            v = taxon_node
                            break
                
                if u is not None and v is not None:
                    edge_taxa_list.append((u, v, frozenset({taxon})))
                else:
                    # Leaf is in blob but has no connection outside - shouldn't happen in valid networks
                    # but handle gracefully by using the leaf node itself
                    edge_taxa_list.append((taxon_node, taxon_node, frozenset({taxon})))
    
    partition = Partition(partition_parts)
    
    if return_edge_taxa:
        return (partition, edge_taxa_list)
    return partition
