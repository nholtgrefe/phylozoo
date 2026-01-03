"""
Network derivations module.

This module provides functions to derive other data structures from semi-directed
and mixed phylogenetic networks (e.g., splits, quartets, distances, blobtrees, subnetworks, etc.).
"""

# TODO: Implement derivation functions here.
# These could include:
# - Split extraction
# - Quartet extraction
# - Distance calculations
# - Blobtree construction
# - Subnetwork extraction
# - Displayed tree extraction
import itertools
from typing import Any, Iterator, Literal

import numpy as np
import networkx as nx

from . import MixedPhyNetwork
from .classifications import is_tree
from .features import blobs
from .transformations import suppress_2_blobs as suppress_2_blobs_fn, identify_parallel_edges as identify_parallel_edges_fn
from ...split import Split, SplitSystem, WeightedSplitSystem
from .sd_phynetwork import SemiDirectedPhyNetwork
from ._utils import _suppress_deg2_nodes
from .conversions import sdnetwork_from_graph
from ...primitives.m_multigraph.transformations import (
    identify_vertices as mm_identify_vertices,
    suppress_degree2_node as mm_suppress_degree2_node,
    subgraph as mm_subgraph,
)
from ...primitives.m_multigraph.features import updown_path_vertices
from ...primitives.m_multigraph import MixedMultiGraph
from .conversions import sdnetwork_from_graph
from ....core.distance import DistanceMatrix


def tree_of_blobs(network: MixedPhyNetwork) -> MixedPhyNetwork:
    """
    Create the tree-of-blobs by suppressing all 2-blobs and collapsing internal blobs.

    This function:
    1. Suppresses all 2-blobs using suppress_2_blobs
    2. Finds all internal blobs (blobs with more than 1 node, excluding leaves)
    3. For each internal blob, identifies all vertices with a single vertex
    4. Returns a new network representing the tree-of-blobs

    Parameters
    ----------
    network : MixedPhyNetwork
        The mixed phylogenetic network to transform.

    Returns
    -------
    MixedPhyNetwork | SemiDirectedPhyNetwork
        A new network where each blob has been collapsed to a single vertex,
        forming a tree structure. Returns a SemiDirectedPhyNetwork if the input
        is a SemiDirectedPhyNetwork, otherwise returns a MixedPhyNetwork.

    Examples
    --------
    >>> # Create a semi-directed network with a hybrid
    >>> from phylozoo.core.network.sdnetwork.classifications import is_tree
    >>> sdnet = SemiDirectedPhyNetwork(
    ...     directed_edges=[
    ...         (5, 4),
    ...         (6, 4)
    ...     ],
    ...     undirected_edges=[
    ...         (5, 3),
    ...         (5, 6),
    ...         (6, 7),
    ...         (4, 8),
    ...         (8, 1),
    ...         (8, 2)
    ...     ],
    ...     nodes=[
    ...         (3, {'label': 'C'}),
    ...         (7, {'label': 'D'}),
    ...         (1, {'label': 'A'}),
    ...         (2, {'label': 'B'})
    ...     ]
    ... )
    >>> tree_net = tree_of_blobs(sdnet)
    >>> is_tree(tree_net)
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
            mm_identify_vertices(working_graph, blob_sorted)

    # Convert back to appropriate network type
    # Preserve the input type (SemiDirectedPhyNetwork or MixedPhyNetwork)
    network_type = 'semi-directed' if isinstance(network, SemiDirectedPhyNetwork) else 'mixed'
    return sdnetwork_from_graph(working_graph, network_type=network_type)


def subnetwork(
    network: SemiDirectedPhyNetwork,
    taxa: list[str],
    suppress_2_blobs: bool = False,
    identify_parallel_edges: bool = False,
) -> SemiDirectedPhyNetwork:
    """
    Extract the subnetwork induced by a subset of taxa (leaf labels).

    The subnetwork is defined as the union of all up-down paths between the
    requested leaves (i.e., all vertices on up-down paths between any pair
    of the requested leaves). The induced subgraph is taken on the underlying
    MixedMultiGraph, then degree-2 internal nodes are suppressed. Optionally,
    the result can be post-processed by suppressing 2-blobs and/or identifying
    parallel edges. After any of these optional steps, degree-2 suppression is
    applied again to clean up artifacts.

    An up-down path between two vertices x and y is a path where no two edges
    are oriented towards each other. Equivalently, it is a path where the first
    k edges can be oriented towards x (where undirected edges can be oriented in
    either way) and the remaining l-k edges can be oriented towards y.

    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        Source network.
    taxa : list[str]
        Subset of taxon labels (leaf labels) to induce the subnetwork on.

    Other parameters (keyword-only)
    -------------------------------
    suppress_2_blobs : bool, default False
        If True, suppress all 2-blobs in the resulting network.
    identify_parallel_edges : bool, default False
        If True, identify/merge parallel edges in the resulting network.

    Returns
    -------
    SemiDirectedPhyNetwork
        The derived subnetwork. Returns an empty network if `taxa` is empty.
        Returns a network with a single leaf if `taxa` contains a single leaf.

    Raises
    ------
    ValueError
        If any of the provided taxa are not found in the network.

    Examples
    --------
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 4)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
    ... )
    >>> subnet = subnetwork(net, ['A', 'B'])
    >>> sorted(subnet.taxa)
    ['A', 'B']
    >>> # Network with hybrid
    >>> net2 = SemiDirectedPhyNetwork(
    ...     directed_edges=[(5, 4), (6, 4)],
    ...     undirected_edges=[(5, 3), (5, 6), (6, 7), (4, 8), (8, 1), (8, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
    ... )
    >>> subnet2 = subnetwork(net2, ['A', 'B'])
    >>> sorted(subnet2.taxa)
    ['A', 'B']
    """
    if not taxa:
        # Return an empty network if no taxa are specified
        return SemiDirectedPhyNetwork(directed_edges=[], undirected_edges=[], nodes=[])
    
    # Map taxa labels to node ids using public API
    leaf_nodes: list[Any] = []
    for t in taxa:
        node_id = network.get_node_id(t)
        if node_id is None:
            raise ValueError(f"Taxon label '{t}' not found in network")
        leaf_nodes.append(node_id)
    
    # Collect all vertices on up-down paths between any pair of leaves
    # Note: updown_path_vertices handles the single-leaf case (returns {x} when x == y)
    nodes_set: set[Any] = set()
    if len(leaf_nodes) == 1:
        # Single leaf case: updown_path_vertices(leaf, leaf) returns {leaf}
        nodes_set = updown_path_vertices(network._graph, leaf_nodes[0], leaf_nodes[0])
    else:
        for leaf1, leaf2 in itertools.combinations(leaf_nodes, 2):
            # Find vertices on up-down paths between this pair of leaves
            path_vertices = updown_path_vertices(network._graph, leaf1, leaf2)
            nodes_set.update(path_vertices)
        # Also include all leaves themselves (they should already be included, but ensure)
        nodes_set.update(leaf_nodes)
    
    # Create induced MixedMultiGraph using existing utility
    induced_mm = mm_subgraph(network._graph, nodes_set)
    
    # Work on a mutable copy for transformations
    working_mm = induced_mm.copy()
    
    # First pass: suppress all degree-2 nodes (excluding leaves)
    leaf_set = set(leaf_nodes)
    _suppress_deg2_nodes(working_mm, exclude_nodes=leaf_set)
    
    # Convert to SemiDirectedPhyNetwork for higher-level transformations
    result_net = sdnetwork_from_graph(working_mm, network_type='semi-directed')
    
    # Optional post-processing steps
    if suppress_2_blobs:
        result_net = suppress_2_blobs_fn(result_net)
    
    if identify_parallel_edges:
        result_net = identify_parallel_edges_fn(result_net)
    
    return result_net


def k_taxon_subnetworks(
    network: SemiDirectedPhyNetwork,
    k: int,
    suppress_2_blobs: bool = False,
    identify_parallel_edges: bool = False,
) -> Iterator[SemiDirectedPhyNetwork]:
    """
    Generate all subnetworks induced by exactly k taxa.
    
    This function yields all possible subnetworks of the network that are
    induced by exactly k taxon labels. For each combination of k taxa,
    the corresponding subnetwork is computed using the `subnetwork` function.
    
    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        Source network.
    k : int
        Number of taxa to include in each subnetwork. Must be between 0 and
        the number of taxa in the network (inclusive).
    
    Other parameters (keyword-only)
    -------------------------------
    suppress_2_blobs : bool, default False
        If True, suppress all 2-blobs in each resulting subnetwork.
    identify_parallel_edges : bool, default False
        If True, identify/merge parallel edges in each resulting subnetwork.
    
    Yields
    ------
    SemiDirectedPhyNetwork
        Subnetworks induced by exactly k taxa. Each subnetwork is generated
        lazily as the iterator is consumed.
    
    Raises
    ------
    ValueError
        If k < 0 or k > number of taxa in the network.
    
    Examples
    --------
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(5, 3), (5, 4), (5, 6), (3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (6, {'label': 'D'})]
    ... )
    >>> # Generate all 2-taxon subnetworks
    >>> subnetworks = list(k_taxon_subnetworks(net, k=2))
    >>> len(subnetworks)
    6  # C(4,2) = 6 combinations
    >>> # Each subnetwork has exactly 2 leaves
    >>> all(len(subnet.taxa) == 2 for subnet in subnetworks)
    True
    >>> # Generate all 1-taxon subnetworks
    >>> single_taxon_subs = list(k_taxon_subnetworks(net, k=1))
    >>> len(single_taxon_subs)
    4  # C(4,1) = 4 combinations
    """
    all_taxa = list(network.taxa)
    num_taxa = len(all_taxa)
    
    # Validate k
    if k < 0:
        raise ValueError(f"k must be non-negative, got {k}")
    if k > num_taxa:
        raise ValueError(
            f"k ({k}) cannot exceed the number of taxa ({num_taxa}) in the network"
        )
    
    # Generate all combinations of k taxa
    for taxa_combination in itertools.combinations(all_taxa, k):
        yield subnetwork(
            network,
            list(taxa_combination),
            suppress_2_blobs=suppress_2_blobs,
            identify_parallel_edges=identify_parallel_edges,
        )


def _switchings(network: SemiDirectedPhyNetwork, probability: bool = False) -> Iterator[MixedMultiGraph]:
    """
    Generate all switchings of a semi-directed phylogenetic network.
    
    A switching is obtained by deleting all but one incident parent edge for each
    hybrid node. This function generates all possible combinations of keeping
    exactly one parent edge per hybrid node. Each switching is a tree (not
    necessarily a phylogenetic tree).
    
    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network.
    probability : bool, optional
        If True, store the probability of the switching in the graph's 'probability'
        attribute. The probability is the product of gamma values for the kept hybrid
        edges. If a hybrid edge has no gamma value, it is taken to be 1/k where k
        is the in-degree of the hybrid node. By default False.
    
    Yields
    ------
    MixedMultiGraph
        A switching of the network (one parent edge kept per hybrid node). Each
        switching is a tree. If probability=True, the graph has a 'probability'
        attribute containing the switching probability.
    
    Examples
    --------
    >>> net = SemiDirectedPhyNetwork(
    ...     directed_edges=[
    ...         (5, 4),
    ...         (6, 4)
    ...     ],
    ...     undirected_edges=[
    ...         (5, 3),
    ...         (5, 6),
    ...         (6, 7),
    ...         (4, 8),
    ...         (8, 1),
    ...         (8, 2)
    ...     ],
    ...     nodes=[
    ...         (3, {'label': 'C'}),
    ...         (7, {'label': 'D'}),
    ...         (1, {'label': 'A'}),
    ...         (2, {'label': 'B'})
    ...     ]
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
    >>> net_with_gamma = SemiDirectedPhyNetwork(
    ...     directed_edges=[
    ...         {'u': 5, 'v': 4, 'gamma': 0.6},
    ...         {'u': 6, 'v': 4, 'gamma': 0.4}
    ...     ],
    ...     undirected_edges=[
    ...         (5, 3), (5, 6), (6, 7), (4, 8), (8, 1), (8, 2)
    ...     ],
    ...     nodes=[
    ...         (3, {'label': 'C'}), (7, {'label': 'D'}),
    ...         (1, {'label': 'A'}), (2, {'label': 'B'})
    ...     ]
    ... )
    >>> switchings_with_prob = list(_switchings(net_with_gamma, probability=True))
    >>> switchings_with_prob[0]._directed.graph.get('probability')
    0.6  # Probability of keeping edge (5,4)
    >>> switchings_with_prob[1]._directed.graph.get('probability')
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


def displayed_trees(network: SemiDirectedPhyNetwork, probability: bool = False) -> Iterator[SemiDirectedPhyNetwork]:
    """
    Generate all displayed trees of a semi-directed phylogenetic network.
    
    A displayed tree is obtained by:
    1. Taking a switching (deleting all but one parent edge per hybrid node)
    2. Exhaustively removing degree-1 nodes that are not leaves
    3. Suppressing all degree-2 nodes
    
    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network.
    probability : bool, optional
        If True, store the probability of the displayed tree in the network's
        'probability' attribute. The probability is inherited from the switching
        and equals the product of gamma values for the kept hybrid edges. If a
        hybrid edge has no gamma value, it is taken to be 1/k where k is the
        in-degree of the hybrid node. If there are no hybrid nodes, the
        probability is 1.0. By default False.
    
    Yields
    ------
    SemiDirectedPhyNetwork
        A displayed tree of the network. If probability=True, the network has a
        'probability' attribute containing the tree probability.
    
    Examples
    --------
    >>> net = SemiDirectedPhyNetwork(
    ...     directed_edges=[(5, 4), (6, 4)],
    ...     undirected_edges=[(5, 3), (5, 6), (6, 7), (4, 8), (8, 1), (8, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
    ... )
    >>> trees = list(displayed_trees(net))
    >>> len(trees)
    2  # Two switchings yield two displayed trees
    """
    # Get original network's leaves for reference
    original_leaves = network.leaves
    
    # Iterate through all switchings
    for switching_graph in _switchings(network, probability=probability):
        # Work on a copy to avoid modifying the switching
        tree_graph = switching_graph.copy()
        
        # Exhaustively remove degree-1 nodes that are not leaves
        while True:
            degree1_nodes = [
                node for node in tree_graph.nodes()
                if tree_graph.degree(node) == 1
                and node not in original_leaves
            ]
            
            if not degree1_nodes:
                break
            
            # Remove all degree-1 nodes (excluding leaves)
            for node in degree1_nodes:
                # Double-check node still exists and is still degree-1
                if node in tree_graph.nodes() and tree_graph.degree(node) == 1:
                    tree_graph.remove_node(node)
        
        # Suppress all degree-2 nodes
        _suppress_deg2_nodes(tree_graph, exclude_nodes=None)
        
        # Convert back to SemiDirectedPhyNetwork
        # Note: sdnetwork_from_graph already copies graph attributes, so probability
        # is automatically preserved from the switching graph.
        displayed_tree = sdnetwork_from_graph(tree_graph, network_type='semi-directed')
        
        yield displayed_tree


def _switching_distance_matrix(
    switching_graph: MixedMultiGraph,
    taxa: list[str],
    original_network: SemiDirectedPhyNetwork,
) -> np.ndarray:
    """
    Compute the full distance matrix for all pairwise distances between taxa in a switching.
    
    Parameters
    ----------
    switching_graph : MixedMultiGraph
        The switching graph (a tree).
    taxa : list[str]
        List of taxon labels.
    original_network : SemiDirectedPhyNetwork
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
            raise ValueError(f"Taxon '{taxon}' not found in network")
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
    network: SemiDirectedPhyNetwork,
    mode: Literal['shortest', 'longest', 'average'] = 'average',
) -> DistanceMatrix:
    """
    Compute pairwise distances between taxa based on switchings.
    
    This function computes distances by considering all switchings of the network.
    For each pair of taxa, the distance is computed in each switching (sum of branch
    lengths along the unique path), and then aggregated according to the specified mode.
    
    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network.
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
    >>> net = SemiDirectedPhyNetwork(
    ...     directed_edges=[(5, 4), (6, 4)],
    ...     undirected_edges=[(5, 3), (5, 6), (6, 7), (4, 8), (8, 1), (8, 2)],
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
        raise ValueError(f"Invalid mode: {mode}. Must be 'shortest', 'longest', or 'average'")
    
    # Iterate through all switchings
    for switching_graph in _switchings(network, probability=(mode == 'average')):
        # Get probability if mode is 'average'
        prob = 1.0
        if mode == 'average':
            prob = switching_graph._directed.graph.get('probability', 1.0)
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


def induced_splits(network: MixedPhyNetwork) -> SplitSystem:
    """
    Extract all splits induced by cut-edges of the network.
    
    This function:
    1. Suppresses all 2-blobs (which don't influence splits)
    2. Finds all cut-edges
    3. For each cut-edge, computes the split it induces (2-partition of taxa)
    
    The split induced by a cut-edge is the 2-partition of taxa obtained when
    removing that edge from the network.
    
    Parameters
    ----------
    network : MixedPhyNetwork
        The mixed phylogenetic network.
    
    Returns
    -------
    SplitSystem
        A split system containing all splits induced by cut-edges.
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> splits = induced_splits(net)
    >>> len(splits) >= 1
    True
    
    Notes
    -----
    The tree-of-blobs is computed, which has the same split system as the original
    network. We can efficiently compute splits using a single DFS traversal of
    the tree structure.
    """
    # Handle empty networks
    if network.number_of_nodes() == 0:
        return SplitSystem()
    
    # Optimization: If network is already a tree, use it directly
    if is_tree(network):
        blob_tree = network
    else:
        # Step 1: Get tree-of-blobs (this has the same splits as the original network)
        blob_tree = tree_of_blobs(network)
    
    # Get all taxa
    all_taxa = list(blob_tree.taxa)
    if len(all_taxa) < 2:
        # Need at least 2 taxa for splits
        return SplitSystem()
    
    # Step 2: Single DFS traversal to compute splits
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


def displayed_splits(network: SemiDirectedPhyNetwork) -> WeightedSplitSystem:
    """
    Compute weighted split system from all displayed trees of the network.
    
    This function iterates through all displayed trees of the network and collects
    their induced splits, weighted by the probability of each displayed tree. If a
    split appears in multiple displayed trees, their probabilities are summed.
    
    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network.
    
    Returns
    -------
    WeightedSplitSystem
        A weighted split system where each split's weight is the sum of probabilities
        of all displayed trees that contain that split.
    
    Examples
    --------
    >>> net = SemiDirectedPhyNetwork(
    ...     directed_edges=[(5, 4), (6, 4)],  # Hybrid edges to hybrid node 4
    ...     undirected_edges=[
    ...         (5, 3), (5, 6), (6, 7),  # Tree edges
    ...         (4, 8), (8, 1), (8, 2)   # Tree edges from hybrid
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
    
    return WeightedSplitSystem(list(split_weights.keys()), weights=split_weights)
