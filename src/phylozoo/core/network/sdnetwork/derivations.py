"""
Network derivations module.

This module provides functions to derive other data structures from semi-directed
and mixed phylogenetic networks (e.g., splits, quartets, distances, blobtrees, subnetworks, etc.).
"""

import itertools
from typing import TYPE_CHECKING, Any, Iterator, Literal

import numpy as np
import networkx as nx

if TYPE_CHECKING:
    from ...network.dnetwork import DirectedPhyNetwork

from . import MixedPhyNetwork
from .classifications import is_tree
from .features import blobs, root_locations, RootLocation
from .transformations import suppress_2_blobs as suppress_2_blobs_fn, identify_parallel_edges as identify_parallel_edges_fn
from ...split import Split, SplitSystem, WeightedSplitSystem
from ...quartet import Quartet, QuartetProfile, QuartetProfileSet
from ...primitives.partition import Partition
from .sd_phynetwork import SemiDirectedPhyNetwork
from ._utils import _suppress_deg2_nodes, _subdivide_edge
from .conversions import sdnetwork_from_graph
from ...primitives.m_multigraph.transformations import (
    identify_vertices as mm_identify_vertices,
    suppress_degree2_node as mm_suppress_degree2_node,
    subgraph as mm_subgraph,
    orient_away_from_vertex,
)
from ...primitives.m_multigraph.features import updown_path_vertices
from ...primitives.m_multigraph import MixedMultiGraph
from .conversions import sdnetwork_from_graph
from ....core.distance import DistanceMatrix
from ....utils.exceptions import PhyloZooValueError, PhyloZooError


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
    PhyloZooValueError
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
            raise PhyloZooValueError(f"Taxon label '{t}' not found in network")
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
    PhyloZooValueError
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
    
    Raises
    ------
    PhyloZooValueError
        If the mode is invalid.
    
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


def split_from_cutedge(
    network: SemiDirectedPhyNetwork,
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
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network.
    u : Any
        First node of the edge.
    v : Any
        Second node of the edge.
    key : int | None, optional
        Edge key for parallel edges. If None and multiple parallel edges exist,
        raises ValueError. If None and exactly one edge exists, that edge is used.
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
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 4)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
    ... )
    >>> split = split_from_cutedge(net, 3, 1)
    >>> 'A' in split.set1 or 'A' in split.set2
    True
    """
    # Create a copy of the underlying graph
    graph_copy = network._graph.copy()
    
    # Check if edge exists and handle keys
    if graph_copy._undirected.has_edge(u, v):
        # Get all edge keys between u and v
        edge_keys = list(graph_copy._undirected[u][v].keys())
        if len(edge_keys) == 0:
            raise PhyloZooValueError(f"Undirected edge ({u}, {v}) does not exist")
        elif len(edge_keys) > 1 and key is None:
            raise PhyloZooValueError(
                f"Multiple parallel undirected edges exist between {u} and {v}. "
                "Must specify 'key' parameter."
            )
        elif key is None:
            key = edge_keys[0]
        elif key not in edge_keys:
            raise PhyloZooValueError(f"Edge ({u}, {v}, key={key}) does not exist")
        
        # Remove the edge using public API
        graph_copy.remove_edge(u, v, key=key)

    elif graph_copy._directed.has_edge(u, v):
        # Get all edge keys between u and v
        edge_keys = list(graph_copy._directed[u][v].keys())
        if len(edge_keys) == 0:
            raise PhyloZooValueError(f"Directed edge ({u}, {v}) does not exist")
        elif len(edge_keys) > 1 and key is None:
            raise PhyloZooValueError(
                f"Multiple parallel directed edges exist between {u} and {v}. "
                "Must specify 'key' parameter."
            )
        elif key is None:
            key = edge_keys[0]
        elif key not in edge_keys:
            raise PhyloZooValueError(f"Edge ({u}, {v}, key={key}) does not exist")
        # Remove the edge using public API
        graph_copy.remove_edge(u, v, key=key)
    elif graph_copy._directed.has_edge(v, u):
        # Get all edge keys between v and u
        edge_keys = list(graph_copy._directed[v][u].keys())
        if len(edge_keys) == 0:
            raise PhyloZooValueError(f"Edge ({u}, {v}) does not exist")
        elif len(edge_keys) > 1 and key is None:
            raise PhyloZooValueError(
                f"Multiple parallel directed edges exist between {v} and {u}. "
                "Must specify 'key' parameter."
            )
        elif key is None:
            key = edge_keys[0]
        elif key not in edge_keys:
            raise PhyloZooValueError(f"Edge ({v}, {u}, key={key}) does not exist")
        # Remove the edge using public API (note: direction is v->u)
        graph_copy.remove_edge(v, u, key=key)
    else:
        raise PhyloZooValueError(f"Edge ({u}, {v}) does not exist in the network")
    
    # Check if removal disconnects the graph
    components = list(nx.connected_components(graph_copy._combined))
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
    
    return WeightedSplitSystem(split_weights)


def displayed_quartets(network: SemiDirectedPhyNetwork) -> QuartetProfileSet:
    """
    Compute quartet profile set from all displayed trees of the network.
    
    For each quartet (4-leaf subnetwork), this function:
    1. Extracts the subnetwork induced by those 4 taxa
    2. Gets all displayed trees of that subnetwork (with probabilities)
    3. Converts each displayed tree to a quartet (4-leaf tree)
    4. Creates a quartet profile where each quartet's weight is the probability
       of the displayed tree that induced it (summing weights if the same quartet
       appears in multiple displayed trees)
    
    The profiles are then returned as a QuartetProfileSet, where each profile
    (one per quartet) has no weight (default weight 1.0).
    
    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network.
    
    Returns
    -------
    QuartetProfileSet
        A quartet profile set where each profile corresponds to a 4-taxon set,
        and contains quartets from displayed trees weighted by their probabilities.
    
    Examples
    --------
    >>> net = SemiDirectedPhyNetwork(
    ...     directed_edges=[(5, 4), (6, 4)],
    ...     undirected_edges=[
    ...         (5, 3), (5, 6), (6, 7),  # Tree edges
    ...         (4, 8), (8, 1), (8, 2)   # Tree edges from hybrid
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
    ... )
    >>> profileset = displayed_quartets(net)
    >>> isinstance(profileset, QuartetProfileSet)
    True
    >>> len(profileset) > 0
    True
    """
    # Handle networks with fewer than 4 taxa
    taxa_list = sorted(network.taxa)
    if len(taxa_list) < 4:
        return QuartetProfileSet()
    
    # Collect profiles for each 4-taxon set
    profiles: list[QuartetProfile] = []
    
    # Iterate through all combinations of 4 taxa
    for four_taxa in itertools.combinations(taxa_list, 4):
        four_taxa_set = frozenset(four_taxa)
        
        # Get subnetwork induced by these 4 taxa
        quartet_subnet = subnetwork(network, list(four_taxa))
        
        # Collect quartets with their weights for this 4-taxon set
        quartet_weights: dict[Quartet, float] = {}
        
        # Get all displayed trees of the subnetwork with probabilities
        for displayed_tree in displayed_trees(quartet_subnet, probability=True):
            # Get probability of this displayed tree
            prob = displayed_tree.get_network_attribute('probability')
            if prob is None:
                prob = 1.0
            
            # Extract quartet from the displayed tree (4-leaf tree)
            # For a 4-leaf tree, induced_splits gives us the splits
            # A resolved 4-leaf tree has exactly one non-trivial split (2|2)
            # A star 4-leaf tree has no non-trivial splits
            tree_splits = induced_splits(displayed_tree)
            
            # Find the non-trivial split (2|2 split) if it exists
            quartet_split: Split | None = None
            for split in tree_splits.splits:
                if not split.is_trivial() and len(split.set1) == 2 and len(split.set2) == 2:
                    quartet_split = split
                    break
            
            # Create quartet from split or as star tree
            if quartet_split is not None:
                quartet = Quartet(quartet_split)
            else:
                # Star tree: no non-trivial split, all 4 taxa are equivalent
                quartet = Quartet(four_taxa_set)
            
            # Add quartet to profile with its weight (sum if duplicate)
            quartet_weights[quartet] = quartet_weights.get(quartet, 0.0) + prob
        
        # Create profile for this 4-taxon set (only if we have quartets)
        if quartet_weights:
            profile = QuartetProfile(quartet_weights)
            profiles.append(profile)
    
    # Create QuartetProfileSet (each profile gets default weight 1.0, i.e., no weight)
    return QuartetProfileSet(profiles=profiles)


def _root_sd_network_at(
    network: SemiDirectedPhyNetwork,
    root_location: RootLocation,
) -> 'DirectedPhyNetwork':
    """
    Root a semi-directed network at the specified location.
    
    This is an internal helper function that performs the actual rooting operation.
    It subdivides edges if needed, orients the graph, and converts to a directed network.
    
    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network to root.
    root_location : RootLocation
        The root location. Can be:
        - A node (T): node in the network
        - An edge (tuple[T, T, int]): edge in the network as (u, v, key)
    
    Returns
    -------
    DirectedPhyNetwork
        A directed phylogenetic network rooted at the specified location.
    
    Raises
    ------
    PhyloZooValueError
        If the root location is not in the network, or if rooting/orientation fails.
    """
    # Local import to avoid circular dependencies
    from ...network.dnetwork.conversions import dnetwork_from_graph
    
    # Step 1: Handle root location and get the graph
    root_vertex: Any
    
    if isinstance(root_location, tuple) and len(root_location) == 3:
        # root_location is an edge: (u, v, key)
        u, v, key = root_location
        
        # Use _subdivide_edge helper to subdivide the edge
        # This returns a new MixedMultiGraph with the subdivided edge
        graph_copy, subdiv_node = _subdivide_edge(network, u, v, key)
        
        # Root vertex is the subdivision node
        root_vertex = subdiv_node
    else:
        # root_location is a node - copy the graph
        graph_copy = network._graph.copy()
        if root_location not in graph_copy.nodes():
            raise PhyloZooValueError(
                f"Node {root_location} not found in the network"
            )
        root_vertex = root_location
    
    # Step 2: Filter attributes
    # Remove all graph attributes
    graph_copy._directed.graph.clear()
    graph_copy._undirected.graph.clear()
    
    # Filter node attributes: only keep 'label'
    allowed_node_attrs = {'label'}
    for node in graph_copy.nodes():
        if node in graph_copy._directed.nodes:
            for attr_key in list(graph_copy._directed.nodes[node].keys()):
                if attr_key not in allowed_node_attrs:
                    del graph_copy._directed.nodes[node][attr_key]
        if node in graph_copy._undirected.nodes:
            for attr_key in list(graph_copy._undirected.nodes[node].keys()):
                if attr_key not in allowed_node_attrs:
                    del graph_copy._undirected.nodes[node][attr_key]
    
    # Filter edge attributes: only keep 'gamma' and 'branch_length'
    allowed_edge_attrs = {'gamma', 'branch_length'}
    for u, v, key, data in graph_copy._directed.edges(keys=True, data=True):
        for attr_key in list(data.keys()):
            if attr_key not in allowed_edge_attrs:
                del graph_copy._directed[u][v][key][attr_key]
    for u, v, key, data in graph_copy._undirected.edges(keys=True, data=True):
        for attr_key in list(data.keys()):
            if attr_key not in allowed_edge_attrs:
                del graph_copy._undirected[u][v][key][attr_key]
    
    # Step 3: Orient the graph away from root vertex
    try:
        oriented_dm = orient_away_from_vertex(graph_copy, root_vertex)
    except PhyloZooError as e:
        raise PhyloZooValueError(
            f"Failed to orient network away from root location {root_location}: {e}"
        )
    
    # Step 4: Convert to DirectedPhyNetwork
    try:
        return dnetwork_from_graph(oriented_dm)
    except PhyloZooError as e:
        raise PhyloZooValueError(
            f"Failed to convert oriented network to DirectedPhyNetwork: {e}"
        )


def to_d_network(
    network: SemiDirectedPhyNetwork,
    root_location: RootLocation | None = None,
) -> 'DirectedPhyNetwork':
    """
    Convert a semi-directed network to a directed network by rooting it.
    
    The network is rooted at the specified location. If no location is provided,
    a default location is chosen from the valid root locations.
    
    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network to convert.
    root_location : RootLocation, optional
        The root location. Can be:
        - A node (T): non-leaf node in the source component
        - An edge (tuple[T, T, int]): edge in the source component as (u, v, key)
        If None, a default location is chosen from valid root locations.
        By default None.
    
    Returns
    -------
    DirectedPhyNetwork
        A directed phylogenetic network rooted at the specified location.
    
    Raises
    ------
    PhyloZooValueError
        If the root location is invalid or not in the network's valid root locations.
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 4)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
    ... )
    >>> d_net = to_d_network(net, root_location=3)
    >>> d_net.root_node
    3
    """
    # If no root location provided, find one using root_locations
    if root_location is None:
        node_locs, undir_edge_locs, dir_edge_locs = root_locations(network)
        all_valid_locations: list[RootLocation] = (
            list(node_locs) + list(undir_edge_locs) + list(dir_edge_locs)
        )
        if not all_valid_locations:
            raise PhyloZooValueError("No valid root locations found for the network")
        root_location = all_valid_locations[0]
    
    # Call the helper function to perform the rooting
    return _root_sd_network_at(network, root_location)


def root_at_outgroup(
    network: SemiDirectedPhyNetwork,
    outgroup: str,
) -> 'DirectedPhyNetwork':
    """
    Root a semi-directed network at the edge leading to the specified outgroup taxon.

    This is a convenience function that finds the edge incident to the leaf node with
    the given taxon label, then roots the network at that edge using `to_d_network`.

    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network to root.
    outgroup : str
        The taxon label of the outgroup leaf.

    Returns
    -------
    DirectedPhyNetwork
        A directed phylogenetic network rooted at the edge leading to the outgroup.

    Raises
    ------
    PhyloZooValueError
        If the outgroup taxon is not found in the network.
        If no edge is found incident to the outgroup leaf.
        If the edge is not a valid root location.

    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> from phylozoo.core.network.sdnetwork.derivations import root_at_outgroup
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 4)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
    ... )
    >>> d_net = root_at_outgroup(net, 'A')
    >>> d_net.root_node is not None
    True
    """
    # Get the node ID for the outgroup taxon
    outgroup_node = network.get_node_id(outgroup)
    if outgroup_node is None:
        raise PhyloZooValueError(f"Outgroup taxon '{outgroup}' not found in network")

    # Find an edge incident to the outgroup leaf node
    # Check undirected edges first (most common case)
    incident_edges: list[tuple[Any, Any, int]] = []

    for u, v, key in network._graph.incident_undirected_edges(outgroup_node, keys=True):
        # The edge is (u, v, key) where one of u or v is outgroup_node
        # We want the edge as (parent, outgroup_node, key) for rooting
        if u == outgroup_node:
            incident_edges.append((v, u, key))
        else:
            incident_edges.append((u, v, key))

    if not incident_edges:
        raise PhyloZooValueError(
            f"No edge found incident to outgroup leaf node {outgroup_node} "
            f"(taxon '{outgroup}')"
        )

    # Use the first incident edge as the root location
    # For undirected edges, we can root on the edge
    root_location = incident_edges[0]

    # Convert to directed network
    return to_d_network(network, root_location=root_location)


def partition_from_blob(
    network: MixedPhyNetwork,
    blob: set[Any],
    return_edge_taxa: bool = False,
) -> Partition | tuple[Partition, list[tuple[Any, Any, frozenset[str]]]]:
    """
    Get the partition of taxa induced by removing a blob from the network.
    
    When all nodes in the blob are removed, the network splits into connected components.
    Each component's taxa form a part of the partition.
    
    Parameters
    ----------
    network : MixedPhyNetwork
        The mixed phylogenetic network.
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
        If removing the blob does not disconnect the network (blob is not a cut-blob).
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 4)],
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
    
    # Find connected components
    components = list(nx.connected_components(graph_copy._combined))
    
    # Check that removing blob disconnects the network (at least 2 components)
    if len(components) < 2:
        raise PhyloZooValueError(
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
                raise PhyloZooValueError(
                    f"Could not find edge connecting component with taxa {component_taxa} to blob"
                )
    
    # Create partition (no need to check coverage - non-leaf blobs guarantee all taxa are covered)
    partition = Partition(partition_parts)
    
    if return_edge_taxa:
        return (partition, edge_taxa_list)
    return partition
