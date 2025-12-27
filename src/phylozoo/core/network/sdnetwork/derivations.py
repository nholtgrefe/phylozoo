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
from typing import Any, Iterator

from . import MixedPhyNetwork
from .features import blobs
from .transformations import suppress_2_blobs as suppress_2_blobs_fn, identify_parallel_edges as identify_parallel_edges_fn
from .sd_phynetwork import SemiDirectedPhyNetwork
from ._utils import _suppress_deg2_nodes
from .io import sdnetwork_from_mmgraph
from ...primitives.m_multigraph.transformations import (
    identify_vertices as mm_identify_vertices,
    suppress_degree2_node as mm_suppress_degree2_node,
    subgraph as mm_subgraph,
)
from ...primitives.m_multigraph.features import updown_path_vertices
from ...primitives.m_multigraph import MixedMultiGraph
from .io import sdnetwork_from_mmgraph


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
    MixedPhyNetwork
        A new network where each blob has been collapsed to a single vertex,
        forming a tree structure.

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
    return sdnetwork_from_mmgraph(working_graph, network_type=network_type)


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
    result_net = sdnetwork_from_mmgraph(working_mm, network_type='semi-directed')
    
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
            switching_graph._directed.graph['probability'] = 1.0
            switching_graph._undirected.graph['probability'] = 1.0
            switching_graph._combined.graph['probability'] = 1.0
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
            switching_graph._directed.graph['probability'] = prob
            switching_graph._undirected.graph['probability'] = prob
            switching_graph._combined.graph['probability'] = prob
        
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
        # Note: sdnetwork_from_mmgraph already copies graph attributes, so probability
        # is automatically preserved from the switching graph.
        displayed_tree = sdnetwork_from_mmgraph(tree_graph, network_type='semi-directed')
        
        yield displayed_tree
