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
import itertools
from typing import Any, Iterator

from . import DirectedPhyNetwork
from .classifications import is_lsa_network
from .features import blobs
from .io import dnetwork_from_dmgraph
from .transformations import (
    to_lsa_network,
    suppress_2_blobs as suppress_2_blobs_fn,
    identify_parallel_edges as identify_parallel_edges_dn,
)
from ._utils import _suppress_deg2_nodes as dm_suppress_deg2_nodes
from ..sdnetwork._utils import _suppress_deg2_nodes as mm_suppress_deg2_nodes
from ...primitives.d_multigraph.transformations import identify_vertices as dm_identify_vertices
from ...primitives.d_multigraph.transformations import subgraph as dm_subgraph
from ...primitives.d_multigraph import DirectedMultiGraph
from ...primitives.m_multigraph import MixedMultiGraph
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
    return dnetwork_from_dmgraph(working_graph)


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

    Other parameters (keyword-only)
    -------------------------------
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
    ValueError
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
            raise ValueError(f"Taxon label '{t}' not found in network")
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
    result_net = dnetwork_from_dmgraph(working_dm)

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

    Other parameters (keyword-only)
    -------------------------------
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
    ValueError
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
            switching_graph._graph.graph['probability'] = 1.0
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
            switching_graph._graph.graph['probability'] = prob
            switching_graph._combined.graph['probability'] = prob
        
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
        # Note: dnetwork_from_dmgraph already copies graph attributes, so probability
        # is automatically preserved from the switching graph.
        displayed_tree = dnetwork_from_dmgraph(tree_graph)
        
        yield displayed_tree
