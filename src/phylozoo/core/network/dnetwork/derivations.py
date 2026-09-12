"""
Network derivations module.

This module provides functions to derive other data structures from directed
phylogenetic networks (e.g., splits, quartets, distances, blobtrees, subnetworks, etc.).
"""

import itertools
from collections import deque
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
from ...triplet import Triplet, TripletProfile, TripletProfileSet
from ...primitives.partition import Partition
from ._utils import _prune_degree1_nodes, _suppress_deg2_nodes as dm_suppress_deg2_nodes
from ..sdnetwork._utils import _suppress_deg2_nodes as mm_suppress_deg2_nodes
from ...primitives.d_multigraph.transformations import identify_vertices as dm_identify_vertices
from ...primitives.d_multigraph.transformations import subgraph as dm_subgraph
from ...primitives.d_multigraph import DirectedMultiGraph
from ...primitives.m_multigraph import MixedMultiGraph
from ..sdnetwork import SemiDirectedPhyNetwork
from ..sdnetwork.conversions import sdnetwork_from_graph
from ....core.distance import DistanceMatrix
from ....utils.exceptions import PhyloZooValueError, PhyloZooAlgorithmError
from ....utils.validation import no_validation


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
    ...     , (1, 100), (2, 101)],
    ...     nodes=[(3, {'label': 'C'}), (100, {'label': 'A'}), (101, {'label': 'B'})]
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
        nodes_list = [(leaf, {"label": label})] if label is not None else None
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
        nodes_list = [(node, {"label": label})] if label else None
        return SemiDirectedPhyNetwork(directed_edges=[], undirected_edges=[], nodes=nodes_list)

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
    mixed: Any = MixedMultiGraph(directed_edges=directed_edges, undirected_edges=undirected_edges)

    # Suppress all degree-2 nodes using the mixed graph utility function
    mm_suppress_deg2_nodes(mixed)

    # Ensure node labels are preserved in the mixed graph
    for node in mixed.nodes():
        if node in working._node_to_label:
            label = working._node_to_label[node]
            mixed._undirected.nodes[node]["label"] = label
            mixed._directed.nodes[node]["label"] = label
            mixed._combined_cache = None

    # Convert the mixed graph to a semi-directed network
    # Valid by construction from a valid input network, so validation is not re-run.
    with no_validation():
        return sdnetwork_from_graph(mixed, network_type="semi-directed", copy=False)


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
    all_blobs: Any = blobs(blob_network, trivial=False, leaves=False)

    # Work on a copy of the internal graph and collapse blobs
    working_graph = blob_network._graph.copy()
    for blob in all_blobs:
        if len(blob) > 1:
            # Sort to ensure deterministic behavior
            blob_sorted = sorted(blob)
            # Identify all vertices in the blob with the first vertex
            dm_identify_vertices(working_graph, blob_sorted)

    # Convert back to DirectedPhyNetwork
    # Valid by construction from a valid input network, so validation is not re-run.
    with no_validation():
        return dnetwork_from_graph(working_graph, copy=False)


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

    # Collect all ancestors (and the leaves themselves). The union of the leaves'
    # ancestor sets is just the set of nodes that reach any of them, so one traversal
    # of the reversed graph suffices -- calling nx.ancestors once per leaf re-walks
    # the shared upper part of the network once for every leaf.
    dag = working_net._graph._graph
    nodes_set: set[Any] = set(leaf_nodes)
    stack: list[Any] = list(leaf_nodes)
    while stack:
        node = stack.pop()
        for parent in dag.predecessors(node):
            if parent not in nodes_set:
                nodes_set.add(parent)
                stack.append(parent)

    # Create induced DirectedMultiGraph using existing utility
    induced_dm = dm_subgraph(working_net._graph, nodes_set)

    # dm_subgraph returns a fresh graph, so it can be reshaped in place.
    working_dm = induced_dm

    # First pass: suppress all degree-2 nodes (directed suppression semantics)
    dm_suppress_deg2_nodes(working_dm, exclude_nodes=None)

    # Convert to DirectedPhyNetwork for higher-level transformations
    # Valid by construction from a valid input network, so validation is not re-run.
    with no_validation():
        result_net = dnetwork_from_graph(working_dm, copy=False)

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
    >>> len(subnetworks)  # C(3,2) = 3 combinations
    3
    >>> # Each subnetwork has exactly 2 leaves
    >>> all(len(subnet.leaves) == 2 for subnet in subnetworks)
    True
    >>> # Generate all 1-taxon subnetworks
    >>> single_taxon_subs = list(k_taxon_subnetworks(net, k=1))
    >>> len(single_taxon_subs)  # C(3,1) = 3 combinations
    3
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


def _switchings(
    network: DirectedPhyNetwork, probability: bool = False
) -> Iterator[DirectedMultiGraph]:
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
    >>> len(switchings)  # Two parent edges for hybrid node 4: (5,4) and (6,4)
    2
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
    >>> switchings_with_prob[0]._graph.graph.get('probability')  # Probability of keeping edge (5,4)
    0.6
    >>> switchings_with_prob[1]._graph.graph.get('probability')  # Probability of keeping edge (6,4)
    0.4
    """
    hybrid_nodes = network.hybrid_nodes

    # If no hybrid nodes, return the original graph
    if not hybrid_nodes:
        switching_graph = network._graph.copy()
        if probability:
            switching_graph.set_graph_attribute("probability", 1.0)
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
            switching_graph.set_graph_attribute("probability", prob)

        yield switching_graph


def _displayed_tree_graphs(
    network: DirectedPhyNetwork, probability: bool = False
) -> Iterator[DirectedMultiGraph]:
    """
    Yield the underlying graph of each displayed tree.

    Shared core of :func:`displayed_trees`: applies a switching, prunes degree-1
    nodes that are neither leaves nor the root, and suppresses degree-2 nodes.
    Callers that only read a displayed tree's topology can use this and skip
    building a network object per tree.

    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network.
    probability : bool, optional
        If True, each yielded graph carries the switching's probability in its
        ``probability`` graph attribute. By default False.

    Yields
    ------
    DirectedMultiGraph
        The graph of one displayed tree. Each is a fresh object owned by the caller.
    """
    original_leaves = network.leaves
    original_root = network.root_node
    keep_nodes = set(original_leaves) | {original_root}

    for tree_graph in _switchings(network, probability=probability):
        # _switchings yields a fresh graph per switching and nothing else sees it,
        # so it can be reshaped in place rather than copied again.

        # Exhaustively remove degree-1 nodes that are not leaves or root
        _prune_degree1_nodes(tree_graph, keep_nodes)

        # Suppress all degree-2 nodes
        dm_suppress_deg2_nodes(tree_graph, exclude_nodes=None)

        yield tree_graph


def displayed_trees(
    network: DirectedPhyNetwork,
    probability: bool = False,
    make_lsa: bool = False,
) -> Iterator[DirectedPhyNetwork]:
    """
    Generate all displayed trees of a directed phylogenetic network.

    A displayed tree is obtained by:

    1. Taking a switching (deleting all but one parent edge per hybrid node)
    2. Exhaustively removing degree-1 nodes that are not leaves or the root
    3. Suppressing all degree-2 nodes

    Note: when the original network root is preserved even though its entire
    subtree on one side is pruned, the resulting tree can have a root of
    out-degree 1. Pass ``make_lsa=True`` to strip that
    redundant root and re-root at the LSA.

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
    make_lsa : bool, optional
        If True, convert each displayed tree to its LSA-network before
        yielding, removing any unary root that may arise when one side of the
        root is entirely pruned during a switching. For trees whose root
        already has out-degree >= 2 the conversion is a no-op. By default
        False.

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
    >>> len(trees)  # Two switchings yield two displayed trees
    2
    """
    for tree_graph in _displayed_tree_graphs(network, probability=probability):

        # Convert back to DirectedPhyNetwork
        # Note: dnetwork_from_graph already copies graph attributes, so probability
        # is automatically preserved from the switching graph.
        # Valid by construction from a valid input network, so validation is not re-run.
        with no_validation():
            displayed_tree = dnetwork_from_graph(tree_graph, copy=False)

        if make_lsa:
            displayed_tree = to_lsa_network(displayed_tree)

        yield displayed_tree


def _hybrid_parent_edges(network: DirectedPhyNetwork) -> dict[Any, list[tuple[Any, Any, int]]]:
    """Map each hybrid node to all of its parent edges as ``(u, v, key)``."""
    return {
        hybrid: list(network.incident_parent_edges(hybrid, keys=True))
        for hybrid in network.hybrid_nodes
    }


def _branch_lengths(network: DirectedPhyNetwork) -> dict[tuple[Any, Any, int], float]:
    """
    Branch length of every edge, keyed ``(u, v, key)`` in both orientations.

    Looked up once per network rather than once per switching; edges without a
    recorded length count as 1.0, as everywhere in the switching-based derivations.
    """
    lengths: dict[tuple[Any, Any, int], float] = {}
    for u, v, key in network._graph._combined.edges(keys=True):
        length = network.get_branch_length(u, v, key)
        if length is None:
            length = network.get_branch_length(v, u, key)
        if length is None:
            length = 1.0
        lengths[(u, v, key)] = length
        lengths[(v, u, key)] = length
    return lengths


def _switching_adjacency(
    network: DirectedPhyNetwork, removed_edges: set[tuple[Any, Any, int]]
) -> dict[Any, list[tuple[Any, int]]]:
    """
    Undirected adjacency of one switching, without building the switching graph.

    A switching is the network minus the parent edges it drops, so its adjacency is
    the network's cached combined graph with ``removed_edges`` skipped. It is
    materialised once per switching (O(V + E)) because the distance computation
    walks it once per leaf; every switching-based traversal (distances, displayed
    splits) uses this instead of copying the graph per switching.

    Parameters
    ----------
    network : DirectedPhyNetwork
        The network.
    removed_edges : set[tuple[Any, Any, int]]
        The hybrid parent edges ``(u, v, key)`` this switching deletes.

    Returns
    -------
    dict[Any, list[tuple[Any, int]]]
        For every node, its ``(neighbour, key)`` pairs over the surviving edges.
    """
    combined = network._graph._combined
    adjacency: dict[Any, list[tuple[Any, int]]] = {node: [] for node in combined.nodes()}
    for u, v, key in combined.edges(keys=True):
        if (u, v, key) in removed_edges or (v, u, key) in removed_edges:
            continue
        adjacency[u].append((v, key))
        if u != v:
            adjacency[v].append((u, key))
    return adjacency


class _BlobSwitchings:
    """
    The covering set of switchings that blob-wise aggregation evaluates.

    A switching keeps one parent edge per hybrid node. Rather than enumerating the
    ``prod_b 2**r_b`` switchings of the whole network, the derivations built on this
    class evaluate one arbitrary *reference* switching plus, for each blob, every
    switching that differs from the reference inside that blob alone -- ``1 + sum_b
    2**r_b`` in total, where ``r_b`` is the number of reticulations in blob ``b``.

    That set suffices because anything a leaf-to-leaf path does inside a blob (its
    length, the splits its edges induce) changes only when that blob's own hybrid
    choices change, and blobs are switched independently. Each evaluated switching
    therefore exposes one blob's contribution in isolation, and the aggregate over all
    switchings follows by combining the per-blob contributions. The work is exponential
    in the network's level, not in its total number of reticulations; for a single-blob
    network the covering set is every switching and nothing is saved.

    Attributes
    ----------
    hybrid_parent_edges : dict
        Each hybrid node's parent edges.
    reference : dict
        The reference switching: the first parent edge of every hybrid.
    groups : list[list[Any]]
        Hybrid nodes grouped by blob (see :func:`_hybrid_blob_groups`).
    masses : list[float]
        Per group, the total gamma weight of its local switchings.
    """

    def __init__(self, network: DirectedPhyNetwork) -> None:
        self.network = network
        self.hybrid_parent_edges = _hybrid_parent_edges(network)
        self.reference = {hybrid: edges[0] for hybrid, edges in self.hybrid_parent_edges.items()}
        self.groups = _hybrid_blob_groups(network)
        self.masses = [
            sum(weight for _choices, weight in self.local_switchings(index))
            for index in range(len(self.groups))
        ]

    def weight(self, group: list[Any], combination: tuple[tuple[Any, Any, int], ...]) -> float:
        """Gamma weight of keeping ``combination`` at the hybrids of ``group``."""
        weight = 1.0
        for hybrid, edge in zip(group, combination):
            weight *= _hybrid_edge_weight(self.network, edge, len(self.hybrid_parent_edges[hybrid]))
        return weight

    def local_switchings(
        self, index: int
    ) -> Iterator[tuple[dict[Any, tuple[Any, Any, int]], float]]:
        """Yield ``(choices, weight)`` for every switching of group ``index`` alone."""
        group = self.groups[index]
        for combination in itertools.product(
            *(self.hybrid_parent_edges[hybrid] for hybrid in group)
        ):
            choices = dict(self.reference)
            choices.update(zip(group, combination))
            yield choices, self.weight(group, combination)

    def removed_edges(self, choices: dict[Any, tuple[Any, Any, int]]) -> set[tuple[Any, Any, int]]:
        """The parent edges the switching ``choices`` deletes."""
        return {
            edge
            for hybrid, edges in self.hybrid_parent_edges.items()
            for edge in edges
            if edge != choices[hybrid]
        }

    @property
    def total_mass(self) -> float:
        """Total gamma weight of all switchings: the product of the group masses."""
        return self.other_mass(-1)

    def other_mass(self, index: int) -> float:
        """Product of the masses of every group except ``index``."""
        mass = 1.0
        for other, group_mass in enumerate(self.masses):
            if other != index:
                mass *= group_mass
        return mass


def _hybrid_blob_groups(network: DirectedPhyNetwork) -> list[list[Any]]:
    """
    Group the hybrid nodes of a network by the blob that contains them.

    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network.

    Returns
    -------
    list[list[Any]]
        One list of hybrid nodes per blob that contains at least one hybrid.

    Notes
    -----
    This is what makes the blob decomposition in :func:`distances` possible. The
    stretch of a path that lies inside a blob depends only on that blob's hybrid
    choices, and different blobs are switched independently, so the groups returned
    here can be enumerated separately instead of jointly.

    Correctness does not require the grouping to be maximal -- only that hybrids in
    different groups never interact. A finer grouping stays correct and merely
    enumerates more combinations than necessary.
    """
    hybrid_nodes = list(network.hybrid_nodes)
    if not hybrid_nodes:
        return []

    remaining = set(hybrid_nodes)
    groups: list[list[Any]] = []
    for blob in blobs(network, trivial=False, leaves=False):
        members = [hybrid for hybrid in hybrid_nodes if hybrid in blob and hybrid in remaining]
        if members:
            groups.append(members)
            remaining.difference_update(members)

    # A hybrid no blob claimed gets a group of its own, which is always safe.
    groups.extend([hybrid] for hybrid in hybrid_nodes if hybrid in remaining)
    return groups


def _switching_splits_by_blob(
    network: DirectedPhyNetwork,
    removed_edges: set[tuple[Any, Any, int]],
    node_to_blob: dict[Any, int],
    all_taxa: frozenset[str],
    only_blob: int | None = None,
) -> tuple[set[Split], dict[int, set[Split]]]:
    """
    Splits induced by the edges of one switching, grouped by the blob owning each edge.

    The switching is the network minus ``removed_edges`` (the parent edges it drops),
    which is a tree, so every remaining edge induces the bipartition of the taxa
    below it against the rest. Edges with no taxa on one side -- the ones a displayed
    tree prunes away -- induce nothing. An edge whose endpoints lie in different blobs
    is a cut edge of the network; otherwise it lies inside the blob both endpoints
    share.

    The traversal uses :func:`_switching_adjacency`, so no switching graph is built.

    Parameters
    ----------
    network : DirectedPhyNetwork
        The network.
    removed_edges : set[tuple[Any, Any, int]]
        The hybrid parent edges ``(u, v, key)`` that this switching deletes.
    node_to_blob : dict
        Maps every node to the index of its blob; blobs partition the nodes.
    all_taxa : frozenset[str]
        All taxon labels of the network.
    only_blob : int | None, optional
        If given, build splits only for the edges inside this blob; every other edge
        is skipped (its taxa are still accumulated). Building a split costs O(n), so
        this is what keeps a blob's local switchings cheap. By default None.

    Returns
    -------
    tuple[set[Split], dict[int, set[Split]]]
        The cut-edge splits (empty when ``only_blob`` is set), and the splits of the
        edges inside each blob considered.
    """
    adjacency = _switching_adjacency(network, removed_edges)
    label_of = network._node_to_label
    leaf_nodes = set(network.leaves)
    total = len(all_taxa)

    # Depth-first order from a leaf; reversing it visits every node after all of its
    # descendants, so taxa sets can be accumulated upwards in one pass.
    start = next(iter(leaf_nodes))
    parent: dict[Any, Any] = {start: None}
    order: list[Any] = []
    stack = [start]
    while stack:
        node = stack.pop()
        order.append(node)
        for neighbour, _key in adjacency[node]:
            if neighbour not in parent:
                parent[neighbour] = node
                stack.append(neighbour)

    cut_splits: set[Split] = set()
    blob_splits: dict[int, set[Split]] = {}
    below: dict[Any, set[str]] = {}
    for node in reversed(order):
        # A node's set is needed only until it is merged into its parent's; dropping
        # it here keeps memory at the pending frontier rather than O(n * V).
        taxa_below = below.pop(node, set())
        if node in leaf_nodes:
            label = label_of.get(node)
            if label is not None:
                taxa_below.add(label)
        up = parent[node]
        if up is None:
            continue
        if taxa_below and len(taxa_below) < total:
            blob = node_to_blob[node]
            inside = node_to_blob[up] == blob
            if only_blob is None or (inside and blob == only_blob):
                side = frozenset(taxa_below)
                split = Split(side, all_taxa - side)
                if inside:
                    blob_splits.setdefault(blob, set()).add(split)
                else:
                    cut_splits.add(split)
        below.setdefault(up, set()).update(taxa_below)
    return cut_splits, blob_splits


def _hybrid_edge_weight(
    network: DirectedPhyNetwork,
    edge: tuple[Any, Any, int],
    indegree: int,
) -> float:
    """
    Probability of keeping a hybrid's parent ``edge``: its gamma, or 1/indegree if unset.
    """
    gamma = network.get_gamma(edge[0], edge[1], edge[2])
    return float(gamma) if gamma is not None else 1.0 / indegree


def _switching_distance_matrix(
    network: DirectedPhyNetwork,
    all_taxa: list[str],
    plan: _BlobSwitchings,
    choices: dict[Any, tuple[Any, Any, int]],
    branch_lengths: dict[tuple[Any, Any, int], float],
) -> np.ndarray:
    """
    Distance matrix of the switching that keeps exactly the parent edges in ``choices``.

    A switching is a tree, so one traversal from each leaf yields its distance to
    every other leaf. The traversal uses :func:`_switching_adjacency`, so no
    switching graph is built.

    Parameters
    ----------
    network : DirectedPhyNetwork
        The network.
    all_taxa : list[str]
        Taxon labels, in the order used for the matrix rows and columns.
    plan : _BlobSwitchings
        The switching machinery of ``network``.
    choices : dict
        Maps each hybrid node to the single parent edge to keep.
    branch_lengths : dict
        Output of :func:`_branch_lengths` for ``network``.

    Returns
    -------
    numpy.ndarray
        The pairwise distance matrix for that switching.
    """
    n = len(all_taxa)
    leaf_nodes: list[Any] = []
    for taxon in all_taxa:
        leaf_node = network._label_to_node.get(taxon)
        if leaf_node is None:
            raise PhyloZooValueError(f"Taxon '{taxon}' not found in network")
        leaf_nodes.append(leaf_node)
    leaf_index = {leaf: i for i, leaf in enumerate(leaf_nodes)}
    adjacency = _switching_adjacency(network, plan.removed_edges(choices))

    distance_matrix = np.zeros((n, n), dtype=np.float64)
    for i, source in enumerate(leaf_nodes):
        visited = {source}
        queue: deque[tuple[Any, float]] = deque([(source, 0.0)])
        while queue:
            node, distance = queue.popleft()
            target = leaf_index.get(node)
            if target is not None:
                distance_matrix[i, target] = distance
            for neighbour, key in adjacency[node]:
                if neighbour not in visited:
                    visited.add(neighbour)
                    queue.append((neighbour, distance + branch_lengths[(node, neighbour, key)]))
    return distance_matrix


def distances(
    network: DirectedPhyNetwork,
    mode: Literal["shortest", "longest", "average"] = "average",
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
    >>> dm.get_distance('A', 'B')  # Example distance
    2.0
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

    if mode not in ("shortest", "longest", "average"):
        raise PhyloZooValueError(
            f"Invalid mode: {mode}. Must be 'shortest', 'longest', or 'average'"
        )

    # Only the covering set of switchings in _BlobSwitchings is evaluated: a reference
    # plus each blob's local switchings. Any switching's distance matrix equals the
    # reference matrix plus one independent difference per blob, so each mode
    # aggregates one blob's differences at a time -- elementwise minimum, maximum or
    # gamma-weighted mean -- and sums the per-blob results onto the reference. Minimum
    # and maximum are valid here for the same reason the mean is: blobs are switched
    # independently, and a sum of independently chosen terms is optimised term by term.
    plan = _BlobSwitchings(network)
    branch_lengths = _branch_lengths(network)
    # ``base`` stays fixed: every per-blob difference is measured against this one
    # reference switching, while ``result`` accumulates them.
    base = _switching_distance_matrix(network, all_taxa, plan, plan.reference, branch_lengths)
    result = base.copy()

    for index in range(len(plan.groups)):
        # Accumulated per element, never stacked: a blob with r reticulations has 2**r
        # local switchings, and holding that many n x n matrices at once is not viable.
        aggregate: np.ndarray | None = None
        weight_sum = 0.0
        for choices, weight in plan.local_switchings(index):
            delta = (
                _switching_distance_matrix(network, all_taxa, plan, choices, branch_lengths) - base
            )
            if mode == "shortest":
                aggregate = delta if aggregate is None else np.minimum(aggregate, delta)
            elif mode == "longest":
                aggregate = delta if aggregate is None else np.maximum(aggregate, delta)
            else:
                weight_sum += weight
                delta *= weight
                aggregate = delta if aggregate is None else aggregate + delta
        if aggregate is None:
            continue
        if mode == "average":
            # Normalising per blob matches normalising globally: the total weight of
            # all switchings factorises over the blobs.
            if weight_sum > 0.0:
                aggregate = aggregate / weight_sum
        result = result + aggregate

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

    # Built once: the complement below is taken for every edge, and rebuilding this
    # set per edge is what made the traversal quadratic in the taxon count.
    all_taxa_set = frozenset(all_taxa)
    total_taxa = len(all_taxa_set)

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
                if neighbor_leaves and len(neighbor_leaves) < total_taxa:
                    other_leaves = all_taxa_set - neighbor_leaves
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
            f"Multiple parallel edges exist between {u} and {v}. " "Must specify 'key' parameter."
        )
    if key is None:
        key = edge_keys[0]
    if key not in edge_keys:
        raise PhyloZooValueError(f"Edge ({u}, {v}, key={key}) does not exist")

    # Remove the edge
    graph_copy.remove_edge(u, v, key)

    # Check if removal disconnects the graph
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

    split: Any = Split(taxa_u, taxa_v)

    if return_node_taxa:
        return (split, (u, taxa_u), (v, taxa_v))
    return split  # type: ignore[no-any-return]


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
    all_taxa = frozenset(network.taxa)
    if len(all_taxa) < 2:
        return WeightedSplitSystem()

    # Only the covering set of switchings in _BlobSwitchings is evaluated. A displayed
    # tree's splits are induced by its edges, and each edge is either a cut edge of the
    # network -- in every displayed tree, always the same split, so it carries the whole
    # gamma mass -- or lies inside one blob, so its split depends only on that blob's
    # switching and its weight is a sum over that blob's local switchings alone, times
    # the mass of every other blob.
    #
    # Two details keep the weights identical to summing over displayed trees. A path
    # of degree-2 nodes created by a switching can make a blob edge repeat an adjacent
    # cut edge's split, and a displayed tree counts each split once, so blob splits
    # equal to a cut-edge split are dropped. And the masses are the raw gamma
    # products, not normalised, exactly as the displayed-tree probabilities were.
    node_to_blob: dict[Any, int] = {}
    for index, blob in enumerate(blobs(network, trivial=True, leaves=True)):
        for node in blob:
            node_to_blob[node] = index

    plan = _BlobSwitchings(network)
    split_weights: dict[Split, float] = {}
    cut_splits, _ = _switching_splits_by_blob(
        network, plan.removed_edges(plan.reference), node_to_blob, all_taxa
    )
    total_mass = plan.total_mass
    for split in cut_splits:
        split_weights[split] = total_mass

    for index, group in enumerate(plan.groups):
        blob_id = node_to_blob[group[0]]
        other_mass = plan.other_mass(index)
        for choices, weight in plan.local_switchings(index):
            _, blob_splits = _switching_splits_by_blob(
                network, plan.removed_edges(choices), node_to_blob, all_taxa, only_blob=blob_id
            )
            for split in blob_splits.get(blob_id, ()):
                if split in cut_splits:
                    continue
                split_weights[split] = split_weights.get(split, 0.0) + weight * other_mass

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


def displayed_triplets(network: DirectedPhyNetwork) -> TripletProfileSet:
    """
    Compute triplet profile set from all displayed trees of the network.

    For each triplet (3-leaf subnetwork), this function:
    1. Extracts the subnetwork induced by those 3 taxa
    2. Gets all displayed trees of that subnetwork (with probabilities)
    3. Converts each displayed tree to a rooted triplet by identifying the outgroup (the leaf that is a direct child of the root) and the cherry (the two remaining leaves sharing an internal parent below the root); if all three leaves are direct children of the root, the triplet is unresolved (star)
    4. Creates a triplet profile where each triplet's weight is the probability of the displayed tree that induced it (summing weights if the same triplet appears in multiple displayed trees)

    The profiles are then returned as a TripletProfileSet, where each profile
    (one per 3-taxon set) has triplet weights that sum to 1.0 (the displayed-tree
    probabilities). Each profile in the set has default profile weight 1.0.

    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network.

    Returns
    -------
    TripletProfileSet
        A triplet profile set where each profile corresponds to a 3-taxon set,
        and contains triplets from displayed trees weighted by their probabilities.

    Examples
    --------
    >>> net = DirectedPhyNetwork(
    ...     edges=[(4, 1), (4, 5), (5, 2), (5, 3)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
    ... )
    >>> profileset = displayed_triplets(net)
    >>> isinstance(profileset, TripletProfileSet)
    True
    >>> len(profileset) > 0
    True
    """
    # Handle networks with fewer than 3 taxa
    taxa_list = sorted(network.taxa)
    if len(taxa_list) < 3:
        return TripletProfileSet()

    # Collect profiles for each 3-taxon set
    profiles: list[TripletProfile] = []

    # Iterate through all combinations of 3 taxa
    for three_taxa in itertools.combinations(taxa_list, 3):
        three_taxa_set = frozenset(three_taxa)

        # Get subnetwork induced by these 3 taxa
        triplet_subnet = subnetwork(network, list(three_taxa))

        # Collect triplets with their weights for this 3-taxon set
        triplet_weights: dict[Triplet, float] = {}

        # Get all displayed trees of the subnetwork with probabilities
        for tree_graph in _displayed_tree_graphs(triplet_subnet, probability=True):
            # Read the topology straight off the graph: a rooted three-leaf tree needs
            # no network object, only its root and which leaves hang directly off it.
            graph = tree_graph._graph
            prob = graph.graph.get("probability")
            if prob is None:
                prob = 1.0

            # Find the effective root by walking down past any single-child stem.
            # The displayed tree may have a single-child path from the topological
            # root down to the first branching node; the rooted topology of the
            # triplet is determined by what hangs below that branching node.
            root = next(node for node in graph.nodes() if graph.in_degree(node) == 0)
            while graph.out_degree(root) == 1:
                root = next(iter(graph.successors(root)))

            # A star has all 3 leaves as direct root-children, and a binary
            # triplet has exactly one leaf as a direct root-child (the outgroup)
            # and the other two as siblings under an internal node.
            leaf_nodes = [node for node in graph.nodes() if graph.out_degree(node) == 0]
            root_children = set(graph.successors(root))
            direct_leaves = [leaf for leaf in leaf_nodes if leaf in root_children]

            def _label(node: Any, _graph: Any = graph) -> str:
                label = _graph.nodes[node].get("label")
                return str(label) if label is not None else str(node)

            if len(direct_leaves) == 3:
                # Star triplet: all 3 leaves are root children
                triplet = Triplet(three_taxa_set)
            elif len(direct_leaves) == 1:
                # Binary triplet: the direct leaf is the outgroup, others form the cherry
                outgroup_label = _label(direct_leaves[0])
                cherry_labels = {_label(leaf) for leaf in leaf_nodes if leaf != direct_leaves[0]}
                triplet = Triplet(Split({outgroup_label}, cherry_labels))
            else:
                # Unexpected structure for a rooted 3-leaf tree
                raise PhyloZooAlgorithmError(
                    f"Displayed tree on 3 taxa has {len(direct_leaves)} leaves as "
                    "direct root children; expected 1 (binary) or 3 (star)."
                )

            # Add triplet to profile with its weight (sum if duplicate)
            triplet_weights[triplet] = triplet_weights.get(triplet, 0.0) + prob

        # Create profile for this 3-taxon set (only if we have triplets)
        if triplet_weights:
            profile = TripletProfile(triplet_weights)
            profiles.append(profile)

    # Create TripletProfileSet (each profile gets default weight 1.0)
    return TripletProfileSet(profiles=profiles)


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
        raise PhyloZooValueError(f"Blob contains nodes not in network: {missing_nodes}")

    # Check that blob is a non-leaf blob (including trivial single-node blobs)
    non_leaf_blobs: Any = blobs(network, trivial=True, leaves=False)
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
