"""
Network derivations module.

This module provides functions to derive other data structures from semi-directed
and mixed phylogenetic networks (e.g., splits, quartets, distances, blobtrees, subnetworks, etc.).
"""

import itertools
from collections import deque
from typing import TYPE_CHECKING, Any, Iterator, Literal

import numpy as np
import networkx as nx

if TYPE_CHECKING:
    from ...network.dnetwork import DirectedPhyNetwork

from . import MixedPhyNetwork
from .classifications import is_tree
from .features import blobs, root_locations, RootLocation
from .transformations import (
    suppress_2_blobs as suppress_2_blobs_fn,
    _identify_parallel_edges_inplace,
)
from ...split import Split, SplitSystem, WeightedSplitSystem
from ...quartet import Quartet, QuartetProfile, QuartetProfileSet
from ...primitives.partition import Partition
from .sd_phynetwork import SemiDirectedPhyNetwork
from ._utils import _RootingContext, _prune_degree1_nodes, _suppress_deg2_nodes, _subdivide_edge
from .conversions import sdnetwork_from_graph
from ...primitives.m_multigraph.transformations import (
    identify_vertices as mm_identify_vertices,
    subgraph as mm_subgraph,
    orient_away_from_vertex,
)
from ...primitives.m_multigraph import MixedMultiGraph
from ....core.distance import DistanceMatrix
from ....utils.exceptions import PhyloZooValueError, PhyloZooError
from ....utils.validation import no_validation


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
    all_blobs: Any = blobs(blob_network, trivial=False, leaves=False)

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
    network_type = "semi-directed" if isinstance(network, SemiDirectedPhyNetwork) else "mixed"
    # Valid by construction from a valid input network, so validation is not re-run.
    with no_validation():
        return sdnetwork_from_graph(working_graph, network_type=network_type, copy=False)


def subnetwork(
    network: SemiDirectedPhyNetwork,
    taxa: list[str],
    suppress_2_blobs: bool = False,
    identify_parallel_edges: bool = False,
    _rooting: "_RootingContext | None" = None,
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

    # Collect the vertices lying on up-down paths between the requested leaves.
    context = _rooting if _rooting is not None else _RootingContext(network)
    nodes_set: set[Any] = context.node_set(leaf_nodes)

    # Create induced MixedMultiGraph — mm_subgraph always returns a fresh object,
    # so no additional copy is needed before mutation.
    working_mm = mm_subgraph(network._graph, nodes_set)

    # Suppress degree-2 nodes (excluding the target leaves)
    leaf_set = set(leaf_nodes)
    _suppress_deg2_nodes(working_mm, exclude_nodes=leaf_set)

    # Parallel-edge identification runs directly on the mutable working_mm so
    # that the subsequent sdnetwork_from_graph call is the only construction.
    if identify_parallel_edges:
        _identify_parallel_edges_inplace(working_mm, exclude_nodes=leaf_set)

    with no_validation():
        result_net = sdnetwork_from_graph(working_mm, network_type="semi-directed", copy=False)

    if suppress_2_blobs:
        result_net = suppress_2_blobs_fn(result_net)

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
    >>> len(subnetworks)  # C(4,2) = 6 combinations
    6
    >>> # Each subnetwork has exactly 2 leaves
    >>> all(len(subnet.taxa) == 2 for subnet in subnetworks)
    True
    >>> # Generate all 1-taxon subnetworks
    >>> single_taxon_subs = list(k_taxon_subnetworks(net, k=1))
    >>> len(single_taxon_subs)  # C(4,1) = 4 combinations
    4
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

    # Root the network (and build its dominator tree) once, then reuse it for every
    # subset: this is what makes the per-subset work proportional to the subset size.
    rooting = _RootingContext(network)

    # Generate all combinations of k taxa
    for taxa_combination in itertools.combinations(all_taxa, k):
        yield subnetwork(
            network,
            list(taxa_combination),
            suppress_2_blobs=suppress_2_blobs,
            identify_parallel_edges=identify_parallel_edges,
            _rooting=rooting,
        )


def _switchings(
    network: SemiDirectedPhyNetwork, probability: bool = False
) -> Iterator[MixedMultiGraph]:
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
    >>> len(switchings)  # Two parent edges for hybrid node 4: (5,4) and (6,4)
    2
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
    >>> switchings_with_prob[0]._directed.graph.get('probability')  # Probability of keeping edge (5,4)
    0.6
    >>> switchings_with_prob[1]._directed.graph.get('probability')  # Probability of keeping edge (6,4)
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


def _undirect_switching(switching_graph: MixedMultiGraph) -> None:
    """
    Undirect every edge of a switching, in place.

    A switching keeps exactly one parent edge per hybrid node, so afterwards no edge
    is a reticulation edge any more and the displayed tree it yields is undirected.
    Leaving the kept edges directed lets a node that was the tail of two hybrid edges
    retain two outgoing directed edges; degree-2 suppression cannot orient such a node
    and raises.

    ``gamma`` is dropped along with the direction: it is an inheritance probability at
    a reticulation, is meaningless once the edge is an ordinary tree edge, and is
    rejected by validation on an undirected edge.

    Parameters
    ----------
    switching_graph : MixedMultiGraph
        The switching to undirect. **Modified in place.**

    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
    >>> G = MixedMultiGraph()
    >>> _ = G.add_directed_edge(1, 2, gamma=0.7)
    >>> _ = G.add_undirected_edge(2, 3)
    >>> _undirect_switching(G)
    >>> list(G.directed_edges_iter())
    []
    >>> sorted(tuple(sorted(e)) for e in G.undirected_edges_iter())
    [(1, 2), (2, 3)]
    """
    for u, v, key, data in list(switching_graph.directed_edges_iter(keys=True, data=True)):
        attrs = {name: value for name, value in (data or {}).items() if name != "gamma"}
        switching_graph.remove_directed_edge(u, v, key=key)
        switching_graph.add_undirected_edge(u, v, key=key, **attrs)


def _displayed_tree_graphs(
    network: SemiDirectedPhyNetwork, probability: bool = False
) -> Iterator[MixedMultiGraph]:
    """
    Yield the underlying graph of each displayed tree.

    This is the shared core of :func:`displayed_trees`: it applies a switching,
    undirects it, prunes non-leaf degree-1 nodes and suppresses degree-2 nodes,
    yielding the resulting graph. Callers that only need to read a displayed tree's
    topology can use this directly and skip building a network object per tree.

    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network.
    probability : bool, optional
        If True, each yielded graph carries the switching's probability in its
        ``probability`` graph attribute. By default False.

    Yields
    ------
    MixedMultiGraph
        The graph of one displayed tree. Each is a fresh object owned by the caller.
    """
    original_leaves = network.leaves
    keep_nodes = set(original_leaves)

    for tree_graph in _switchings(network, probability=probability):
        # _switchings yields a fresh graph per switching and nothing else sees it,
        # so it can be reshaped in place rather than copied again.
        _undirect_switching(tree_graph)

        # Exhaustively remove degree-1 nodes that are not leaves
        _prune_degree1_nodes(tree_graph, keep_nodes)

        # Suppress all degree-2 nodes
        _suppress_deg2_nodes(tree_graph, exclude_nodes=None)

        yield tree_graph


def displayed_trees(
    network: SemiDirectedPhyNetwork, probability: bool = False
) -> Iterator[SemiDirectedPhyNetwork]:
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
    >>> len(trees)  # Two switchings yield two displayed trees
    2
    """
    for tree_graph in _displayed_tree_graphs(network, probability=probability):

        # Convert back to SemiDirectedPhyNetwork
        # Note: sdnetwork_from_graph already copies graph attributes, so probability
        # is automatically preserved from the switching graph.
        # Valid by construction from a valid input network, so validation is not re-run.
        with no_validation():
            displayed_tree = sdnetwork_from_graph(
                tree_graph, network_type="semi-directed", copy=False
            )

        yield displayed_tree


def _hybrid_parent_edges(network: SemiDirectedPhyNetwork) -> dict[Any, list[tuple[Any, Any, int]]]:
    """Map each hybrid node to all of its parent edges as ``(u, v, key)``."""
    return {
        hybrid: list(network.incident_parent_edges(hybrid, keys=True))
        for hybrid in network.hybrid_nodes
    }


def _branch_lengths(network: SemiDirectedPhyNetwork) -> dict[tuple[Any, Any, int], float]:
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
    network: SemiDirectedPhyNetwork, removed_edges: set[tuple[Any, Any, int]]
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
    network : SemiDirectedPhyNetwork
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

    def __init__(self, network: SemiDirectedPhyNetwork) -> None:
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


def _hybrid_blob_groups(network: SemiDirectedPhyNetwork) -> list[list[Any]]:
    """
    Group the hybrid nodes of a network by the blob that contains them.

    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network.

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
    network: SemiDirectedPhyNetwork,
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
    network : SemiDirectedPhyNetwork
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
    network: SemiDirectedPhyNetwork,
    edge: tuple[Any, Any, int],
    indegree: int,
) -> float:
    """
    Probability of keeping a hybrid's parent ``edge``: its gamma, or 1/indegree if unset.
    """
    gamma = network.get_gamma(edge[0], edge[1], edge[2])
    return float(gamma) if gamma is not None else 1.0 / indegree


def _switching_distance_matrix(
    network: SemiDirectedPhyNetwork,
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
    network : SemiDirectedPhyNetwork
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
    network: SemiDirectedPhyNetwork,
    mode: Literal["shortest", "longest", "average"] = "average",
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
    ...     undirected_edges=[(3, 1), (3, 2), (3, 100)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (100, {'label': 'C'})]
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

    # Built once: the complement below is taken for every edge, and rebuilding this
    # set per edge is what made the traversal quadratic in the taxon count.
    all_taxa_set = frozenset(all_taxa)
    total_taxa = len(all_taxa_set)

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

    split: Any = Split(taxa_u, taxa_v)

    if return_node_taxa:
        return (split, (u, taxa_u), (v, taxa_v))
    return split  # type: ignore[no-any-return]


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


def _quartet_split_from_tree_graph(tree_graph: MixedMultiGraph) -> "Split | None":
    """
    Read the 2|2 split of a four-leaf tree directly from its graph.

    A four-leaf unrooted tree is either resolved -- exactly one internal edge, whose
    removal separates the leaves 2|2 -- or a star, which has none.

    Parameters
    ----------
    tree_graph : MixedMultiGraph
        Graph of a displayed tree on four leaves.

    Returns
    -------
    Split | None
        The 2|2 split, or None when the tree is a star.
    """
    labels: dict[Any, str] = {}
    for node in tree_graph.nodes():
        if tree_graph.degree(node) == 1:
            label = tree_graph._undirected.nodes.get(node, {}).get("label")
            labels[node] = str(label) if label is not None else str(node)
    if len(labels) != 4:
        return None

    adjacency: dict[Any, list[Any]] = {node: [] for node in tree_graph.nodes()}
    for u, v, _key in tree_graph.undirected_edges_iter(keys=True):
        adjacency[u].append(v)
        adjacency[v].append(u)

    for u, v, _key in tree_graph.undirected_edges_iter(keys=True):
        if u in labels or v in labels:
            continue  # a pendant edge can only cut off a single leaf
        # Walk the component containing u without crossing the edge (u, v)
        seen = {u}
        stack = [u]
        while stack:
            node = stack.pop()
            for neighbour in adjacency[node]:
                if neighbour in seen or (node == u and neighbour == v):
                    continue
                seen.add(neighbour)
                stack.append(neighbour)
        if v in seen:
            continue  # not a cut edge (parallel edges keep both sides connected)
        side = {labels[leaf] for leaf in labels if leaf in seen}
        if len(side) == 2:
            return Split(side, set(labels.values()) - side)
    return None


def displayed_quartets(network: SemiDirectedPhyNetwork) -> QuartetProfileSet:
    """
    Compute quartet profile set from all displayed trees of the network.

    For each quartet (4-leaf subnetwork), this function:
    1. Extracts the subnetwork induced by those 4 taxa
    2. Gets all displayed trees of that subnetwork (with probabilities)
    3. Converts each displayed tree to a quartet (4-leaf tree)
    4. Creates a quartet profile where each quartet's weight is the probability of the displayed tree that induced it (summing weights if the same quartet appears in multiple displayed trees)

    The profiles are then returned as a QuartetProfileSet, where each profile
    (one per 4-taxon set) has quartet weights that sum to 1.0 (the displayed-tree
    probabilities). Each profile in the set has default profile weight 1.0.

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

    # Root the network (and build its dominator tree) once and reuse it for every
    # 4-taxon set; otherwise each subnetwork call re-roots the whole network.
    rooting = _RootingContext(network)

    # Iterate through all combinations of 4 taxa
    for four_taxa in itertools.combinations(taxa_list, 4):
        four_taxa_set = frozenset(four_taxa)

        # Get subnetwork induced by these 4 taxa
        quartet_subnet = subnetwork(network, list(four_taxa), _rooting=rooting)

        # Collect quartets with their weights for this 4-taxon set
        quartet_weights: dict[Quartet, float] = {}

        # Read each displayed tree's quartet straight off its graph: a 4-leaf tree has
        # at most one 2|2 split, so building a network object per tree (and deriving
        # its whole split system) is unnecessary here.
        for tree_graph in _displayed_tree_graphs(quartet_subnet, probability=True):
            prob = tree_graph._undirected.graph.get("probability")
            if prob is None:
                prob = 1.0

            quartet_split = _quartet_split_from_tree_graph(tree_graph)

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
    validate: bool = True,
) -> "DirectedPhyNetwork":
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
    validate : bool, optional
        Whether to validate the resulting directed network. A user-supplied root
        location may be invalid, and that is detected by this validation; a location
        taken from :func:`root_locations` yields a valid network by construction, so
        callers may pass False. By default True.

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
        source_graph, subdiv_node = _subdivide_edge(network, u, v, key)

        # Root vertex is the subdivision node
        root_vertex = subdiv_node
    else:
        # root_location is a node. No copy is needed: orient_away_from_vertex only
        # reads its input and builds a fresh DirectedMultiGraph, so the network's own
        # graph is never touched. Attributes are filtered on that result below.
        source_graph = network._graph
        if not source_graph.has_node(root_location):
            raise PhyloZooValueError(f"Node {root_location} not found in the network")
        root_vertex = root_location

    # Step 2: Orient the graph away from root vertex
    try:
        oriented_dm = orient_away_from_vertex(source_graph, root_vertex)
    except PhyloZooError as e:
        raise PhyloZooValueError(
            f"Failed to orient network away from root location {root_location}: {e}"
        )

    # Step 3: Filter attributes on the oriented graph. It is a fresh object that only
    # this function holds, and orient_away_from_vertex re-packs attributes into new
    # dicts, so editing them here cannot affect the source network.
    oriented_dm._graph.graph.clear()
    for _node, node_attrs in oriented_dm._graph.nodes(data=True):
        for attr_key in [key for key in node_attrs if key != "label"]:
            del node_attrs[attr_key]
    allowed_edge_attrs = {"gamma", "branch_length"}
    for _u, _v, _key, edge_attrs in oriented_dm._graph.edges(keys=True, data=True):
        for attr_key in [key for key in edge_attrs if key not in allowed_edge_attrs]:
            del edge_attrs[attr_key]
    oriented_dm._combined_cache = None

    # Step 4: Convert to DirectedPhyNetwork
    try:
        if validate:
            return dnetwork_from_graph(oriented_dm, copy=False)
        with no_validation():
            return dnetwork_from_graph(oriented_dm, copy=False)
    except PhyloZooError as e:
        raise PhyloZooValueError(f"Failed to convert oriented network to DirectedPhyNetwork: {e}")


def to_d_network(
    network: SemiDirectedPhyNetwork,
    root_location: RootLocation | None = None,
) -> "DirectedPhyNetwork":
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
    # A root location chosen here comes from root_locations and gives a valid network
    # by construction; only a user-supplied location needs the result validated.
    user_supplied = root_location is not None
    if root_location is None:
        node_locs: Any
        undir_edge_locs: Any
        dir_edge_locs: Any
        node_locs, undir_edge_locs, dir_edge_locs = root_locations(network)
        all_valid_locations: list[RootLocation] = (
            list(node_locs) + list(undir_edge_locs) + list(dir_edge_locs)
        )
        if not all_valid_locations:
            raise PhyloZooValueError("No valid root locations found for the network")
        root_location = all_valid_locations[0]

    return _root_sd_network_at(network, root_location, validate=user_supplied)


def root_at_outgroup(
    network: SemiDirectedPhyNetwork,
    outgroup: str,
) -> "DirectedPhyNetwork":
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
            f"No edge found incident to outgroup leaf node {outgroup_node} " f"(taxon '{outgroup}')"
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
