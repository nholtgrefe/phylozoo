"""
Split system algorithms module.

This module provides algorithms for working with split systems, including
conversion to phylogenetic networks and computation of distance matrices.
"""

import itertools
from collections import defaultdict
from typing import TYPE_CHECKING, Any

import numpy as np

from ...utils.exceptions import PhyloZooValueError
from ..distance import DistanceMatrix
from ..network.sdnetwork import SemiDirectedPhyNetwork
from ..network.sdnetwork.conversions import sdnetwork_from_graph
from ..primitives.m_multigraph import MixedMultiGraph
from ..quartet import Quartet, QuartetProfile, QuartetProfileSet
from .base import Split
from .classifications import is_tree_compatible
from .splitsystem import SplitSystem
from .weighted_splitsystem import WeightedSplitSystem

if TYPE_CHECKING:
    pass


def induced_quartetsplits(split: Split, include_trivial: bool = False) -> set[Split]:
    """
    Return a set of all subsplits of size 4 of the split.

    Generates all quartet splits (2|2 splits) that can be induced from this
    split by selecting 2 elements from each side.

    Parameters
    ----------
    split : Split[T]
        The split to generate quartet splits from.
    include_trivial : bool, optional
        If True, also include trivial quartet splits (1|3 splits).
        By default False.

    Returns
    -------
    set[Split[T]]
        A set of quartet splits induced from this split.

    Examples
    --------
    >>> split = Split({1, 2, 3}, {4, 5, 6})
    >>> quartets = induced_quartetsplits(split)
    >>> len(quartets) > 0
    True

    """
    res: list[Split] = []

    # Generate 2|2 splits
    for s1 in itertools.combinations(split.set1, 2):
        for s2 in itertools.combinations(split.set2, 2):
            quartet_split = Split(set(s1), set(s2))
            res.append(quartet_split)

    # Optionally include trivial splits (1|3)
    if include_trivial:
        for s1 in itertools.combinations(split.set1, 1):
            for s2 in itertools.combinations(split.set2, 3):
                quartet_split = Split(set(s1), set(s2))
                res.append(quartet_split)

        for s1 in itertools.combinations(split.set1, 3):
            for s2 in itertools.combinations(split.set2, 1):
                quartet_split = Split(set(s1), set(s2))
                res.append(quartet_split)

    return set(res)


def tree_from_splitsystem(
    system: SplitSystem,
    check_compatibility: bool = True,
) -> SemiDirectedPhyNetwork:
    """
    Convert a split system to a tree (SemiDirectedPhyNetwork).

    Builds the tree by clusters: with a reference taxon fixed, every non-trivial
    split is read as the set of taxa on the side away from it. Compatible clusters
    are nested or disjoint, so, inserted from smallest to largest, each one's taxa
    already hang under a single node as complete subtrees, which are gathered
    under a new internal node. Starting from a star tree this produces the tree
    inducing all splits in O(k n) for k splits on n taxa.

    Parameters
    ----------
    system : SplitSystem
        The split system to convert to a tree.
    check_compatibility : bool, optional
        Whether to check if the system is compatible with a tree before
        building. If False, assumes compatibility (e.g., if known by construction).
        By default True.

    Returns
    -------
    SemiDirectedPhyNetwork
        A tree network displaying all splits in the system.

    Raises
    ------
    PhyloZooValueError
        If check_compatibility is True and the system is not tree-compatible.
        If a split cannot be created (indicating incompatibility).

    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> from phylozoo.core.network.sdnetwork.classifications import is_tree
    >>> split1 = Split({1, 2}, {3, 4})
    >>> split2 = Split({1}, {2, 3, 4})
    >>> split3 = Split({2}, {1, 3, 4})
    >>> split4 = Split({3}, {1, 2, 4})
    >>> split5 = Split({4}, {1, 2, 3})
    >>> system = SplitSystem([split1, split2, split3, split4, split5])
    >>> tree = tree_from_splitsystem(system)
    >>> is_tree(tree)
    True
    """
    if check_compatibility:
        if not is_tree_compatible(system):
            raise PhyloZooValueError("Split system is not compatible with a tree")

    # Handle empty system
    if len(system.elements) == 0:
        return SemiDirectedPhyNetwork()

    # Handle single element
    if len(system.elements) == 1:
        element = next(iter(system.elements))
        return SemiDirectedPhyNetwork(nodes=[(element, {"label": str(element)})])

    taxa = frozenset(system.elements)

    # Root the construction at a reference taxon and read every non-trivial split as
    # the cluster on the side away from it. Compatible clusters are nested or
    # disjoint, so processed smallest first, the taxa of a cluster always hang under
    # one node as a few complete subtrees: gather those subtrees under a new node.
    # Each split costs O(|cluster|), with no cut-vertex search or graph copying.
    reference = min(taxa, key=str)
    clusters: list[tuple[Split, frozenset[Any]]] = []
    for split in system.splits:
        if split.is_trivial:
            continue
        side = split.set2 if reference in split.set1 else split.set1
        clusters.append((split, frozenset(side)))
    clusters.sort(key=lambda item: (len(item[1]), sorted(map(str, item[1]))))

    center_node = "_center"
    parent: dict[Any, Any] = {taxon: center_node for taxon in taxa}
    # top[t]: the highest already-created node whose cluster contains taxon t.
    top: dict[Any, Any] = {taxon: taxon for taxon in taxa}
    cluster_size: dict[Any, int] = {taxon: 1 for taxon in taxa}
    internal_nodes: list[Any] = []

    for split, cluster in clusters:
        tops = {top[taxon] for taxon in cluster}
        parents = {parent[node] for node in tops}
        # Every gathered subtree must lie wholly inside the cluster and all must
        # share a parent; otherwise this split conflicts with one already placed.
        if len(parents) != 1 or sum(cluster_size[node] for node in tops) != len(cluster):
            raise PhyloZooValueError(
                f"Could not find a cut-vertex to create split {split}. "
                "This indicates the split system may not be compatible."
            )
        new_node = f"_i{len(internal_nodes)}"
        internal_nodes.append(new_node)
        parent[new_node] = next(iter(parents))
        cluster_size[new_node] = len(cluster)
        for node in tops:
            parent[node] = new_node
        for taxon in cluster:
            top[taxon] = new_node

    # Materialise the tree.
    T: Any = MixedMultiGraph()
    T.add_node(center_node)
    for node in internal_nodes:
        T.add_node(node)
    for taxon in taxa:
        T.add_node(taxon)
        if taxon not in T._undirected:
            T._undirected.add_node(taxon)
        T._undirected.nodes[taxon]["label"] = str(taxon)
    for child, node in parent.items():
        T.add_undirected_edge(node, child)

    # Convert MixedMultiGraph to SemiDirectedPhyNetwork
    return sdnetwork_from_graph(T, network_type="semi-directed", copy=False)


def distances_from_splitsystem(system: SplitSystem | WeightedSplitSystem) -> DistanceMatrix:
    """
    Compute distance matrix from a split system.

    The distance between two elements x and y is the sum of weights of all splits
    that separate x and y. A split separates x and y if one element is in set1
    and the other is in set2.

    All splits are folded into one weighted indicator matrix and the distances come
    out of a single matrix product, so the work runs in BLAS rather than as one
    n x n update per split.

    Parameters
    ----------
    system : SplitSystem | WeightedSplitSystem
        The split system. If WeightedSplitSystem, split weights are used.
        If SplitSystem, each split has implicit weight 1.0.

    Returns
    -------
    DistanceMatrix
        A distance matrix on the elements of the split system, where the distance
        between x and y is the sum of weights of splits that separate them.

    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> split1 = Split({1, 2}, {3, 4})
    >>> split2 = Split({1, 3}, {2, 4})
    >>> weights = {split1: 2.0, split2: 1.5}
    >>> system = WeightedSplitSystem(weights)
    >>> dm = distances_from_splitsystem(system)
    >>> dm.get_distance(1, 2)  # Separated by split2 only
    1.5
    >>> dm.get_distance(1, 3)  # Separated by split1 only
    2.0
    >>> dm.get_distance(1, 4)  # Separated by both splits
    3.5
    >>> dm.get_distance(2, 3)  # Separated by both splits
    3.5
    >>> # Unweighted split system (each split has weight 1.0)
    >>> system2 = SplitSystem([split1, split2])
    >>> dm2 = distances_from_splitsystem(system2)
    >>> dm2.get_distance(1, 4)  # Separated by both splits, each with weight 1.0
    2.0
    """
    # Handle empty system
    if len(system.elements) == 0:
        return DistanceMatrix(np.zeros((0, 0), dtype=np.float64), labels=[])

    # Check if system is weighted
    is_weighted = isinstance(system, WeightedSplitSystem)

    # Get all elements as a sorted list for consistent ordering
    elements_list = sorted(system.elements)
    n = len(elements_list)
    index = {element: position for position, element in enumerate(elements_list)}

    # With b_s(x) = 1 if x is on set1 of split s, a pair is separated by s exactly
    # when b_s(x) != b_s(y), i.e. b_s(x) + b_s(y) - 2 b_s(x) b_s(y) for 0/1 values. So
    #     d(x, y) = r(x) + r(y) - 2 * sum_s w_s b_s(x) b_s(y),   r(x) = sum_s w_s b_s(x)
    # which is one matrix product over an indicator matrix instead of an n x n update
    # per split. Splits are processed in blocks to bound the indicator's size.
    gram = np.zeros((n, n), dtype=np.float64)
    row_weight = np.zeros(n, dtype=np.float64)
    splits = list(system.splits)
    block_rows = max(1, min(len(splits), 4_000_000 // max(n, 1)))

    for start in range(0, len(splits), block_rows):
        block = splits[start : start + block_rows]
        indicator = np.zeros((len(block), n), dtype=np.float64)
        weights = np.empty(len(block), dtype=np.float64)
        for row, split in enumerate(block):
            indicator[row, [index[element] for element in split.set1]] = 1.0
            weights[row] = system.get_weight(split) if is_weighted else 1.0
        weighted = indicator * weights[:, None]
        row_weight += weighted.sum(axis=0)
        gram += indicator.T @ weighted

    distance_matrix = row_weight[:, None] + row_weight[None, :] - 2.0 * gram
    # Exact symmetry and zero diagonal (the product is symmetric only up to rounding).
    distance_matrix = np.triu(distance_matrix, k=1)
    distance_matrix = distance_matrix + distance_matrix.T

    # Create DistanceMatrix
    return DistanceMatrix(distance_matrix, labels=elements_list)


def quartets_from_splitsystem(system: SplitSystem | WeightedSplitSystem) -> QuartetProfileSet:
    """
    Compute quartet profile set from a split system.

    For each split in the system, this function extracts all quartets induced by it
    (all 2|2 splits: 2 elements from one side, 2 from the other). Quartets are then
    grouped by their 4-taxon set into profiles, with weights equal to how often each
    quartet appeared (summing weights if the system is weighted).

    Parameters
    ----------
    system : SplitSystem | WeightedSplitSystem
        The split system. If WeightedSplitSystem, split weights are used.
        If SplitSystem, each split has implicit weight 1.0.

    Returns
    -------
    QuartetProfileSet
        A quartet profile set where each profile corresponds to a 4-taxon set,
        and contains quartets weighted by how often they appeared in the splits.

    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> split1 = Split({1, 2, 3}, {4, 5, 6})
    >>> split2 = Split({1, 2}, {3, 4, 5, 6})
    >>> system = SplitSystem([split1, split2])
    >>> profileset = quartets_from_splitsystem(system)
    >>> len(profileset) > 0
    True
    """
    # Handle empty system or system with fewer than 4 elements
    if len(system.elements) < 4:
        return QuartetProfileSet()

    # Check if system is weighted
    is_weighted = isinstance(system, WeightedSplitSystem)

    # Collect quartets with their weights, grouped by 4-taxon set
    # Use defaultdict to avoid nested if checks for efficiency
    # Structure: {frozenset(taxa): {Quartet: total_weight}}
    profile_data: dict[frozenset[Any], dict[Quartet, float]] = defaultdict(dict)

    # Process each split
    for split in system.splits:
        # Skip splits that can't produce quartets (need at least 2 elements on each side)
        if split.is_trivial:
            continue

        # Get weight for this split (1.0 if not weighted, or actual weight if weighted)
        split_weight = system.get_weight(split) if is_weighted else 1.0

        # Get all quartet splits induced by this split
        quartet_splits = induced_quartetsplits(split, include_trivial=False)

        # Convert each quartet split to a Quartet and add to profile
        for quartet_split in quartet_splits:
            # Create Quartet from the split
            quartet = Quartet(quartet_split)
            quartet_taxa = quartet.taxa

            # Add weight (sum if quartet appears multiple times)
            # Using defaultdict, we don't need to check if quartet_taxa exists
            profile_data[quartet_taxa][quartet] = (
                profile_data[quartet_taxa].get(quartet, 0.0) + split_weight
            )

    # Create QuartetProfile objects from the grouped data.
    # QuartetProfile requires weights to sum to 1.0; normalize accumulated weights per profile.
    profiles: list[QuartetProfile] = []
    for quartet_taxa, quartets_dict in profile_data.items():
        if quartets_dict:  # Only create profile if it has quartets
            total = sum(quartets_dict.values())
            normalized = {q: w / total for q, w in quartets_dict.items()}
            profile = QuartetProfile(normalized)
            profiles.append(profile)

    # Create QuartetProfileSet (each profile gets default weight 1.0)
    return QuartetProfileSet(profiles=profiles)
