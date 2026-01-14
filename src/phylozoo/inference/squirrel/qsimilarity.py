"""
Quartet similarity module.

This module provides functions for computing similarity between squirrel quartet
profile sets and for creating profile sets from networks.
"""

from typing import TYPE_CHECKING, Any, Iterator

from ...utils.exceptions import PhyloZooValueError, PhyloZooAlgorithmError

from ...core.quartet.base import Quartet
from ...core.network.sdnetwork.derivations import induced_splits, partition_from_blob
from ...core.network.sdnetwork.classifications import level, is_binary, has_parallel_edges
from ...core.network.sdnetwork.features import blobs
from ...core.split.base import Split, induced_quartetsplits
from ...core.primitives.circular_ordering import CircularSetOrdering

from .sqprofile import SqQuartetProfile
from .sqprofileset import SqQuartetProfileSet



if TYPE_CHECKING:
    from ...core.network.sdnetwork import MixedPhyNetwork, SemiDirectedPhyNetwork


def _circular_orders_from_cycles(
    network: 'MixedPhyNetwork',
) -> Iterator[tuple[CircularSetOrdering, frozenset[str] | None]]:
    """
    Yield circular set ordering and reticulation set for each cycle in the network.
    
    This function finds all cycles in the network using NetworkX, then for each cycle,
    computes the circular set ordering of leaves induced by the cycle and the
    reticulation set (if the cycle contains a hybrid node).
    
    Parameters
    ----------
    network : MixedPhyNetwork
        The mixed or semi-directed network. Must be binary and level-1 or lower.
    
    Yields
    ------
    tuple[CircularSetOrdering, frozenset[str] | None]
        For each cycle, yields a tuple containing:
        - The circular set ordering of leaves induced by the cycle
        - The set of leaves below the hybrid (reticulation set), or None if no hybrid
    
    Raises
    ------
    PhyloZooValueError
        If the network is not binary.
        If the network has parallel edges.
        If the network is not level-1 or lower (i.e., level > 1).
        If a cycle has fewer than 3 nodes.
        If a partition cannot be created for a cycle.
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> 
    >>> # Create a simple network with a cycle
    >>> net = SemiDirectedPhyNetwork(
    ...     directed_edges=[(5, 4), (6, 4)],
    ...     undirected_edges=[(5, 3), (5, 6), (6, 7), (4, 8), (8, 1), (8, 2)],
    ...     nodes=[(3, {'label': 'C'}), (7, {'label': 'D'}), (1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> 
    >>> # Get circular set orderings for all cycles
    >>> results = list(_circular_orders_from_cycles(net))
    >>> len(results)
    1
    >>> set_ordering, ret_set = results[0]
    >>> len(set_ordering)
    3
    >>> ret_set is not None
    True
    """
    # Check if network is binary
    if not is_binary(network):
        raise PhyloZooValueError(
            "Network must be binary for _circular_orders_from_cycles"
        )
    
    # Check if network has parallel edges
    if has_parallel_edges(network):
        raise PhyloZooValueError(
            "Network must not have parallel edges for _circular_orders_from_cycles"
        )
    
    # Check if network is level-1 or lower
    net_level = level(network)
    if net_level > 1:
        raise PhyloZooValueError(
            f"Network must be level-1 or lower for _circular_orders_from_cycles, "
            f"but got level {net_level}"
        )
    
    def _find_cycle_in_blob(
        graph: Any,
        blob: set[Any],
    ) -> list[Any]:
        """
        Find a cycle in a blob by traversing neighbors.
        
        Parameters
        ----------
        graph : Any
            The combined graph (undirected view).
        blob : set[Any]
            Set of nodes forming the blob.
        
        Returns
        -------
        list[Any]
            List of nodes in cycle order.
        """
        # Start at one node and traverse neighbors to find the cycle
        start_node = next(iter(blob))
        cycle_nodes: list[Any] = [start_node]
        visited: set[Any] = {start_node}
        
        # Traverse the blob to find the cycle
        current = start_node
        previous = None
        
        while len(cycle_nodes) < len(blob):
            # Find a neighbor within the blob that we haven't visited yet
            # and is not the previous node
            next_node: Any | None = None
            
            for neighbor in graph.neighbors(current):
                if neighbor in blob and neighbor != previous:
                    next_node = neighbor
                    break
            
            if next_node is None:
                # This should not happen in a valid cycle blob
                raise PhyloZooAlgorithmError(
                    f"Could not find a valid next node in blob {blob} from {current}. "
                    "Blob may not form a valid cycle."
                )
            
            cycle_nodes.append(next_node)
            visited.add(next_node)
            previous = current
            current = next_node
        
        # After traversing all nodes, the last node must connect back to the start_node
        if start_node not in graph.neighbors(current):
            raise PhyloZooAlgorithmError(
                f"Cycle in blob {blob} is not closed. Last node {current} does not connect to start node {start_node}."
            )
        
        return cycle_nodes
    
    # Get all non-trivial, non-leaf blobs
    # In level-1 networks, each blob forms a cycle
    network_blobs = blobs(network, trivial=False, leaves=False)
    
    # Process each blob to extract its cycle
    combined_graph = network._graph._combined
    
    for blob in network_blobs:
        if len(blob) < 3:
            raise PhyloZooAlgorithmError(
                f"Blob {blob} has fewer than 3 nodes. "
                "Level-1 networks should not have such blobs."
            )

        # Process blobs with 3 or more nodes
        # Note: 3-node cycles are returned but should be filtered when used for cycle quartets

        # Find cycle in the blob
        cycle_nodes = _find_cycle_in_blob(combined_graph, blob)
        
        # Find the hybrid node in the cycle (if any)
        network_hybrid_nodes = network.hybrid_nodes
        cycle_hybrids = [node for node in cycle_nodes if node in network_hybrid_nodes]
        
        hybrid: Any | None = cycle_hybrids[0] if cycle_hybrids else None
        
        # Get partition and edge_taxa_list from the cycle blob
        _, edge_taxa_list = partition_from_blob(
            network, blob, return_edge_taxa=True
        )
        
        # Create mapping from blob node (v) to position in cycle
        cycle_order_dict = {node: i for i, node in enumerate(cycle_nodes)}
        
        # Sort edge_taxa_list by the blob node's position in the cycle
        # edge_taxa_list contains (u, v, taxa_set) where v is in the blob (cycle)
        sorted_edge_taxa = sorted(
            edge_taxa_list,
            key=lambda x: cycle_order_dict.get(x[1], -1)
        )
        
        # Extract taxa sets in cycle order
        taxa_sets = [taxa_set for _, _, taxa_set in sorted_edge_taxa]
        
        # Find the reticulation set (set below the hybrid) if hybrid exists
        ret_set: frozenset[str] | None = None
        
        if hybrid is not None:
            for (_, v, taxa_set) in sorted_edge_taxa:
                if v == hybrid:
                    ret_set = taxa_set
                    break
            
            if ret_set is None:
                raise PhyloZooAlgorithmError(
                    f"Could not find reticulation set for hybrid {hybrid} in cycle {blob}"
                )
        
        # Create circular set ordering
        set_ordering = CircularSetOrdering(taxa_sets)
        
        yield set_ordering, ret_set


def sqprofileset_similarity(
    profileset1: SqQuartetProfileSet,
    profileset2: SqQuartetProfileSet,
    weighted: bool = True,
) -> float:
    """
    Compute similarity between two squirrel quartet profile sets.
    
    Returns the consistency measure (C-measure), which is the ratio of consistent
    profiles. Two profiles are consistent if they are equal (same quartets with same
    weights). The measure can be weighted (using profile weights) or unweighted.
    
    Parameters
    ----------
    profileset1 : SqQuartetProfileSet
        First squirrel quartet profile set.
    profileset2 : SqQuartetProfileSet
        Second squirrel quartet profile set.
    weighted : bool, optional
        Whether to use weighted consistency (using profile weights).
        By default True.
    
    Returns
    -------
    float
        Similarity measure between 0.0 and 1.0. Returns 1.0 if profileset1 is empty
        (to avoid division by zero).
    
    Raises
    ------
    PhyloZooValueError
        If the two profile sets have different taxa sets.
    
    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> from phylozoo.core.quartet.base import Quartet
    >>> from phylozoo.inference.squirrel.sqprofile import SqQuartetProfile
    >>> 
    >>> # Create two profile sets
    >>> q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
    >>> q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
    >>> 
    >>> profile1 = SqQuartetProfile([q1])
    >>> profile2 = SqQuartetProfile([q1])  # Same as profile1
    >>> profile3 = SqQuartetProfile([q2])  # Different from profile1
    >>> 
    >>> set1 = SqQuartetProfileSet(profiles=[profile1, profile3])
    >>> set2 = SqQuartetProfileSet(profiles=[profile2, profile3])
    >>> 
    >>> # Compute similarity
    >>> similarity = sqprofileset_similarity(set1, set2)
    >>> 0.0 <= similarity <= 1.0
    True
    """
    # Validate that both sets have the same taxa
    if profileset1.taxa != profileset2.taxa:
        raise PhyloZooValueError(
            f"Profile sets must have the same taxa. "
            f"Set 1 taxa: {profileset1.taxa}, Set 2 taxa: {profileset2.taxa}"
        )
    
    # If profileset1 is empty, return 1.0 (to avoid division by zero)
    if len(profileset1) == 0:
        return 1.0
    
    # Collect consistent profiles
    if weighted:
        # Weighted consistency: sum weights of consistent profiles
        consistent_weight = 0.0
        total_weight = 0.0
        
        for taxa_set, (profile1, weight1) in profileset1.profiles.items():
            total_weight += weight1
            
            # Check if profileset2 has a profile for the same taxa
            profile2_result = profileset2.get_profile(taxa_set)
            if profile2_result is not None:
                profile2, weight2 = profileset2.profiles[taxa_set]
                
                # Check if profiles are equal (compare quartets and reticulation_leaf)
                if (
                    profile1._quartets == profile2._quartets
                    and profile1.reticulation_leaf == profile2.reticulation_leaf
                ):
                    # Use weight from profileset1 (matching old implementation)
                    consistent_weight += weight1
        
        if total_weight == 0:
            return 1.0
        
        return consistent_weight / total_weight
    else:
        # Unweighted consistency: count consistent profiles
        consistent_count = 0
        total_count = len(profileset1)
        
        for taxa_set, (profile1, _) in profileset1.profiles.items():
            # Check if profileset2 has a profile for the same taxa
            profile2 = profileset2.get_profile(taxa_set)
            if profile2 is not None and (
                profile1._quartets == profile2._quartets
                and profile1.reticulation_leaf == profile2.reticulation_leaf
            ):
                consistent_count += 1
        
        return consistent_count / total_count if total_count > 0 else 1.0


def sqprofileset_from_network(
    network: 'SemiDirectedPhyNetwork',
) -> SqQuartetProfileSet:
    """
    Compute squirrel quartet profile set from a semi-directed network.
    
    This function requires the network to be level-1 and uses an optimized approach:
    it extracts displayed splits and converts them to quartets, which is faster
    than computing all displayed trees for each 4-taxon set.
    
    Each profile contains either 1 or 2 resolved quartets (no unresolved quartets).
    For profiles with 2 quartets, the reticulation_leaf is determined from cycles
    in the network when possible.
    
    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network. Must be level-1 or lower.
    
    Returns
    -------
    SqQuartetProfileSet
        A squirrel quartet profile set containing profiles for each 4-taxon set
        that has at least one resolved quartet.
    
    Raises
    ------
    PhyloZooValueError
        If the network is not binary.
        If the network has parallel edges.
        If the network is not level-1 or lower (i.e., level > 1).
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> 
    >>> # Create a simple network
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 4)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
    ... )
    >>> 
    >>> # Network has fewer than 4 taxa, so result is empty
    >>> profileset = sqprofileset_from_network(net)
    >>> len(profileset)
    0
    """
    # Handle networks with fewer than 4 taxa
    taxa_list = sorted(network.taxa)
    if len(taxa_list) < 4:
        return SqQuartetProfileSet()
    
    # Check if network is binary
    if not is_binary(network):
        raise PhyloZooValueError(
            "Network must be binary for sqprofileset_from_network"
        )
    
    # Check if network has parallel edges
    if has_parallel_edges(network):
        raise PhyloZooValueError(
            "Network must not have parallel edges for sqprofileset_from_network"
        )
    
    # Check if network is level-1 (required for this optimization)
    net_level = level(network)
    if net_level > 1:
        raise PhyloZooValueError(
            f"Network must be level-1 or lower for sqprofileset_from_network, "
            f"but got level {net_level}"
        )
    
    # Collect quartets for each 4-taxon set
    # Build a list of SqQuartetProfile objects
    sq_quartet_profiles: list[SqQuartetProfile] = []
    
    # Step 1: Process cycles to get quartets from cycles
    # Get circular set orderings and reticulation sets for all cycles
    # Only consider cycles with 4+ sets (needed for cycle quartets)
    for set_ordering, ret_set in _circular_orders_from_cycles(network):
        # Skip cycles with fewer than 4 sets (can't form cycle quartets)
        if len(set_ordering) < 4:
            continue
        # Get 4-suborderings and create quartets
        for sub_order in set_ordering.suborderings(4):
            # If no reticulation set, create a four-cycle (i.e. 2 quartets)
            # without a reticulation leaf for each representative ordering
            if ret_set is None:
                for repr_order in sub_order.representative_orderings():
                    # Create both quartets from the circular ordering
                    order_list = list(repr_order.order)

                    # For a circular ordering [a, b, c, d], the quartets are:
                    # (a b | c d) and (a d | b c)
                    a, b, c, d = order_list
                    q1 = Quartet(Split({a, b}, {c, d}))
                    q2 = Quartet(Split({a, d}, {b, c}))
                    
                    sq_profile = SqQuartetProfile([q1, q2])
                    sq_quartet_profiles.append(sq_profile)
            
            # If there is a reticulation set and it is in the sub-ordering, create a 
            # four-cycle with a reticulation leaf for each representative ordering
            # Check if ret_set is one of the sets in sub_order
            elif ret_set is not None and any(ret_set == part for part in sub_order.parts):
                for repr_order in sub_order.representative_orderings():
                    # Create both quartets from the circular ordering
                    order_list = list(repr_order.order)

                    # Find the reticulation leaf in the sub-ordering
                    # repr_order is a CircularOrdering, which has an elements property from Partition
                    ret_leaf = next(iter(ret_set & repr_order.elements))

                    # Create both quartets from the circular ordering
                    a, b, c, d = order_list
                    q1 = Quartet(Split({a, b}, {c, d}))
                    q2 = Quartet(Split({a, d}, {b, c}))
                    sq_profile = SqQuartetProfile([q1, q2], reticulation_leaf=ret_leaf)
                    sq_quartet_profiles.append(sq_profile)
            
            # If there is a reticulation set and it is not in the sub-ordering, create a 
            # quartet tree for each representative ordering
            elif ret_set is not None and not any(ret_set == part for part in sub_order.parts):
                for repr_order in sub_order.representative_orderings():
                    # Get the split from the sub-ordering (i.e. the 2|2 split)
                    # by taking the first two and last two elements of the ordering
                    a, b, c, d = list(repr_order.order)
                    split = Split({a, b}, {c, d})
                    quartet = Quartet(split)
                    sq_profile = SqQuartetProfile([quartet])
                    sq_quartet_profiles.append(sq_profile)
    
    # Step 2: Process induced splits to get quartets from splits
    # Get all induced splits with weights
    split_system = induced_splits(network)
    seen_quartet_splits: set[Split] = set()
    # Extract all 2|2 splits on 4 taxa (induced quartetsplits)    
    for split in split_system.splits:
        # Skip trivial splits and splits we have already seen
        if split.is_trivial():
            continue

        # Get all quartet splits induced by this split
        for quartet_split in induced_quartetsplits(split):
            if quartet_split in seen_quartet_splits:
                continue
            seen_quartet_splits.add(quartet_split)
            quartet = Quartet(quartet_split)
            sq_profile = SqQuartetProfile([quartet])
            sq_quartet_profiles.append(sq_profile)
    
    # Create SqQuartetProfileSet
    return SqQuartetProfileSet(profiles=sq_quartet_profiles)

