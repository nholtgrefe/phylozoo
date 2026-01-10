"""
Quartet similarity module.

This module provides functions for computing similarity between squirrel quartet
profile sets and for creating profile sets from networks.
"""

from typing import TYPE_CHECKING, Any, Iterator

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
    ValueError
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
        raise ValueError(
            "Network must be binary for _circular_orders_from_cycles"
        )
    
    # Check if network has parallel edges
    if has_parallel_edges(network):
        raise ValueError(
            "Network must not have parallel edges for _circular_orders_from_cycles"
        )
    
    # Check if network is level-1 or lower
    net_level = level(network)
    if net_level > 1:
        raise ValueError(
            f"Network must be level-1 or lower for _circular_orders_from_cycles, "
            f"but got level {net_level}"
        )
    
    # Get all non-trivial, non-leaf blobs
    # In level-1 networks, each blob forms a cycle
    network_blobs = blobs(network, trivial=False, leaves=False)
    
    # Process each blob to extract its cycle
    cycles: list[list[Any]] = []
    combined_graph = network._graph._combined
    
    for blob in network_blobs:
        if len(blob) < 3:
            raise ValueError(
                f"Blob {blob} has fewer than 3 nodes. "
                "Level-1 networks should not have such blobs."
            )
        
        # Find a cycle in the blob by starting at one node and traversing neighbors
        # In a cycle, each node has exactly 2 neighbors within the cycle
        # Start at an arbitrary node in the blob
        start_node = next(iter(blob))
        cycle_nodes: list[Any] = [start_node]
        prev_node: Any | None = None
        current = start_node
        
        # Traverse the blob to find the cycle
        # At each step, pick the neighbor that's not the previous node
        while len(cycle_nodes) < len(blob):
            # Get neighbors of current node within the blob
            neighbors = [
                n for n in combined_graph.neighbors(current)
                if n in blob
            ]
            
            # In a cycle, each node should have exactly 2 neighbors
            if len(neighbors) != 2:
                raise ValueError(
                    f"Node {current} in blob {blob} has {len(neighbors)} neighbors, "
                    "expected 2 for a cycle in a level-1 network."
                )
            
            # Pick the neighbor that's not the previous node
            if prev_node is None:
                # First step: pick either neighbor
                next_node = neighbors[0]
            else:
                # Subsequent steps: pick the neighbor that's not prev_node
                next_node = neighbors[1] if neighbors[0] == prev_node else neighbors[0]
            
            cycle_nodes.append(next_node)
            prev_node = current
            current = next_node
        
        # Verify the cycle is closed (last node connects to first)
        if cycle_nodes[-1] not in combined_graph.neighbors(cycle_nodes[0]):
            raise ValueError(
                f"Blob {blob} does not form a closed cycle. "
                f"Last node {cycle_nodes[-1]} does not connect to first node {cycle_nodes[0]}."
            )
        
        cycles.append(cycle_nodes)
    
    # Process each cycle
    for cycle_nodes in cycles:
        
        cycle_vertices = set(cycle_nodes)
        
        # Find the hybrid node in the cycle (if any)
        network_hybrid_nodes = network.hybrid_nodes
        cycle_hybrids = [node for node in cycle_nodes if node in network_hybrid_nodes]
        
        hybrid: Any | None = cycle_hybrids[0] if cycle_hybrids else None
        
        # Get partition and edge_taxa_list from the cycle blob
        _, edge_taxa_list = partition_from_blob(
            network, cycle_vertices, return_edge_taxa=True
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
                raise ValueError(
                    f"Could not find reticulation set for hybrid {hybrid} in cycle {cycle_vertices}"
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
    ValueError
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
        raise ValueError(
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
    ValueError
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
        raise ValueError(
            "Network must be binary for sqprofileset_from_network"
        )
    
    # Check if network has parallel edges
    if has_parallel_edges(network):
        raise ValueError(
            "Network must not have parallel edges for sqprofileset_from_network"
        )
    
    # Check if network is level-1 (required for this optimization)
    net_level = level(network)
    if net_level > 1:
        raise ValueError(
            f"Network must be level-1 or lower for sqprofileset_from_network, "
            f"but got level {net_level}"
        )
    
    # Collect quartets for each 4-taxon set
    # Build a list of SqQuartetProfile objects
    sq_quartet_profiles: list[SqQuartetProfile] = []
    
    # Step 1: Process cycles to get quartets from cycles
    # Get circular set orderings and reticulation sets for all cycles
    for set_ordering, ret_set in _circular_orders_from_cycles(network):
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
            elif ret_set is not None and ret_set in sub_order:
                for repr_order in sub_order.representative_orderings():
                    # Create both quartets from the circular ordering
                    order_list = list(repr_order.order)

                    # Find the reticulation leaf in the sub-ordering
                    ret_leaf = next(iter(ret_set & repr_order.elements))

                    # Create both quartets from the circular ordering
                    a, b, c, d = order_list
                    q1 = Quartet(Split({a, b}, {c, d}))
                    q2 = Quartet(Split({a, d}, {b, c}))
                    sq_profile = SqQuartetProfile([q1, q2], reticulation_leaf=ret_leaf)
                    sq_quartet_profiles.append(sq_profile)
            
            # If there is a reticulation set and it is not in the sub-ordering, create a 
            # quartet tree for each representative ordering
            elif ret_set is not None and ret_set not in sub_order:
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
    
    # Extract all 2|2 splits on 4 taxa (induced quartetsplits)    
    for split in split_system.splits:
        # Skip trivial splits
        if split.is_trivial():
            continue

        # Get all quartet splits induced by this split
        for quartet_split in induced_quartetsplits(split):
            quartet = Quartet(quartet_split)
            sq_profile = SqQuartetProfile([quartet])
            sq_quartet_profiles.append(sq_profile)
    
    # Create SqQuartetProfileSet
    return SqQuartetProfileSet(profiles=sq_quartet_profiles)

