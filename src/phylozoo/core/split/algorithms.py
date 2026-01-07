"""
Split system algorithms module.

This module provides algorithms for working with split systems, such as
converting split systems to trees.

TODO: 
- Use refinement approach clusters for splitsystem to tree function.
- circular splitsystem classificaiton
"""

import networkx as nx
from collections import defaultdict
from typing import TYPE_CHECKING, Any

import numpy as np

from ..distance import DistanceMatrix
from ..network.sdnetwork import SemiDirectedPhyNetwork
from ..network.sdnetwork.conversions import sdnetwork_from_graph
from ..primitives.m_multigraph import MixedMultiGraph
from ..primitives.m_multigraph.features import cut_vertices
from ..primitives.partition import Partition
from ..quartet import Quartet, QuartetProfile, QuartetProfileSet
from .classifications import is_tree_compatible
from .splitsystem import SplitSystem
from .weighted_splitsystem import WeightedSplitSystem

if TYPE_CHECKING:
    pass


def splitsystem_to_tree(
    system: SplitSystem,
    check_compatibility: bool = True,
) -> SemiDirectedPhyNetwork:
    """
    Convert a split system to a tree (SemiDirectedPhyNetwork).
    
    Builds a tree that induces all splits in the system using a star tree approach:
    1. Start with a star tree (center node connected to all leaves)
    2. For each non-trivial split S = A|B:
       - Find a cut-vertex v whose partition is a refinement of S
       - Replace v with a cut-edge (two internal nodes u and w connected by an edge)
       - Reconnect components: parts in A connect to u, parts in B connect to w
    
    This approach iteratively refines the star tree by splitting cut-vertices into
    cut-edges, creating the final tree structure that displays induces all splits.
    
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
    ValueError
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
    >>> tree = splitsystem_to_tree(system)
    >>> is_tree(tree)
    True
    """
    if check_compatibility:
        if not is_tree_compatible(system):
            raise ValueError("Split system is not compatible with a tree")
    
    # Handle empty system
    if len(system.elements) == 0:
        return SemiDirectedPhyNetwork()
    
    # Handle single element
    if len(system.elements) == 1:
        element = next(iter(system.elements))
        return SemiDirectedPhyNetwork(
            nodes=[(element, {'label': str(element)})]
        )
    
    # Get all taxa and splits
    taxa = frozenset(system.elements)
    splits_list = list(system.splits)
    
    # Initialize the graph as a star tree
    T = MixedMultiGraph()
    
    # Create center node
    center_node = "_center"
    T.add_node(center_node)
    
    # Add all taxa as nodes with labels
    for taxon in taxa:
        T.add_node(taxon)
        if taxon not in T._undirected:
            T._undirected.add_node(taxon)
        T._undirected.nodes[taxon]['label'] = str(taxon)
        # Connect each taxon to the center
        T.add_undirected_edge(center_node, taxon)
    
    # Store the next available node index for internal nodes
    next_node_index = [0]  # Use list to allow modification in nested function
    
    def _get_partition_from_cutvertex(v: Any, graph: MixedMultiGraph) -> tuple[Partition, dict[frozenset[Any], tuple[Any, Any]]]:
        """
        Get the partition induced by a cut-vertex v.
        
        Returns a tuple of:
        - The Partition object representing the partition of taxa
        - A dictionary mapping each part (frozenset of taxa) to the neighbor (u, v)
          where u is the neighbor of v in that component
        """
        partition_parts: list[set[Any]] = []
        neighbor_dict: dict[frozenset[Any], tuple[Any, Any]] = {}
        
        # Remove v temporarily to find connected components
        temp_graph = graph._combined.copy()
        temp_graph.remove_node(v)
        components = list(nx.connected_components(temp_graph))
        
        # For each component, find the taxa and the neighbor connecting it to v
        for component in components:
            # Get taxa in this component
            component_taxa = frozenset(t for t in component if t in taxa)
            if not component_taxa:
                continue
            
            # Find the neighbor of v in this component
            neighbor = None
            for u in component:
                if graph.has_edge(u, v) or graph.has_edge(v, u):
                    neighbor = u
                    break
            
            if neighbor is not None:
                partition_parts.append(set(component_taxa))
                neighbor_dict[component_taxa] = (neighbor, v)
        
        # Ensure partition covers all taxa (add missing as singletons)
        partition_taxa = set()
        for part in partition_parts:
            partition_taxa.update(part)
        missing_taxa = taxa - partition_taxa
        for taxon in missing_taxa:
            partition_parts.append({taxon})
        
        vertex_partition = Partition(partition_parts)
        return vertex_partition, neighbor_dict
    
    
    # Process each non-trivial split
    for split in splits_list:
        if split.is_trivial():
            continue
        
        # Find a cut-vertex whose partition is a refinement of this split
        found = False
        cut_vertices_set = cut_vertices(T, data=False)
        
        for v in cut_vertices_set:
            # Skip if v is a leaf (taxon)
            if v in taxa:
                continue
            
            # Get partition from cut-vertex
            vertex_partition, neighbor_dict = _get_partition_from_cutvertex(v, T)
            
            # Check if partition is a refinement of the split
            # Split is already a Partition object, so we can use is_refinement directly
            if vertex_partition.is_refinement(split):
                # Found a suitable cut-vertex! Replace it with a cut-edge
                found = True
                
                # Create two new internal nodes
                u = f"_i{next_node_index[0]}"
                next_node_index[0] += 1
                w = f"_i{next_node_index[0]}"
                next_node_index[0] += 1
                
                T.add_node(u)
                T.add_node(w)
                
                # Create the cut-edge (u, w)
                T.add_undirected_edge(u, w)
                
                # Remove the old cut-vertex v
                T.remove_node(v)
                
                # Reconnect components
                for part_frozen, (neighbor, _) in neighbor_dict.items():
                    part = set(part_frozen)
                    # Determine which side of the split this part belongs to
                    if part.issubset(split.set1):
                        # Connect to u
                        T.add_undirected_edge(u, neighbor)
                    else:  # part.issubset(split.set2)
                        # Connect to w
                        T.add_undirected_edge(w, neighbor)
                
                break
        
        if not found:
            # This should not happen if the system is compatible
            raise ValueError(
                f"Could not find a cut-vertex to create split {split}. "
                "This indicates the split system may not be compatible."
            )
    
    # Convert MixedMultiGraph to SemiDirectedPhyNetwork
    return sdnetwork_from_graph(T, network_type='semi-directed')


def distances_from_splitsystem(system: SplitSystem | WeightedSplitSystem) -> DistanceMatrix:
    """
    Compute distance matrix from a split system.
    
    The distance between two elements x and y is the sum of weights of all splits
    that separate x and y. A split separates x and y if one element is in set1
    and the other is in set2.
    
    This function uses vectorized numpy operations for efficiency, avoiding nested
    Python loops by using boolean arrays and broadcasting.
    
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
    >>> system = WeightedSplitSystem([split1, split2], weights=weights)
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
    
    # Initialize distance matrix
    distance_matrix = np.zeros((n, n), dtype=np.float64)
    
    # For each split, use vectorized operations to add weights
    for split in system.splits:
        # Get weight for this split (1.0 if not weighted, or actual weight if weighted)
        weight = system.get_weight(split) if is_weighted else 1.0
        
        # Create boolean array: True if element is in set1, False if in set2
        # This is O(n) instead of O(n^2) for checking each pair
        in_set1 = np.array([elem in split.set1 for elem in elements_list], dtype=bool)
        
        # Use broadcasting to create a boolean matrix where True means the pair is separated
        # in_set1[:, None] creates a column vector, in_set1[None, :] creates a row vector
        # The != operation broadcasts to create an n x n matrix
        separated = in_set1[:, None] != in_set1[None, :]
        
        # Add weight to all separated pairs (vectorized operation)
        distance_matrix += separated.astype(np.float64) * weight
    
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
        if split.is_trivial():
            continue
        
        # Get weight for this split (1.0 if not weighted, or actual weight if weighted)
        split_weight = system.get_weight(split) if is_weighted else 1.0
        
        # Get all quartet splits induced by this split
        quartet_splits = split.induced_quartetsplits(include_trivial=False)
        
        # Convert each quartet split to a Quartet and add to profile
        for quartet_split in quartet_splits:
            # Create Quartet from the split
            quartet = Quartet(quartet_split)
            quartet_taxa = quartet.taxa
            
            # Add weight (sum if quartet appears multiple times)
            # Using defaultdict, we don't need to check if quartet_taxa exists
            profile_data[quartet_taxa][quartet] = profile_data[quartet_taxa].get(quartet, 0.0) + split_weight
    
    # Create QuartetProfile objects from the grouped data
    profiles: list[QuartetProfile] = []
    for quartet_taxa, quartets_dict in profile_data.items():
        if quartets_dict:  # Only create profile if it has quartets
            profile = QuartetProfile(quartets_dict)
            profiles.append(profile)
    
    # Create QuartetProfileSet (each profile gets default weight 1.0)
    return QuartetProfileSet(profiles=profiles)

