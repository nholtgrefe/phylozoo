"""
Split system algorithms module.

This module provides algorithms for working with split systems, such as
converting split systems to trees.

TODO: Use refinement approach clusters for splitsystem to tree function.
"""

import networkx as nx
from typing import TYPE_CHECKING, Any

from ..network.sdnetwork import SemiDirectedPhyNetwork
from ..network.sdnetwork.conversions import sdnetwork_from_graph
from ..primitives.m_multigraph import MixedMultiGraph
from ..primitives.m_multigraph.features import cut_vertices
from ..primitives.partition import Partition
from .classifications import is_tree_compatible
from .splitsystem import SplitSystem

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

