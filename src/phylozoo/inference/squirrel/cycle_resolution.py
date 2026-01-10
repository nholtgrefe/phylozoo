"""
Cycle resolution module.

This module provides functions for resolving cycles in phylogenetic networks
by computing circular set orderings based on quartet profiles.
"""

from typing import TYPE_CHECKING, Any, Literal

from ...core.distance.operations import approximate_tsp_tour, optimal_tsp_tour
from ...core.primitives.circular_ordering import CircularSetOrdering
from ...core.primitives.partition import Partition
from ...core.quartet.qdistance import quartet_distance_with_partition

if TYPE_CHECKING:
    from ...core.quartet.qprofileset import QuartetProfileSet
    from ...core.network.sdnetwork import SemiDirectedPhyNetwork, MixedPhyNetwork
    from ...core.primitives.m_multigraph import MixedMultiGraph
    from .sqprofileset import SqQuartetProfileSet


def _qprofiles_to_circular_ordering(
    profileset: 'QuartetProfileSet',
    partition: Partition,
    rho: tuple[float, float, float, float] = (0.5, 1.0, 0.5, 1.0),
    tsp_method: Literal['optimal', 'simulated_annealing', 'greedy', 'christofides'] = 'optimal',
) -> 'CircularSetOrdering':
    """
    Compute the optimal circular set ordering from quartet profiles.
    
    This function computes a distance matrix between partition sets based on
    quartet profiles, then solves a TSP to find the optimal circular ordering
    of the sets.
    
    Parameters
    ----------
    profileset : QuartetProfileSet
        The quartet profile set to compute distances from. Must be dense.
    partition : Partition
        A partition of the taxa. The circular ordering will be computed
        for the sets (parts) of this partition.
    rho : tuple[float, float, float, float], optional
        Rho vector (rho_c, rho_s, rho_a, rho_o) specifying distance contributions.
        By default (0.5, 1.0, 0.5, 1.0) (Squirrel/MONAD).
    tsp_method : Literal['optimal', 'simulated_annealing', 'greedy', 'christofides'], optional
        TSP method to use. By default 'optimal'.
    
    Returns
    -------
    CircularSetOrdering
        A circular ordering of the partition sets.
    
    Raises
    ------
    ValueError
        If partition size is less than 3 (TSP requires at least 3 nodes).
    
    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> from phylozoo.core.quartet.base import Quartet
    >>> from phylozoo.core.quartet.qprofileset import QuartetProfileSet
    >>> from phylozoo.core.primitives.partition import Partition
    >>> 
    >>> # Create a simple profile set
    >>> q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
    >>> profileset = QuartetProfileSet(profiles=[q1])
    >>> 
    >>> # Create a partition
    >>> partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}])
    >>> 
    >>> # Compute circular ordering
    >>> ordering = _qprofiles_to_circular_ordering(profileset, partition)
    >>> len(ordering)
    4
    """
    if len(partition) < 3:
        raise ValueError(
            f"Partition must have at least 3 sets for TSP, got {len(partition)}"
        )
    
    # Compute distance matrix between partition sets
    dist_matrix = quartet_distance_with_partition(
        profileset=profileset,
        partition=partition,
        rho=rho,
    )
    
    # Solve TSP
    if tsp_method == 'optimal':
        tour = optimal_tsp_tour(dist_matrix)
    elif tsp_method == 'simulated_annealing':
        tour = approximate_tsp_tour(dist_matrix, method='simulated_annealing')
    elif tsp_method == 'greedy':
        tour = approximate_tsp_tour(dist_matrix, method='greedy')
    elif tsp_method == 'christofides':
        tour = approximate_tsp_tour(dist_matrix, method='christofides')
    else:
        raise ValueError(
            f"Invalid tsp_method: {tsp_method}. "
            "Must be one of ['optimal', 'simulated_annealing', 'greedy', 'christofides']"
        )
    
    # Convert tour (CircularOrdering of frozensets) to CircularSetOrdering
    # The tour.order contains the frozensets (which are the partition sets) in order
    # Convert frozensets to sets for CircularSetOrdering
    set_order = [set(fs) for fs in tour.order]
    return CircularSetOrdering(set_order)


def _qprofiles_to_hybrid_ranking(
    profileset: 'SqQuartetProfileSet',
    partition: Partition,
    weights: bool = True,
) -> list[frozenset[str]]:
    """
    Return an ordered list of partition sets by likelihood of being the hybrid/reticulation set.
    
    The first set in the returned list is the most likely to be the hybrid set,
    the second is the second most likely, etc. This ranking is based on voting from
    quartet profiles:
    - For 4-set partitions: uses reticulation_leaf from 4-cycle profiles
    - For larger partitions: aggregates cycle percentages from 4-subpartitions
    
    Parameters
    ----------
    profileset : SqQuartetProfileSet
        The squirrel quartet profile set to analyze.
    partition : Partition
        A partition of the taxa. The ranking will be computed for the sets
        (parts) of this partition.
    weights : bool, optional
        Whether to use quartet weights in the voting. By default True.
    
    Returns
    -------
    list[frozenset[str]]
        Ordered list of partition sets (frozensets), sorted by hybrid likelihood
        (highest first).
    
    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> from phylozoo.core.quartet.base import Quartet
    >>> from phylozoo.inference.squirrel.sqprofileset import SqQuartetProfileSet
    >>> from phylozoo.core.primitives.partition import Partition
    >>> 
    >>> # Create a simple profile set
    >>> q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
    >>> profileset = SqQuartetProfileSet(profiles=[q1])
    >>> 
    >>> # Create a partition
    >>> partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}])
    >>> 
    >>> # Get hybrid ranking
    >>> ranking = _qprofiles_to_hybrid_ranking(profileset, partition)
    >>> len(ranking)
    4
    """
    # Initialize voting dictionary for each partition set
    set_voting: dict[frozenset[str], float] = {frozenset(part): 0.0 for part in partition.parts}
    
    if len(partition) == 4:
        # For 4-set partitions: use reticulation_leaf from profiles with 2 quartets
        # Get all representative partitions (one leaf per set)
        for four_leaf_partition in partition.representative_partitions():
            four_taxa_set = four_leaf_partition.elements
            
            # Get profile for this 4-taxon set
            profile = profileset.get_profile(four_taxa_set)
            if profile is None:
                continue  # Skip if no profile exists
            
            # Check if profile has 2 quartets (cycle) with reticulation_leaf
            if len(profile.quartets) == 2 and profile.reticulation_leaf is not None:
                ret_leaf = profile.reticulation_leaf
                
                # Find which partition set contains the reticulation leaf
                for part in partition.parts:
                    if ret_leaf in part:
                        X = frozenset(part)
                        # Add vote (weighted if requested)
                        if weights:
                            # Use the profile weight (sum of quartet weights)
                            profile_weight = profileset._profiles[four_taxa_set][1]
                            set_voting[X] += profile_weight
                        else:
                            set_voting[X] += 1.0
                        break
    else:
        # For larger partitions: iterate over 4-subpartitions
        for four_sub_partition in partition.subpartitions(4):
            # Get all representative partitions (one leaf per set in the 4-subpartition)
            splits_count = 0.0
            cycles_count = 0.0
            
            for four_leaf_partition in four_sub_partition.representative_partitions():
                four_taxa_set = four_leaf_partition.elements
                
                # Get profile for this 4-taxon set
                profile = profileset.get_profile(four_taxa_set)
                if profile is None:
                    continue  # Skip if no profile exists
                
                # Count splits (1 quartet) vs cycles (2 quartets)
                if weights:
                    profile_weight = profileset._profiles[four_taxa_set][1]
                    if len(profile.quartets) == 1:
                        splits_count += profile_weight
                    elif len(profile.quartets) == 2:
                        cycles_count += profile_weight
                else:
                    if len(profile.quartets) == 1:
                        splits_count += 1.0
                    elif len(profile.quartets) == 2:
                        cycles_count += 1.0
            
            # Calculate cycle percentage
            total = splits_count + cycles_count
            if total > 0:
                cycle_perc = cycles_count / total
                
                # Add cycle_perc to voting for each set in the 4-subpartition
                for part in four_sub_partition.parts:
                    X = frozenset(part)
                    set_voting[X] += cycle_perc
    
    # Sort the sets by their voting score, highest first
    ordered_sets = sorted(set_voting, key=set_voting.get, reverse=True)
    
    return ordered_sets


def _insert_cycle(
    network: 'SemiDirectedPhyNetwork',
    vertex: Any,
    circular_setorder: CircularSetOrdering,
    reticulation_ranking: list[frozenset[str]] | None = None,
    outgroup: str | None = None,
) -> 'SemiDirectedPhyNetwork | MixedPhyNetwork':
    """
    Insert a cycle into the network by replacing a vertex.
    
    This function creates a new network by:
    1. Removing the specified vertex
    2. Creating an undirected cycle of nodes (one per set in the circular ordering)
    3. Reconnecting cut-edges to the appropriate cycle nodes
    4. If reticulation_ranking is provided, iteratively tries to make each set a hybrid
       by converting cycle edges to directed edges, checking that source components = 1
    5. If reticulation_ranking is None, returns a MixedPhyNetwork with no hybrid
    
    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network.
    vertex : Any
        The vertex to replace with a cycle. Must be a cut-vertex in the network.
    circular_setorder : CircularSetOrdering
        The circular ordering of partition sets that determines the cycle structure.
        Must match the partition induced by the vertex.
    reticulation_ranking : list[frozenset[str]] | None, optional
        Ordered list of partition sets to try as hybrid sets (most likely first).
        If None, no hybrid is created and a MixedPhyNetwork is returned.
        By default None.
    
    Returns
    -------
    SemiDirectedPhyNetwork | MixedPhyNetwork
        A new network with the vertex replaced by a cycle. Returns SemiDirectedPhyNetwork
        if a valid hybrid was found, MixedPhyNetwork if reticulation_ranking is None.
    
    Raises
    ------
    ValueError
        If vertex is not a cut-vertex.
        If circular_setorder does not match the partition induced by vertex.
        If no valid hybrid configuration is found (when reticulation_ranking is provided).
    
    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> from phylozoo.core.primitives.circular_ordering import CircularSetOrdering
    >>> 
    >>> # Create a simple network
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 4)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
    ... )
    >>> 
    >>> # Create circular set ordering
    >>> ordering = CircularSetOrdering([{'A'}, {'B'}, {'C'}])
    >>> 
    >>> # Insert cycle with hybrid ranking
    >>> ranking = [frozenset({'A'}), frozenset({'B'}), frozenset({'C'})]
    >>> new_net = _insert_cycle(net, 3, ordering, reticulation_ranking=ranking)
    >>> isinstance(new_net, SemiDirectedPhyNetwork)
    True
    """
    from ...core.network.sdnetwork.derivations import partition_from_blob
    from ...core.network.sdnetwork.features import cut_vertices
    
    # Validate that vertex is a cut-vertex
    cut_verts = cut_vertices(network)
    if vertex not in cut_verts:
        raise ValueError(
            f"Vertex {vertex} is not a cut-vertex in the network"
        )
    
    # Get partition from cut-vertex
    induced_partition, edge_taxa_list = partition_from_blob(
        network, {vertex}, return_edge_taxa=True
    )
    
    # Validate that circular_setorder matches the partition
    # Use the parent Partition.__eq__ method directly (CircularSetOrdering is a Partition)
    from ...core.primitives.partition import Partition
    if not Partition.__eq__(circular_setorder, induced_partition):
        raise ValueError(
            "Circular set ordering does not match partition induced by vertex. "
            f"Partition: {induced_partition.parts}, Ordering: {circular_setorder.parts}"
        )
    
    # Validate reticulation_ranking if provided
    if reticulation_ranking is not None:
        for ret_set in reticulation_ranking:
            ret_frozen = frozenset(ret_set) if not isinstance(ret_set, frozenset) else ret_set
            if ret_frozen not in circular_setorder:
                raise ValueError(
                    f"Reticulation set {ret_set} in ranking is not part of the circular set ordering"
                )
    
    # Create a copy of the graph to modify
    graph_copy = network._graph.copy()
    
    def replace_vertex_with_cycle(
        graph: 'MixedMultiGraph',
        vertex: Any,
        circular_setorder: CircularSetOrdering,
        edge_taxa_list: list[tuple[Any, Any, frozenset[str]]],
        original_graph: 'MixedMultiGraph',
    ) -> tuple['MixedMultiGraph', list[Any], dict[frozenset[str], Any]]:
        """Replace a vertex with an undirected cycle in the graph."""
        # Remove the vertex
        graph.remove_node(vertex)
        
        # Generate new cycle node IDs
        existing_nodes = set(graph.nodes)
        max_node = max(existing_nodes, default=-1) if existing_nodes else -1
        n_cycle_nodes = len(circular_setorder)
        cycle_nodes = list(range(max_node + 1, max_node + 1 + n_cycle_nodes))
        
        # Add cycle nodes to the graph
        for cycle_node in cycle_nodes:
            graph.add_node(cycle_node)
        
        # Create cycle edges (undirected) - connect consecutive nodes in the cycle
        cycle_edges: list[tuple[Any, Any]] = []
        for i in range(n_cycle_nodes):
            u = cycle_nodes[i]
            v = cycle_nodes[(i + 1) % n_cycle_nodes]
            cycle_edges.append((u, v))
        
        # Add cycle edges (undirected)
        for u, v in cycle_edges:
            graph.add_undirected_edge(u, v)
        
        # Map each partition set to its corresponding cycle node
        set_to_cycle_node: dict[frozenset[str], Any] = {}
        for i, part_set in enumerate(circular_setorder.setorder):
            set_to_cycle_node[frozenset(part_set)] = cycle_nodes[i]
        
        # Reconnect edges from edge_taxa_list
        for u, v, taxa_set in edge_taxa_list:
            # v was the vertex, u is the node in the component
            # Find which cycle node corresponds to taxa_set
            cycle_node = set_to_cycle_node.get(taxa_set)
            if cycle_node is None:
                raise ValueError(
                    f"Could not find cycle node for taxa set {taxa_set}"
                )
            
            # Add edge from u to cycle_node (preserve edge attributes if any)
            # Check if original edge was directed or undirected
            attrs: dict[str, Any] = {}
            if original_graph.has_edge(u, vertex):
                # Get edge attributes from undirected or directed edge
                for key in original_graph[u][vertex]:
                    attrs = original_graph[u][vertex][key].copy()
                    break
                # Remove gamma if present (cycle edges don't have gamma)
                attrs.pop('gamma', None)
            elif original_graph.has_edge(vertex, u):
                # Directed edge from vertex to u
                for key in original_graph[vertex][u]:
                    attrs = original_graph[vertex][u][key].copy()
                    break
                # Remove gamma if present
                attrs.pop('gamma', None)
            
            # Add undirected edge with preserved attributes
            graph.add_undirected_edge(u, cycle_node, **attrs)
        
        return graph, cycle_nodes, set_to_cycle_node
    
    def try_hybrid(
        graph: 'MixedMultiGraph',
        ret_set: frozenset[str],
        circular_setorder: CircularSetOrdering,
        cycle_nodes: list[Any],
        outgroup: str | None = None,
    ) -> 'MixedMultiGraph | None':
        """Try to make a hybrid node from a given partition set."""
        from ...core.primitives.m_multigraph.features import source_components
        
        # Find the index of this set in the circular_setorder
        setorder_list = list(circular_setorder.setorder)
        ret_frozen = frozenset(ret_set) if not isinstance(ret_set, frozenset) else ret_set
        
        ret_idx = None
        for i, s in enumerate(setorder_list):
            if frozenset(s) == ret_frozen:
                ret_idx = i
                break
        
        if ret_idx is None:
            return None  # Set not found
        
        # Create a copy of the graph for this attempt
        test_graph = graph.copy()
        
        n_cycle_nodes = len(cycle_nodes)
        
        # Convert the two cycle edges adjacent to ret_idx to directed edges pointing to the hybrid node
        prev_idx = (ret_idx - 1) % n_cycle_nodes
        next_idx = (ret_idx + 1) % n_cycle_nodes
        
        hybrid_node = cycle_nodes[ret_idx]
        parent1_node = cycle_nodes[prev_idx]
        parent2_node = cycle_nodes[next_idx]
        
        # Remove the undirected edges
        test_graph.remove_edge(parent1_node, hybrid_node)
        test_graph.remove_edge(parent2_node, hybrid_node)
        
        # Add directed edges with gamma values
        test_graph.add_directed_edge(parent1_node, hybrid_node, gamma=0.5)
        test_graph.add_directed_edge(parent2_node, hybrid_node, gamma=0.5)
        
        # Check if number of source components is 1
        components = source_components(test_graph)
        if len(components) != 1:
            return None
        
        # If outgroup is specified, check if it's in the source component
        if outgroup is not None:
            nodes_in_component, _, _ = components[0]
            # Get the node corresponding to the outgroup taxon
            outgroup_node = network._label_to_node.get(outgroup)
            if outgroup_node is None or outgroup_node not in nodes_in_component:
                return None
        
        return test_graph
    
    # Replace vertex with undirected cycle
    graph_with_cycle, cycle_nodes, set_to_cycle_node = replace_vertex_with_cycle(
        graph_copy, vertex, circular_setorder, edge_taxa_list, network._graph
    )
    
    # If no reticulation_ranking provided, return MixedPhyNetwork
    if reticulation_ranking is None:
        from ...core.network.sdnetwork.conversions import sdnetwork_from_graph
        return sdnetwork_from_graph(graph_with_cycle, network_type='mixed')
    
    # Try each set in the ranking as a hybrid
    for ret_set in reticulation_ranking:
        ret_frozen = frozenset(ret_set) if not isinstance(ret_set, frozenset) else ret_set
        
        # Try to make this set a hybrid
        hybrid_graph = try_hybrid(
            graph_with_cycle, ret_frozen, circular_setorder, cycle_nodes, outgroup
        )
        
        if hybrid_graph is not None:
            # Valid configuration found
            from ...core.network.sdnetwork.conversions import sdnetwork_from_graph
            return sdnetwork_from_graph(hybrid_graph, network_type='semi-directed')
    
    # No valid hybrid configuration found
    raise ValueError(
        "No valid hybrid configuration found from reticulation_ranking that "
        "maintains exactly one source component"
    )


def resolve_cycles(
    profileset: 'SqQuartetProfileSet',
    tree: 'SemiDirectedPhyNetwork',
    outgroup: str | None = None,
    rho: tuple[float, float, float, float] = (0.5, 1.0, 0.5, 1.0),
    tsp_threshold: int | None = 13,
) -> 'SemiDirectedPhyNetwork':
    """
    Resolve cycles in a tree by replacing high-degree vertices with cycles.
    
    This function converts a tree into a semi-directed level-1 triangle-free
    network by iteratively replacing high-degree vertices with cycles. For each vertex
    with degree > 3 (processed from highest to lowest degree):
    1. Compute the partition induced by that vertex
    2. Compute quartet distances and solve TSP to get circular set ordering
    3. Get hybrid ranking from quartet profiles
    4. Insert a cycle with the hybrid configuration
    
    Parameters
    ----------
    profileset : SqQuartetProfileSet
        The squirrel quartet profile set to use for distance computation.
        Must be dense.
    tree : SemiDirectedPhyNetwork
        The starting tree. Must be a tree.
    outgroup : str | None, optional
        Optional outgroup taxon. If specified, the outgroup must remain in the
        source component when inserting cycles. By default None.
    rho : tuple[float, float, float, float], optional
        Rho vector (rho_c, rho_s, rho_a, rho_o) for distance computation.
        By default (0.5, 1.0, 0.5, 1.0) (Squirrel/MONAD).
    tsp_threshold : int | None, optional
        Maximum partition size for which to solve TSP optimally. If partition size
        is larger, use the approximate method. If None, always use optimal method.
        By default 13.
    
    Returns
    -------
    SemiDirectedPhyNetwork
        A semi-directed level-1 triangle-free network with cycles resolved.
    
    Raises
    ------
    ValueError
        If tree is not a tree.
        If profileset is not dense.
        If no valid network can be constructed.
    
    Examples
    --------
    >>> from phylozoo.inference.squirrel.sqprofileset import SqQuartetProfileSet
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> 
    >>> # Create a profile set and tree
    >>> profileset = SqQuartetProfileSet(...)
    >>> tree = SemiDirectedPhyNetwork(...)
    >>> 
    >>> # Resolve cycles
    >>> network = resolve_cycles(profileset, tree)
    """
    from ...core.network.sdnetwork.derivations import partition_from_blob
    from ...core.network.sdnetwork.features import cut_vertices
    
    # Start with the tree
    network = tree
    
    # Get all cut-vertices (internal nodes) sorted by degree (highest first)
    cut_verts = cut_vertices(network)
    if not cut_verts:
        # No cut-vertices, return tree as-is
        return network
    
    # Sort by degree (highest first)
    nodes_by_degree = sorted(
        cut_verts,
        key=lambda v: network.degree(v),
        reverse=True
    )
    
    # Process nodes with degree > 3
    for vertex in nodes_by_degree:
        if network.degree(vertex) <= 3:
            continue  # Skip nodes with degree <= 3
        
        # Get partition from cut-vertex
        induced_partition, _ = partition_from_blob(
            network, {vertex}, return_edge_taxa=False
        )
        
        # Compute circular set ordering
        tsp_method = 'optimal' if (tsp_threshold is None or len(induced_partition) <= tsp_threshold) else 'simulated_annealing'
        circular_setorder = _qprofiles_to_circular_ordering(
            profileset, induced_partition, rho=rho, tsp_method=tsp_method
        )
        
        # Get hybrid ranking
        hybrid_ranking = _qprofiles_to_hybrid_ranking(
            profileset, induced_partition, weights=True
        )
        
        # Insert cycle
        network = _insert_cycle(
            network, vertex, circular_setorder,
            reticulation_ranking=hybrid_ranking,
            outgroup=outgroup
        )
        
        # Ensure we still have a SemiDirectedPhyNetwork
        if not isinstance(network, SemiDirectedPhyNetwork):
            raise ValueError(
                f"Network became MixedPhyNetwork after inserting cycle at vertex {vertex}"
            )
    
    return network

