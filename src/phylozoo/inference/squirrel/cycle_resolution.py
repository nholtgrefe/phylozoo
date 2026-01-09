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
    from ...core.network.sdnetwork import SemiDirectedPhyNetwork
    from .sqprofileset import SqQuartetProfileSet


def resolve_cycle_ordering(
    profileset: 'QuartetProfileSet',
    partition: Partition,
    rho: tuple[float, float, float, float] = (0.5, 1.0, 0.5, 1.0),
    tsp_method: Literal['optimal', 'simulated_annealing', 'greedy', 'christofides'] = 'optimal',
) -> 'CircularSetOrdering':
    """
    Resolve a cycle by computing the optimal circular set ordering.
    
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
    >>> # Resolve cycle ordering
    >>> ordering = resolve_cycle_ordering(profileset, partition)
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


def _reticulation_order(
    profileset: 'SqQuartetProfileSet',
    partition: Partition,
    weights: bool = True,
) -> list[frozenset[str]]:
    """
    Return an ordered list of partition sets by likelihood of being the reticulation set.
    
    The first set in the returned list is the most likely to be the reticulation set,
    the second is the second most likely, etc. This ordering is based on voting from
    quartet profiles:
    - For 4-set partitions: uses reticulation_leaf from 4-cycle profiles
    - For larger partitions: aggregates cycle percentages from 4-subpartitions
    
    Parameters
    ----------
    profileset : SqQuartetProfileSet
        The squirrel quartet profile set to analyze.
    partition : Partition
        A partition of the taxa. The ordering will be computed for the sets
        (parts) of this partition.
    weights : bool, optional
        Whether to use quartet weights in the voting. By default True.
    
    Returns
    -------
    list[frozenset[str]]
        Ordered list of partition sets (frozensets), sorted by reticulation likelihood
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
    >>> # Get reticulation order
    >>> order = _reticulation_order(profileset, partition)
    >>> len(order)
    4
    """
    # TODO: Implement _reticulation_order
    raise NotImplementedError("_reticulation_order is not yet implemented")


def cutvertex_to_cycle(
    network: 'SemiDirectedPhyNetwork',
    vertex: Any,
    circular_setorder: CircularSetOrdering,
    reticulation_set: frozenset[str],
) -> None:
    """
    Replace a cut-vertex with a cycle in the network.
    
    This function modifies the network in-place by:
    1. Removing the specified cut-vertex
    2. Creating a cycle of nodes (one per set in the circular ordering)
    3. Adding cycle edges (undirected) and reticulation edges (directed)
    4. Reconnecting cut-edges to the appropriate cycle nodes
    
    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The semi-directed phylogenetic network to modify (modified in-place).
    vertex : Any
        The cut-vertex to replace with a cycle. Must be a cut-vertex in the network.
    circular_setorder : CircularSetOrdering
        The circular ordering of partition sets that determines the cycle structure.
        Must match the partition induced by the cut-vertex.
    reticulation_set : frozenset[str]
        The partition set that will be below the reticulation. Must be one of the
        sets in circular_setorder.
    
    Raises
    ------
    ValueError
        If vertex is not a cut-vertex.
        If circular_setorder does not match the partition induced by vertex.
        If reticulation_set is not in circular_setorder.
    
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
    >>> # Replace cut-vertex with cycle
    >>> cutvertex_to_cycle(net, 3, ordering, frozenset({'A'}))
    """
    # TODO: Implement cutvertex_to_cycle
    raise NotImplementedError("cutvertex_to_cycle is not yet implemented")


def reconstruct_network_from_tree(
    profileset: 'SqQuartetProfileSet',
    blobtree: 'SemiDirectedPhyNetwork | None' = None,
    outgroup: str | None = None,
    rho: tuple[float, float, float, float] = (0.5, 1.0, 0.5, 1.0),
    tsp_method: Literal['optimal', 'simulated_annealing', 'greedy', 'christofides'] = 'optimal',
    tsp_threshold: int | None = 13,
) -> 'SemiDirectedPhyNetwork':
    """
    Reconstruct a level-1 network from a tree by adding cycles.
    
    This function converts a tree (blob tree) into a semi-directed level-1 triangle-free
    network by iteratively replacing high-degree vertices with cycles. For each vertex
    with degree > 3:
    1. Compute the partition induced by that vertex
    2. Compute quartet distances and solve TSP to get circular set ordering
    3. Try different reticulation sets in order of likelihood
    4. Replace the vertex with a cycle that makes the network rootable
    
    Parameters
    ----------
    profileset : SqQuartetProfileSet
        The squirrel quartet profile set to use for distance computation.
        Must be dense.
    blobtree : SemiDirectedPhyNetwork | None, optional
        The starting tree. If None, the T*-tree is used. Must be a tree.
        By default None.
    outgroup : str | None, optional
        If specified, the network must be rootable at this outgroup taxon.
        By default None.
    rho : tuple[float, float, float, float], optional
        Rho vector (rho_c, rho_s, rho_a, rho_o) for distance computation.
        By default (0.5, 1.0, 0.5, 1.0) (Squirrel/MONAD).
    tsp_method : Literal['optimal', 'simulated_annealing', 'greedy', 'christofides'], optional
        TSP method to use for finding circular orderings. By default 'optimal'.
    tsp_threshold : int | None, optional
        Maximum partition size for which to solve TSP optimally. If partition size
        is larger, use the approximate method. If None, always use optimal method.
        By default 13.
    
    Returns
    -------
    SemiDirectedPhyNetwork
        A semi-directed level-1 triangle-free network reconstructed from the tree.
    
    Raises
    ------
    ValueError
        If blobtree is not a tree.
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
    >>> # Reconstruct network
    >>> network = reconstruct_network_from_tree(profileset, blobtree=tree)
    """
    # TODO: Implement reconstruct_network_from_tree
    raise NotImplementedError("reconstruct_network_from_tree is not yet implemented")

