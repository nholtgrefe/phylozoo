"""
Cycle resolution module.

This module provides functions for resolving cycles in phylogenetic networks
by computing circular set orderings based on quartet profiles.
"""

from typing import TYPE_CHECKING, Literal

from ...core.distance.operations import approximate_tsp_tour, optimal_tsp_tour
from ...core.primitives.circular_ordering import CircularSetOrdering
from ...core.primitives.partition import Partition
from ...core.quartet.qdistance import quartet_distance_with_partition

if TYPE_CHECKING:
    from ...core.quartet.qprofileset import QuartetProfileSet


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

