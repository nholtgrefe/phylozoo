"""
Distance matrix operations module.

This module provides operations on distance matrices, including TSP (Traveling Salesman Problem)
solving and other optimization algorithms.

TODO: optimize TSP
"""

from __future__ import annotations

import itertools
from functools import lru_cache
from typing import TypeVar

import networkx as nx
import numpy as np
from numba import njit

from ..primitives.circular_ordering import CircularOrdering
from .base import DistanceMatrix

T = TypeVar('T')


@njit(cache=True)
def _held_karp_numba(
    matrix: np.ndarray, n: int, memo: np.ndarray, next_node: np.ndarray
) -> float:
    """
    Numba-accelerated Held-Karp algorithm using bitmasks.
    
    Parameters
    ----------
    matrix : np.ndarray
        Distance matrix (n x n).
    n : int
        Number of nodes.
    memo : np.ndarray
        Memoization table for distances (n x 2^(n-1)).
    next_node : np.ndarray
        Table to reconstruct path (n x 2^(n-1)).
    
    Returns
    -------
    float
        Optimal tour distance.
    """
    # Initialize: dist(node, empty_set) = distance from node to 0
    for i in range(n):
        memo[i, 0] = matrix[i, 0]
    
    # Fill memo table using dynamic programming
    # mask represents the set of nodes still to visit (excluding node 0)
    # mask is a bitmask where bit i-1 represents node i (since node 0 is always visited)
    # memo[ni, mask] = min distance from ni through all nodes in mask back to 0
    # So ni should NOT be in mask (we're already at ni)
    max_mask = 1 << (n - 1)  # 2^(n-1)
    
    for mask in range(1, max_mask):
        for ni in range(n):
            if ni == 0:
                continue
            
            # Only compute memo[ni, mask] if ni is NOT in mask
            # (because we're already at ni, so we don't need to visit it again)
            if mask & (1 << (ni - 1)):
                continue
            
            # If mask is empty, go directly to node 0
            if mask == 0:
                min_cost = matrix[ni, 0]
                best_nj = 0
            else:
                # Find minimum over all nodes nj in mask
                min_cost = np.inf
                best_nj = -1
                
                for nj in range(1, n):
                    if nj == ni:
                        continue
                    if not (mask & (1 << (nj - 1))):
                        continue
                    
                    # Remove nj from mask
                    mask_without_nj = mask & ~(1 << (nj - 1))
                    cost = matrix[ni, nj] + memo[nj, mask_without_nj]
                    if cost < min_cost:
                        min_cost = cost
                        best_nj = nj
            
            memo[ni, mask] = min_cost
            next_node[ni, mask] = best_nj
    
    # Find optimal tour starting from node 0
    full_mask = max_mask - 1  # All nodes except 0
    min_cost = np.inf
    best_first = -1
    
    for nj in range(1, n):
        cost = matrix[0, nj] + memo[nj, full_mask & ~(1 << (nj - 1))]
        if cost < min_cost:
            min_cost = cost
            best_first = nj
    
    return min_cost


def optimal_tsp_tour(distance_matrix: DistanceMatrix) -> CircularOrdering:
    """
    Solve TSP to optimality using dynamic programming (Held-Karp algorithm).
    
    This function finds the optimal traveling salesman tour that visits all labels
    exactly once and returns to the starting point, minimizing the total distance.
    
    Parameters
    ----------
    distance_matrix : DistanceMatrix
        The distance matrix to solve TSP for.
    
    Returns
    -------
    CircularOrdering
        A circular ordering (tour) of all labels. The ordering is in canonical form.
    
    Notes
    -----
    This implementation uses the Held-Karp algorithm with dynamic programming,
    optimized with Numba JIT compilation and bitmask-based set operations.
    The time complexity is O(n^2 * 2^n), so it's only practical for small
    instances (typically n <= 20).
    
    The algorithm uses bitmasks to represent sets of nodes, which is more
    efficient than Python sets and compatible with Numba acceleration.
    
    References
    ----------
    .. [1] Held, M., & Karp, R. M. (1962). A dynamic programming approach to sequencing
           problems. Journal of the Society for Industrial and Applied Mathematics,
           10(1), 196-210.
    .. [2] https://en.wikipedia.org/wiki/Held%E2%80%93Karp_algorithm
    
    Examples
    --------
    >>> import numpy as np
    >>> from phylozoo.core.distance import DistanceMatrix
    >>> from phylozoo.core.distance.operations import optimal_tsp_tour
    >>> 
    >>> # Small example
    >>> matrix = np.array([
    ...     [0, 1, 2, 3],
    ...     [1, 0, 1, 2],
    ...     [2, 1, 0, 1],
    ...     [3, 2, 1, 0]
    ... ])
    >>> dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
    >>> tour = optimal_tsp_tour(dm)
    >>> len(tour) == len(dm.labels)
    True
    >>> set(tour.order) == set(dm.labels)
    True
    """
    n = len(distance_matrix)
    if n == 0:
        return CircularOrdering([])
    if n == 1:
        return CircularOrdering([distance_matrix.labels[0]])
    
    # Convert matrix to contiguous float64 array for Numba
    matrix = np.ascontiguousarray(distance_matrix._matrix, dtype=np.float64)
    
    # Initialize memoization tables
    max_mask = 1 << (n - 1)  # 2^(n-1)
    memo = np.full((n, max_mask), np.inf, dtype=np.float64)
    next_node = np.full((n, max_mask), -1, dtype=np.int32)
    
    # Compute optimal distance using Numba-accelerated function
    _held_karp_numba(matrix, n, memo, next_node)
    
    # Reconstruct the path
    solution = [0]
    full_mask = max_mask - 1  # All nodes except 0
    
    # Find first node after 0
    min_cost = np.inf
    best_first = -1
    for nj in range(1, n):
        mask_without_nj = full_mask & ~(1 << (nj - 1))
        cost = matrix[0, nj] + memo[nj, mask_without_nj]
        if cost < min_cost:
            min_cost = cost
            best_first = nj
    
    # Reconstruct path starting from best_first
    ni = best_first
    current_mask = full_mask & ~(1 << (ni - 1))
    solution.append(ni)
    
    # Reconstruct rest of path using next_node table
    while current_mask != 0:
        next_ni = int(next_node[ni, current_mask])
        if next_ni == 0:
            # We've visited all nodes, return to start
            break
        solution.append(next_ni)
        ni = next_ni
        current_mask = current_mask & ~(1 << (ni - 1))
    
    # Convert solution indices to labels
    label_tour = [distance_matrix.labels[i] for i in solution]
    
    # Return as CircularOrdering (automatically canonicalized)
    return CircularOrdering(label_tour)


def approximate_tsp_tour(
    distance_matrix: DistanceMatrix,
    method: str = 'simulated_annealing'
) -> CircularOrdering:
    """
    Find an approximate TSP tour using heuristic methods.
    
    This function uses approximation algorithms to find a good (but not necessarily
    optimal) traveling salesman tour. These methods are much faster than the optimal
    algorithm and can handle larger instances.
    
    Parameters
    ----------
    distance_matrix : DistanceMatrix
        The distance matrix to solve TSP for.
    method : str, optional
        Heuristic method to use. Must be one of:
        - 'simulated_annealing': Simulated annealing heuristic (default)
        - 'greedy': Greedy nearest-neighbor heuristic
        - 'christofides': Christofides algorithm (for metric distances)
        By default 'simulated_annealing'.
    
    Returns
    -------
    CircularOrdering
        A circular ordering (tour) of all labels. The ordering is in canonical form.
    
    Raises
    ------
    ValueError
        If method is not one of the supported methods.
    
    Notes
    -----
    - **simulated_annealing**: Uses simulated annealing with a greedy initialization.
      Generally produces good solutions but is slower than greedy.
    - **greedy**: Simple nearest-neighbor heuristic. Fast but may produce poor solutions.
    - **christofides**: Provides a 3/2-approximation for metric distances. Slower than
      greedy but guarantees better worst-case performance.
    
    Examples
    --------
    >>> import numpy as np
    >>> from phylozoo.core.distance import DistanceMatrix
    >>> from phylozoo.core.distance.operations import approximate_tsp_tour
    >>> 
    >>> matrix = np.array([
    ...     [0, 1, 2, 3],
    ...     [1, 0, 1, 2],
    ...     [2, 1, 0, 1],
    ...     [3, 2, 1, 0]
    ... ])
    >>> dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
    >>> 
    >>> # Use simulated annealing
    >>> tour1 = approximate_tsp_tour(dm, method='simulated_annealing')
    >>> len(tour1) == len(dm.labels)
    True
    >>> 
    >>> # Use greedy heuristic
    >>> tour2 = approximate_tsp_tour(dm, method='greedy')
    >>> len(tour2) == len(dm.labels)
    True
    """
    if method not in ['simulated_annealing', 'greedy', 'christofides']:
        raise ValueError(
            f"Method must be one of ['simulated_annealing', 'greedy', 'christofides'], "
            f"got '{method}'"
        )
    
    n = len(distance_matrix)
    if n == 0:
        return CircularOrdering([])
    if n == 1:
        return CircularOrdering([distance_matrix.labels[0]])
    
    # Build complete graph with weights
    complete_graph = nx.Graph()
    complete_graph.add_nodes_from(distance_matrix.indices)
    
    # Add weighted edges
    weight_edges = [
        (i, j, distance_matrix._matrix[i, j])
        for i, j in itertools.combinations(distance_matrix.indices, 2)
    ]
    complete_graph.add_weighted_edges_from(weight_edges)
    
    # Solve TSP using NetworkX approximation
    tsp_func = nx.approximation.traveling_salesman_problem
    
    if method == "simulated_annealing":
        method_func = nx.approximation.simulated_annealing_tsp
        tour = tsp_func(
            complete_graph,
            cycle=False,
            method=method_func,
            init_cycle='greedy',
            seed=0
        )
    elif method == "greedy":
        method_func = nx.approximation.greedy_tsp
        tour = tsp_func(complete_graph, cycle=False, method=method_func)
    elif method == "christofides":
        method_func = nx.approximation.christofides
        tour = tsp_func(complete_graph, cycle=False, method=method_func)
    
    # Convert indices to labels
    label_tour = [distance_matrix.labels[i] for i in tour]
    
    # Return as CircularOrdering (automatically canonicalized)
    return CircularOrdering(label_tour)

