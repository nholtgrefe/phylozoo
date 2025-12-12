"""
Distance matrix operations module.

This module provides operations on distance matrices, including TSP (Traveling Salesman Problem)
solving and other optimization algorithms.
"""

from __future__ import annotations

import itertools
from functools import lru_cache
from typing import TypeVar

import networkx as nx

from .base import DistanceMatrix

T = TypeVar('T')


def optimal_tsp_tour(distance_matrix: DistanceMatrix) -> list[T]:
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
    list[T]
        A circular ordering (tour) of all labels starting and ending at the first label.
        The tour is represented as a list: [l0, l1, l2, ..., ln-1, l0].
    
    Notes
    -----
    This implementation uses the Held-Karp algorithm with dynamic programming and
    memoization. The time complexity is O(n^2 * 2^n), so it's only practical for
    small instances (typically n <= 20).
    
    The algorithm is adapted from the 'python_tsp' package source code (MIT license).
    
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
    >>> set(tour) == set(dm.labels)
    True
    """
    n = len(distance_matrix)
    if n == 0:
        return []
    if n == 1:
        return [distance_matrix.labels[0]]
    
    # Get initial set {1, 2, ..., n-1} as a frozenset (required for @lru_cache)
    N = frozenset(range(1, n))
    memo: dict[tuple[int, frozenset[int]], int] = {}
    matrix = distance_matrix._matrix
    
    # Step 1: Compute minimum distance using dynamic programming
    @lru_cache(maxsize=None)
    def dist(ni: int, N_set: frozenset[int]) -> float:
        """
        Compute minimum distance from node ni through all nodes in N_set back to node 0.
        
        Parameters
        ----------
        ni : int
            Current node index.
        N_set : frozenset[int]
            Set of remaining nodes to visit.
        
        Returns
        -------
        float
            Minimum distance.
        """
        if not N_set:
            return float(matrix[ni, 0])
        
        # Store costs in the form (nj, dist(nj, N_set - {nj}))
        costs = [
            (nj, matrix[ni, nj] + dist(nj, N_set.difference({nj})))
            for nj in N_set
        ]
        nmin, min_cost = min(costs, key=lambda x: x[1])
        memo[(ni, N_set)] = nmin
        
        return min_cost
    
    # Compute optimal distance
    best_distance = dist(0, N)
    
    # Step 2: Reconstruct the path with minimum distance
    ni = 0  # Start at the origin
    solution = [0]
    N_current = N
    
    while N_current:
        ni = memo[(ni, N_current)]
        solution.append(ni)
        N_current = N_current.difference({ni})
    
    # Convert solution indices to labels
    label_tour = [distance_matrix.labels[i] for i in solution]
    
    return label_tour


def approximate_tsp_tour(
    distance_matrix: DistanceMatrix,
    method: str = 'simulated_annealing'
) -> list[T]:
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
    list[T]
        A circular ordering (tour) of all labels. The tour is represented as
        a list: [l0, l1, l2, ..., ln-1].
    
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
        return []
    if n == 1:
        return [distance_matrix.labels[0]]
    
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
    
    return label_tour

