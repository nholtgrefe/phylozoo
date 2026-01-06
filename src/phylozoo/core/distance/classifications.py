"""
Distance matrix classification module.

This module provides functions for classifying distance matrices (metric, pseudo-metric, Kalmanson, etc.).

TODO: kalmanson should return a circular ordering?
"""

from __future__ import annotations

from typing import TypeVar

import numpy as np
from numba import njit

from .base import DistanceMatrix

T = TypeVar('T')


def satisfies_triangle_inequality(distance_matrix: DistanceMatrix) -> bool:
    """
    Check if the distance matrix satisfies the triangle inequality.
    
    A distance matrix satisfies the triangle inequality if:
    d(i,k) <= d(i,j) + d(j,k) for all i, j, k.
    
    Parameters
    ----------
    distance_matrix : DistanceMatrix
        The distance matrix to check.
    
    Returns
    -------
    bool
        True if triangle inequality holds, False otherwise.
    
    Examples
    --------
    >>> import numpy as np
    >>> from phylozoo.core.distance import DistanceMatrix
    >>> from phylozoo.core.distance.classifications import satisfies_triangle_inequality
    >>> 
    >>> # Matrix satisfying triangle inequality
    >>> matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
    >>> dm = DistanceMatrix(matrix)
    >>> satisfies_triangle_inequality(dm)
    True
    >>> 
    >>> # Matrix violating triangle inequality
    >>> bad_matrix = np.array([[0, 1, 5], [1, 0, 1], [5, 1, 0]])
    >>> bad_dm = DistanceMatrix(bad_matrix)
    >>> satisfies_triangle_inequality(bad_dm)
    False
    """
    @njit(cache=True)
    def _check_triangle_inequality_numba(matrix: np.ndarray, n: int) -> bool:
        """Numba-accelerated triangle inequality check."""
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                for k in range(n):
                    if i == k or j == k:
                        continue
                    # Check triangle inequality: d(i,k) <= d(i,j) + d(j,k)
                    if matrix[i, k] > matrix[i, j] + matrix[j, k]:
                        return False
        return True
    
    return _check_triangle_inequality_numba(distance_matrix._matrix, len(distance_matrix))


def has_zero_diagonal(distance_matrix: DistanceMatrix) -> bool:
    """
    Check if the diagonal of the distance matrix is zero.
    
    Parameters
    ----------
    distance_matrix : DistanceMatrix
        The distance matrix to check.
    
    Returns
    -------
    bool
        True if diagonal is zero, False otherwise.
    
    Examples
    --------
    >>> import numpy as np
    >>> from phylozoo.core.distance import DistanceMatrix
    >>> from phylozoo.core.distance.classifications import has_zero_diagonal
    >>> 
    >>> matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
    >>> dm = DistanceMatrix(matrix)
    >>> has_zero_diagonal(dm)
    True
    """
    return bool(np.all(np.diag(distance_matrix._matrix) == 0))


def is_nonnegative(distance_matrix: DistanceMatrix) -> bool:
    """
    Check if all distances are non-negative.
    
    Parameters
    ----------
    distance_matrix : DistanceMatrix
        The distance matrix to check.
    
    Returns
    -------
    bool
        True if all distances are non-negative, False otherwise.
    
    Examples
    --------
    >>> import numpy as np
    >>> from phylozoo.core.distance import DistanceMatrix
    >>> from phylozoo.core.distance.classifications import is_nonnegative
    >>> 
    >>> matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
    >>> dm = DistanceMatrix(matrix)
    >>> is_nonnegative(dm)
    True
    """
    return bool(np.all(distance_matrix._matrix >= 0))


def is_metric(distance_matrix: DistanceMatrix) -> bool:
    """
    Check if the distance matrix is a metric.
    
    A metric distance matrix satisfies:
    1. Non-negativity: d(x, y) >= 0 for all x, y
    2. Triangle inequality: d(x, z) <= d(x, y) + d(y, z) for all x, y, z
    3. Zero diagonal: d(x, x) = 0 for all x
    4. Symmetry: d(x, y) = d(y, x) for all x, y (already enforced in constructor)
    
    Parameters
    ----------
    distance_matrix : DistanceMatrix
        The distance matrix to check.
    
    Returns
    -------
    bool
        True if the matrix is a metric, False otherwise.
    
    Examples
    --------
    >>> import numpy as np
    >>> from phylozoo.core.distance import DistanceMatrix
    >>> from phylozoo.core.distance.classifications import is_metric
    >>> 
    >>> # Euclidean distance matrix (metric)
    >>> matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
    >>> dm = DistanceMatrix(matrix)
    >>> is_metric(dm)
    True
    >>> 
    >>> # Non-metric (violates triangle inequality)
    >>> bad_matrix = np.array([[0, 1, 5], [1, 0, 1], [5, 1, 0]])
    >>> bad_dm = DistanceMatrix(bad_matrix)
    >>> is_metric(bad_dm)
    False
    """
    return (
        is_nonnegative(distance_matrix) and
        satisfies_triangle_inequality(distance_matrix) and
        has_zero_diagonal(distance_matrix)
    )


def is_pseudo_metric(distance_matrix: DistanceMatrix) -> bool:
    """
    Check if the distance matrix is a pseudo-metric.
    
    A pseudo-metric distance matrix satisfies:
    1. Non-negativity: d(x, y) >= 0 for all x, y
    2. Triangle inequality: d(x, z) <= d(x, y) + d(y, z) for all x, y, z
    
    Note: Unlike a metric, a pseudo-metric does not require d(x, x) = 0
    (though it may still hold).
    
    Parameters
    ----------
    distance_matrix : DistanceMatrix
        The distance matrix to check.
    
    Returns
    -------
    bool
        True if the matrix is a pseudo-metric, False otherwise.
    
    Examples
    --------
    >>> import numpy as np
    >>> from phylozoo.core.distance import DistanceMatrix
    >>> from phylozoo.core.distance.classifications import is_pseudo_metric
    >>> 
    >>> # Pseudo-metric (satisfies non-negativity and triangle inequality)
    >>> matrix = np.array([[0.1, 1, 2], [1, 0.1, 1], [2, 1, 0.1]])
    >>> dm = DistanceMatrix(matrix)
    >>> is_pseudo_metric(dm)
    True
    >>> 
    >>> # Not a pseudo-metric (violates triangle inequality)
    >>> bad_matrix = np.array([[0, 1, 5], [1, 0, 1], [5, 1, 0]])
    >>> bad_dm = DistanceMatrix(bad_matrix)
    >>> is_pseudo_metric(bad_dm)
    False
    """
    return (
        is_nonnegative(distance_matrix) and
        satisfies_triangle_inequality(distance_matrix)
    )


def is_kalmanson(
    distance_matrix: DistanceMatrix,
    circular_order: list[T]
) -> bool:
    """
    Check if the distance matrix is Kalmanson with respect to a circular order.
    
    A distance matrix is Kalmanson with respect to a circular order if it satisfies
    the Kalmanson inequalities for all quadruples of labels in that order.
    
    For a circular order (l1, l2, ..., ln), the Kalmanson conditions are:
    - d(ei, ej) + d(ek, el) <= d(ei, ek) + d(ej, el) for all i < j < k < l
    - d(ei, el) + d(ej, ek) <= d(ei, ek) + d(ej, el) for all i < j < k < l
    
    Parameters
    ----------
    distance_matrix : DistanceMatrix
        The distance matrix to check.
    circular_order : list[T]
        A circular ordering of all labels in the distance matrix.
    
    Returns
    -------
    bool
        True if the matrix is Kalmanson with respect to the given order, False otherwise.
    
    Raises
    ------
    ValueError
        If circular_order does not contain all labels, or if the matrix is not
        a pseudo-metric.
    TypeError
        If circular_order is not a list.
    
    Examples
    --------
    >>> import numpy as np
    >>> from phylozoo.core.distance import DistanceMatrix
    >>> from phylozoo.core.distance.classifications import is_kalmanson
    >>> 
    >>> # Kalmanson matrix (e.g., from a circular network)
    >>> matrix = np.array([
    ...     [0, 1, 2, 2, 1],
    ...     [1, 0, 1, 2, 2],
    ...     [2, 1, 0, 1, 2],
    ...     [2, 2, 1, 0, 1],
    ...     [1, 2, 2, 1, 0]
    ... ])
    >>> dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D', 'E'])
    >>> is_kalmanson(dm, ['A', 'B', 'C', 'D', 'E'])
    True
    """
    # Input validation
    if not isinstance(circular_order, list):
        raise TypeError("circular_order must be a list")
    
    if len(circular_order) == 0:
        raise ValueError("circular_order cannot be empty")
    
    if not set(circular_order) == set(distance_matrix.labels):
        raise ValueError(
            "circular_order must contain all labels of the distance matrix"
        )
    
    if len(circular_order) != len(distance_matrix.labels):
        raise ValueError(
            f"circular_order has {len(circular_order)} elements, but distance matrix "
            f"has {len(distance_matrix.labels)} labels"
        )
    
    if not is_pseudo_metric(distance_matrix):
        raise ValueError(
            "Distance matrix must be pseudo-metric to check Kalmanson property"
        )
    
    # Early exit for small matrices (need at least 4 elements for Kalmanson check)
    n = len(distance_matrix)
    if n < 4:
        # Trivially Kalmanson if less than 4 elements
        return True
    
    # Get ordered indices
    try:
        ordered_indices = np.array(
            [distance_matrix.get_index(label) for label in circular_order],
            dtype=np.int64
        )
    except ValueError as e:
        raise ValueError(
            f"circular_order contains invalid labels: {e}"
        ) from e
    
    @njit(cache=True)
    def _check_kalmanson_conditions(
        matrix: np.ndarray,
        ordered_indices: np.ndarray,
        n: int
    ) -> bool:
        """Numba-accelerated Kalmanson condition check."""
        # Check all combinations of 4 indices (i < j < k < l)
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(j + 1, n):
                    for l in range(k + 1, n):
                        ii = ordered_indices[i]
                        jj = ordered_indices[j]
                        kk = ordered_indices[k]
                        ll = ordered_indices[l]
                        
                        # Bounds checking
                        if ii < 0 or ii >= n or jj < 0 or jj >= n or \
                           kk < 0 or kk >= n or ll < 0 or ll >= n:
                            return False
                        
                        d_ij = matrix[ii, jj]
                        d_kl = matrix[kk, ll]
                        d_ik = matrix[ii, kk]
                        d_il = matrix[ii, ll]
                        d_jk = matrix[jj, kk]
                        d_jl = matrix[jj, ll]
                        
                        # Kalmanson conditions
                        cond1 = d_ij + d_kl - d_ik - d_jl
                        cond2 = d_il + d_jk - d_ik - d_jl
                        
                        if cond1 > 0 or cond2 > 0:
                            return False
        return True
    
    return _check_kalmanson_conditions(distance_matrix._matrix, ordered_indices, n)

