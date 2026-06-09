"""
Distance matrix classification module.

This module provides functions for classifying distance matrices based on mathematical
properties: triangle inequality, metric properties (triangle inequality, symmetry,
non-negativity), Kalmanson conditions (circular ordering constraints), and split-
decomposition properties (tree metrics, total decomposability).
"""

from __future__ import annotations

from typing import TypeVar

import numpy as np
from numba import njit

from ...utils.exceptions import PhyloZooValueError
from ..primitives.circular_ordering import CircularOrdering
from .base import DistanceMatrix

T = TypeVar("T")


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

    return _check_triangle_inequality_numba(distance_matrix._matrix, len(distance_matrix))  # type: ignore[no-any-return]


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
        is_nonnegative(distance_matrix)
        and satisfies_triangle_inequality(distance_matrix)
        and has_zero_diagonal(distance_matrix)
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
    return is_nonnegative(distance_matrix) and satisfies_triangle_inequality(distance_matrix)


def is_kalmanson(distance_matrix: DistanceMatrix, circular_order: CircularOrdering[T]) -> bool:
    """
    Check if the distance matrix is Kalmanson with respect to a circular order.

    A distance matrix is Kalmanson with respect to a circular order if it satisfies
    the Kalmanson inequalities for all quadruples of labels in that order.

    The Kalmanson conditions are classical inequalities for circular metrics :cite:`Kalmanson1975`.

    For a circular order (l1, l2, ..., ln), the Kalmanson conditions are:

    - d(ei, ej) + d(ek, el) <= d(ei, ek) + d(ej, el) for all i < j < k < l
    - d(ei, el) + d(ej, ek) <= d(ei, ek) + d(ej, el) for all i < j < k < l

    Parameters
    ----------
    distance_matrix : DistanceMatrix
        The distance matrix to check.
    circular_order : CircularOrdering[T]
        A circular ordering of all labels in the distance matrix. Must contain
        the same elements as the distance matrix labels.

    Returns
    -------
    bool
        True if the matrix is Kalmanson with respect to the given order, False otherwise.

    Raises
    ------
    PhyloZooValueError
        If circular_order is empty, does not contain all labels, or if the matrix is not
        a pseudo-metric.
    TypeError
        If circular_order is not a CircularOrdering.

    Examples
    --------
    >>> import numpy as np
    >>> from phylozoo.core.distance import DistanceMatrix
    >>> from phylozoo.core.distance.classifications import is_kalmanson
    >>> from phylozoo.core.primitives.circular_ordering import CircularOrdering
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
    >>> co = CircularOrdering(['A', 'B', 'C', 'D', 'E'])
    >>> is_kalmanson(dm, co)
    True
    """
    # Input validation
    if not isinstance(circular_order, CircularOrdering):
        raise TypeError(
            f"circular_order must be a CircularOrdering, got {type(circular_order).__name__}"
        )

    if len(circular_order) == 0:
        raise PhyloZooValueError("circular_order cannot be empty")

    # Extract order list from CircularOrdering
    order_list = list(circular_order.order)

    if not set(order_list) == set(distance_matrix.labels):
        raise PhyloZooValueError("circular_order must contain all labels of the distance matrix")

    if len(order_list) != len(distance_matrix.labels):
        raise PhyloZooValueError(
            f"circular_order has {len(order_list)} elements, but distance matrix "
            f"has {len(distance_matrix.labels)} labels"
        )

    if not is_pseudo_metric(distance_matrix):
        raise PhyloZooValueError(
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
            [distance_matrix.get_index(label) for label in order_list], dtype=np.int64
        )
    except PhyloZooValueError as e:
        raise PhyloZooValueError(f"circular_order contains invalid labels: {e}") from e

    return _check_kalmanson_conditions(distance_matrix._matrix, ordered_indices, n)  # type: ignore[no-any-return]


@njit(cache=True)
def _check_kalmanson_conditions(matrix: np.ndarray, ordered_indices: np.ndarray, n: int) -> bool:
    """Numba-accelerated Kalmanson condition check."""
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                for m in range(k + 1, n):
                    ii = ordered_indices[i]
                    jj = ordered_indices[j]
                    kk = ordered_indices[k]
                    ll = ordered_indices[m]

                    if (
                        ii < 0
                        or ii >= n
                        or jj < 0
                        or jj >= n
                        or kk < 0
                        or kk >= n
                        or ll < 0
                        or ll >= n
                    ):
                        return False

                    d_ij = matrix[ii, jj]
                    d_kl = matrix[kk, ll]
                    d_ik = matrix[ii, kk]
                    d_il = matrix[ii, ll]
                    d_jk = matrix[jj, kk]
                    d_jl = matrix[jj, ll]

                    cond1 = d_ij + d_kl - d_ik - d_jl
                    cond2 = d_il + d_jk - d_ik - d_jl

                    if cond1 > 0 or cond2 > 0:
                        return False
    return True


@njit(cache=True)
def _check_four_point_condition(matrix: np.ndarray, n: int, atol: float) -> bool:
    """Numba-accelerated four-point condition check over all C(n,4) quartets."""
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                for ll in range(k + 1, n):
                    s1 = matrix[i, j] + matrix[k, ll]
                    s2 = matrix[i, k] + matrix[j, ll]
                    s3 = matrix[i, ll] + matrix[j, k]
                    max_s = max(s1, max(s2, s3))
                    count = 0
                    if max_s - s1 <= atol:
                        count += 1
                    if max_s - s2 <= atol:
                        count += 1
                    if max_s - s3 <= atol:
                        count += 1
                    if count < 2:
                        return False
    return True


def is_tree_metric(distance_matrix: DistanceMatrix, atol: float = 1e-10) -> bool:
    """
    Check if the distance matrix is a tree metric.

    A distance matrix is a tree metric if and only if it satisfies the
    *four-point condition*: for every choice of four elements i, j, k, l the
    maximum of the three sums

        {d_ij + d_kl,  d_ik + d_jl,  d_il + d_jk}

    is attained by at least two of them :cite:`Bandelt1992`.  Equivalently,
    the split decomposition has zero residual and all d-splits are pairwise
    compatible.

    Parameters
    ----------
    distance_matrix : DistanceMatrix
        The distance matrix to check.
    atol : float, optional
        Absolute tolerance for floating-point comparisons.  By default 1e-10.

    Returns
    -------
    bool
        True if the four-point condition holds for every quartet, False otherwise.

    Examples
    --------
    >>> import numpy as np
    >>> from phylozoo.core.distance import DistanceMatrix
    >>> from phylozoo.core.distance.classifications import is_tree_metric
    >>>
    >>> # Path-tree metric on four taxa (1-2-3-4)
    >>> matrix = np.array([
    ...     [0, 1, 2, 3],
    ...     [1, 0, 1, 2],
    ...     [2, 1, 0, 1],
    ...     [3, 2, 1, 0],
    ... ], dtype=float)
    >>> dm = DistanceMatrix(matrix)
    >>> is_tree_metric(dm)
    True
    >>>
    >>> # Not a tree metric (incompatible splits)
    >>> bad = np.array([[0, 1, 2, 2], [1, 0, 2, 2], [2, 2, 0, 1], [2, 2, 1, 0]], dtype=float)
    >>> dm2 = DistanceMatrix(bad)
    >>> is_tree_metric(dm2)
    False
    """
    n = len(distance_matrix)
    if n < 4:
        return True
    return bool(_check_four_point_condition(distance_matrix._matrix, n, atol))


def is_totally_decomposable(distance_matrix: DistanceMatrix, atol: float = 1e-10) -> bool:
    """
    Check if the distance matrix is totally decomposable.

    A distance matrix is totally decomposable if its split-prime residual d^0
    is zero, i.e., d can be expressed *exactly* as a weighted sum of split metrics::

        d = sum_S alpha_S * delta_S

    This is equivalent to saying that all pairwise distances are fully explained
    by the d-splits (no indecomposable noise remains).

    Parameters
    ----------
    distance_matrix : DistanceMatrix
        The distance matrix to check.
    atol : float, optional
        Absolute tolerance used to compare the residual to zero.
        By default 1e-10.

    Returns
    -------
    bool
        True if the residual is zero within ``atol``, False otherwise.

    See Also
    --------
    split_decomposition : Returns the residual directly.
    is_tree_metric : Stronger condition (totally decomposable + compatible splits).

    Examples
    --------
    >>> import numpy as np
    >>> from phylozoo.core.distance import DistanceMatrix
    >>> from phylozoo.core.distance.classifications import is_totally_decomposable
    >>>
    >>> # Table 1 from :cite:`Bandelt1992` has zero residual
    >>> d = np.array([
    ...     [ 0,  4,  5,  7, 13,  8,  6],
    ...     [ 4,  0,  1,  3,  9, 12, 10],
    ...     [ 5,  1,  0,  2,  8, 13, 11],
    ...     [ 7,  3,  2,  0,  6, 11, 13],
    ...     [13,  9,  8,  6,  0,  5,  7],
    ...     [ 8, 12, 13, 11,  5,  0,  2],
    ...     [ 6, 10, 11, 13,  7,  2,  0],
    ... ], dtype=float)
    >>> dm = DistanceMatrix(d, labels=list("ABCDEFG"))
    >>> is_totally_decomposable(dm)
    True
    """
    from .decomposition import split_decomposition

    _, residual = split_decomposition(distance_matrix)
    return bool(np.allclose(residual.np_array, 0.0, atol=atol))
