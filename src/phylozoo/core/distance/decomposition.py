"""
Split decomposition module for distance matrices.

Implements the canonical split decomposition of :cite:`Bandelt1992`:

    d = d^0 + sum_S alpha_S * delta_S

where S ranges over all d-splits, alpha_S is the isolation index, delta_S is the
split metric, and d^0 is the split-prime residual.
"""

from __future__ import annotations

import numpy as np
from numba import njit

from ..split.base import Split
from ..split.weighted_splitsystem import WeightedSplitSystem
from ...utils.exceptions import PhyloZooValueError
from .base import DistanceMatrix

_ATOL = 1e-10


@njit(cache=True)
def _compute_isolation_index_nb(
    matrix: np.ndarray,
    idx1: np.ndarray,
    idx2: np.ndarray,
) -> float:
    """
    Numba-accelerated computation of the isolation index of a split.

    Computes::

        (1/2) * min_{i,j in idx1, k,l in idx2}
            (max{d_ij+d_kl, d_ik+d_jl, d_il+d_jk} - d_ij - d_kl)

    For trivial splits (one side has a single element) the single element is used
    as both i and j (or k and l), so the formula remains well-defined.

    Parameters
    ----------
    matrix : np.ndarray
        Square distance matrix (float64).
    idx1 : np.ndarray
        Row/column indices of elements on one side of the split (int64).
    idx2 : np.ndarray
        Row/column indices of elements on the other side of the split (int64).

    Returns
    -------
    float
        The isolation index (>= 0).  Zero means the partition is not a d-split.
    """
    n1 = len(idx1)
    n2 = len(idx2)
    min_val = np.inf

    for a in range(n1):
        i = idx1[a]
        for b in range(a, n1):
            j = idx1[b]
            if n1 > 1 and a == b:
                continue  # require i != j for non-trivial side
            for c in range(n2):
                k = idx2[c]
                for d_idx in range(c, n2):
                    ll = idx2[d_idx]
                    if n2 > 1 and c == d_idx:
                        continue  # require k != ll for non-trivial side
                    dij = matrix[i, j]
                    dkl = matrix[k, ll]
                    dik = matrix[i, k]
                    djl = matrix[j, ll]
                    dil = matrix[i, ll]
                    djk = matrix[j, k]
                    val = max(dij + dkl, max(dik + djl, dil + djk)) - dij - dkl
                    if val < min_val:
                        min_val = val

    if min_val == np.inf:
        return 0.0
    return 0.5 * min_val


def isolation_index(distance_matrix: DistanceMatrix, split: Split) -> float:
    """
    Compute the isolation index of a split with respect to a distance matrix.

    The isolation index of a split (A, B) is defined as::

        alpha_{A,B} = (1/2) * min_{i,j in A, k,l in B}
            (max{d_ij+d_kl, d_ik+d_jl, d_il+d_jk} - d_ij - d_kl)

    The index is always non-negative.  It is strictly positive if and only if
    (A, B) is a *d-split* of ``distance_matrix`` — i.e., for every choice of
    i, j in A and k, l in B the sum d_ij + d_kl is not the largest of the three
    quartet sums (:cite:`Bandelt1992`, Eq. 1-2).

    For trivial splits (:math:`|A| = 1` or :math:`|B| = 1`) the formula reduces to::

        alpha = (1/2) * min_{k != l in B} (d_ak + d_al - d_kl)

    where a is the single element of A (and symmetrically for :math:`|B| = 1`).

    Parameters
    ----------
    distance_matrix : DistanceMatrix
        The distance matrix to evaluate against.
    split : Split
        The split whose isolation index is to be computed.  All elements of the
        split must be present as labels in ``distance_matrix``.

    Returns
    -------
    float
        The isolation index (>= 0).

    Raises
    ------
    PhyloZooValueError
        If any element of ``split`` is not a label in ``distance_matrix``.

    Examples
    --------
    >>> import numpy as np
    >>> from phylozoo.core.distance import DistanceMatrix
    >>> from phylozoo.core.distance.decomposition import isolation_index
    >>> from phylozoo.core.split import Split
    >>>
    >>> # Table 1 from :cite:`Bandelt1992`: split EFG|ABCD has index 6
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
    >>> s = Split({"E", "F", "G"}, {"A", "B", "C", "D"})
    >>> isolation_index(dm, s)
    6.0
    """
    labels_set = set(distance_matrix.labels)
    if not split.elements.issubset(labels_set):
        missing = split.elements - labels_set
        raise PhyloZooValueError(
            f"Split contains elements not present in the distance matrix: {missing}"
        )

    idx1 = np.array(
        [distance_matrix.get_index(lbl) for lbl in split.set1],
        dtype=np.int64,
    )
    idx2 = np.array(
        [distance_matrix.get_index(lbl) for lbl in split.set2],
        dtype=np.int64,
    )
    return float(_compute_isolation_index_nb(distance_matrix.np_array, idx1, idx2))


def split_decomposition(
    distance_matrix: DistanceMatrix,
) -> tuple[WeightedSplitSystem, DistanceMatrix]:
    """
    Compute the canonical split decomposition of a distance matrix.

    Decomposes the distance matrix *d* as::

        d = d^0 + sum_S alpha_S * delta_S

    where S ranges over all *d-splits* (partitions with positive isolation
    index), alpha_S is the isolation index of S, delta_S is the corresponding
    split metric (delta_S(i,j) = 1 iff i and j are on opposite sides of S,
    0 otherwise), and d^0 is the split-prime residual — a metric that admits no
    further splits with positive isolation index.

    The algorithm is the recursive procedure of :cite:`Bandelt1992` (Section
    "Finding the d-Splits"): taxa are added one at a time; at each step the
    existing d-splits of the current subset are extended by placing the new
    taxon on either side, and the new trivial split (all prior taxa | new taxon)
    is also tested.

    Parameters
    ----------
    distance_matrix : DistanceMatrix
        The distance matrix to decompose.

    Returns
    -------
    weighted_system : WeightedSplitSystem
        All d-splits with their isolation indices as weights.  Includes trivial
        splits (:math:`|A| = 1` or :math:`|B| = 1`) if their isolation index is positive.
        Empty if no d-splits exist (only possible for n < 2).
    residual : DistanceMatrix
        The split-prime residual d^0 = d - d^1, where
        d^1 = sum_S alpha_S * delta_S.  Has the same labels as the input.

    Examples
    --------
    >>> import numpy as np
    >>> from phylozoo.core.distance import DistanceMatrix
    >>> from phylozoo.core.distance.decomposition import split_decomposition
    >>>
    >>> # Table 1 from :cite:`Bandelt1992`: four d-splits, zero residual
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
    >>> system, residual = split_decomposition(dm)
    >>> len(system.splits)
    4
    >>> np.allclose(residual.np_array, 0)
    True
    """
    n = len(distance_matrix)
    labels = list(distance_matrix.labels)
    matrix = distance_matrix.np_array

    if n < 2:
        residual = DistanceMatrix(np.zeros((n, n), dtype=np.float64), labels=labels)
        return WeightedSplitSystem(), residual

    # Internal representation: frozensets of matrix row indices.
    # Canonical split key: frozenset({s1, s2}) so that s1|s2 == s2|s1.
    # current_splits maps canonical_key -> (s1, s2, alpha) with s1, s2 frozensets of
    # row indices and alpha the split's isolation index on the current subset. The
    # index is carried along so the final splits need not be re-evaluated.
    current_splits: dict[frozenset, tuple[frozenset, frozenset, float]] = {}

    for i in range(1, n):
        new_splits: dict[frozenset, tuple[frozenset, frozenset, float]] = {}

        # --- Test trivial split {0,...,i-1} | {i} ---
        s1_triv = frozenset(range(i))
        s2_triv = frozenset({i})
        key_triv = frozenset([s1_triv, s2_triv])
        alpha = _compute_isolation_index_nb(
            matrix,
            np.array(sorted(s1_triv), dtype=np.int64),
            np.array(sorted(s2_triv), dtype=np.int64),
        )
        if alpha > _ATOL:
            new_splits[key_triv] = (s1_triv, s2_triv, float(alpha))

        # --- Extend each d-split of the previous subset ---
        for s1, s2, _previous_alpha in current_splits.values():
            # Option A: add i to s1
            ext1 = s1 | frozenset({i})
            key_a = frozenset([ext1, s2])
            if key_a not in new_splits:
                alpha = _compute_isolation_index_nb(
                    matrix,
                    np.array(sorted(ext1), dtype=np.int64),
                    np.array(sorted(s2), dtype=np.int64),
                )
                if alpha > _ATOL:
                    new_splits[key_a] = (ext1, s2, float(alpha))

            # Option B: add i to s2
            ext2 = s2 | frozenset({i})
            key_b = frozenset([s1, ext2])
            if key_b not in new_splits:
                alpha = _compute_isolation_index_nb(
                    matrix,
                    np.array(sorted(s1), dtype=np.int64),
                    np.array(sorted(ext2), dtype=np.int64),
                )
                if alpha > _ATOL:
                    new_splits[key_b] = (s1, ext2, float(alpha))

        current_splits = new_splits

    # Build WeightedSplitSystem from the d-splits of the full set
    weighted_splits: dict[Split, float] = {}
    for s1_idx, s2_idx, alpha in current_splits.values():
        set1_labels = {labels[j] for j in s1_idx}
        set2_labels = {labels[j] for j in s2_idx}
        weighted_splits[Split(set1_labels, set2_labels)] = alpha

    if weighted_splits:
        from ..split.algorithms import distances_from_splitsystem

        system = WeightedSplitSystem(weighted_splits)
        d1_raw = distances_from_splitsystem(system)
        # Reorder d1 from its own (sorted) label order to ours with one fancy index
        # rather than one get_distance call per matrix entry.
        order = np.array([d1_raw.get_index(label) for label in labels], dtype=np.intp)
        d1_matrix = d1_raw.np_array[np.ix_(order, order)]
    else:
        system = WeightedSplitSystem()
        d1_matrix = np.zeros((n, n), dtype=np.float64)

    residual_arr = matrix - d1_matrix
    # Symmetrize and clamp sub-zero floating-point noise
    residual_arr = (residual_arr + residual_arr.T) / 2.0
    np.clip(residual_arr, 0.0, None, out=residual_arr)
    np.fill_diagonal(residual_arr, 0.0)

    return system, DistanceMatrix(residual_arr, labels=labels)
