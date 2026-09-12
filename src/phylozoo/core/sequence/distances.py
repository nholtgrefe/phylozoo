"""
Distance computation module for MSA.

This module provides functions for computing distance matrices from MSAs,
including normalized Hamming distances.
"""

from __future__ import annotations

import numpy as np

from ..distance import DistanceMatrix
from .base import MSA


def hamming_distances(msa: MSA) -> DistanceMatrix:
    """
    Compute normalized Hamming distances between all pairs of taxa in an MSA.

    The Hamming distance between two sequences is the number of positions where
    they differ, normalized by the number of positions where both sequences have
    valid nucleotides (i.e., excluding positions with gaps or unknown characters).

    Parameters
    ----------
    msa : MSA
        The multiple sequence alignment.

    Returns
    -------
    DistanceMatrix
        A symmetric distance matrix with normalized Hamming distances between
        all pairs of taxa. The matrix has zero diagonal and is normalized to
        [0, 1] range.

    Examples
    --------
    >>> sequences = {
    ...     "taxon1": "ACGTACGT",
    ...     "taxon2": "ACGTACGT",
    ...     "taxon3": "ACGTTTAA"
    ... }
    >>> msa = MSA(sequences)
    >>> dm = hamming_distances(msa)
    >>> dm.get_distance("taxon1", "taxon2")
    0.0
    >>> dm.get_distance("taxon1", "taxon3")  # 4 differences out of 8 positions
    0.5

    Notes
    -----
    The normalization excludes positions where either sequence has a gap (-) or
    unknown character (N). Only positions where both sequences have valid
    nucleotides (A, C, G, T) are considered.

    Site columns are processed in blocks so the float indicator arrays stay small
    even for very long alignments.

    All pairwise counts are computed as matrix products of 0/1 indicator arrays,
    so the O(n^2 L) work runs in BLAS rather than in a Python loop over pairs.
    """
    # Get coded array (shape: num_taxa, sequence_length)
    coded_array = msa.coded_array
    num_taxa, sequence_length = coded_array.shape

    # For a pair (i, j) the normalised distance is
    #     (both_valid - matches) / both_valid
    # where both_valid counts sites at which both sequences carry a valid code and
    # matches counts the sites at which they carry the *same* valid code. Both counts
    # are inner products of 0/1 indicator rows, so every pair at once is a matrix
    # product: V @ V.T for validity and, per code c, I_c @ I_c.T for identity. This
    # hands the O(n^2 L) work to BLAS instead of a Python loop over pairs.
    # Columns are processed in blocks to bound the size of the indicator arrays. The
    # indicators are float32 for bandwidth; every partial sum inside a block is an
    # integer at most block_columns <= 2**24, so the products are exact. The running
    # totals are float64.
    matches = np.zeros((num_taxa, num_taxa), dtype=np.float64)
    both_valid = np.zeros((num_taxa, num_taxa), dtype=np.float64)
    codes_present = np.unique(coded_array[coded_array >= 0])
    block_columns = max(1024, min(sequence_length, 2**24, 8_000_000 // max(num_taxa, 1)))

    for start in range(0, sequence_length, block_columns):
        block = coded_array[:, start : start + block_columns]
        valid = (block >= 0).astype(np.float32)
        both_valid += valid @ valid.T
        for code in codes_present:
            indicator = (block == code).astype(np.float32)
            matches += indicator @ indicator.T

    # Pairs with no jointly valid site get distance 0, as before.
    with np.errstate(divide="ignore", invalid="ignore"):
        distance_matrix = np.where(both_valid > 0, (both_valid - matches) / both_valid, 0.0)
    # Make the result exactly symmetric with a zero diagonal (BLAS products are
    # symmetric only up to rounding), mirroring the pairwise fill of the loop version.
    distance_matrix = np.triu(distance_matrix, k=1)
    distance_matrix = distance_matrix + distance_matrix.T

    # Create DistanceMatrix with taxa labels
    return DistanceMatrix(distance_matrix, labels=list(msa.taxa_order))
