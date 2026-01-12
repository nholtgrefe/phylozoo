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
    
    This function uses vectorized numpy operations for efficiency, making it
    suitable for large alignments.
    
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
    >>> dm.get_distance("taxon1", "taxon3")
    0.5  # 4 differences out of 8 positions
    
    Notes
    -----
    The normalization excludes positions where either sequence has a gap (-) or
    unknown character (N). Only positions where both sequences have valid
    nucleotides (A, C, G, T) are considered.
    
    This implementation uses vectorized numpy operations for maximum efficiency
    on large alignments.
    """
    # Get coded array (shape: num_taxa, sequence_length)
    coded_array = msa.coded_array
    num_taxa = msa.num_taxa
    
    # Create mask for valid positions (>= 0 means valid nucleotide, < 0 means gap/unknown)
    valid_mask = coded_array >= 0  # Shape: (num_taxa, sequence_length)
    
    # Initialize distance matrix
    distance_matrix = np.zeros((num_taxa, num_taxa), dtype=np.float64)
    
    # Compute distances using vectorized operations
    # For each pair (i, j), we need:
    # 1. Positions where both sequences are valid
    # 2. Count of differences at those positions
    # 3. Normalize by number of valid positions
    
    # Use broadcasting for efficiency: compare all pairs at once
    for i in range(num_taxa):
        seq_i = coded_array[i, :]
        valid_i = valid_mask[i, :]
        
        # Compare with all sequences j > i
        for j in range(i + 1, num_taxa):
            seq_j = coded_array[j, :]
            valid_j = valid_mask[j, :]
            
            # Positions where both sequences have valid nucleotides
            both_valid = valid_i & valid_j
            
            if not np.any(both_valid):
                # No valid positions to compare
                distance_matrix[i, j] = 0.0
                distance_matrix[j, i] = 0.0
                continue
            
            # Count differences at valid positions (vectorized)
            differences = np.sum(seq_i[both_valid] != seq_j[both_valid])
            valid_count = np.sum(both_valid)
            
            # Normalized Hamming distance
            normalized_distance = differences / valid_count if valid_count > 0 else 0.0
            
            distance_matrix[i, j] = normalized_distance
            distance_matrix[j, i] = normalized_distance
    
    # Create DistanceMatrix with taxa labels
    return DistanceMatrix(distance_matrix, labels=list(msa.taxa_order))

