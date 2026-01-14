"""
Bootstrap module for Multiple Sequence Alignments.

This module provides functions for generating bootstrap replicates of MSAs,
including standard bootstrap and per-gene bootstrap methods.
"""

from __future__ import annotations

from typing import Iterator

import numpy as np

from .base import MSA
from ...utils.exceptions import PhyloZooValueError


def bootstrap(
    msa: MSA,
    n_bootstrap: int,
    length: int | None = None,
    seed: int | None = None,
) -> Iterator[MSA]:
    """
    Generate bootstrap replicates of an MSA.
    
    This function yields bootstrapped versions of the input MSA by sampling
    columns with replacement. Each bootstrap replicate has the same number
    of taxa and the same sequence length (or specified length).
    
    Parameters
    ----------
    msa : MSA
        The multiple sequence alignment to bootstrap.
    n_bootstrap : int
        Number of bootstrap replicates to generate. Must be positive.
    length : int | None, optional
        Length of the bootstrapped sequences. If None, uses the original
        sequence length. By default None.
    seed : int | None, optional
        Random seed for reproducibility. If None, uses random state.
        By default None.
    
    Yields
    ------
    MSA
        A bootstrapped MSA replicate.
    
    Examples
    --------
    >>> sequences = {
    ...     "taxon1": "ACGTACGT",
    ...     "taxon2": "ACGTACGT",
    ...     "taxon3": "ACGTACGT"
    ... }
    >>> msa = MSA(sequences)
    >>> bootstrapped = next(bootstrap(msa, n_bootstrap=1, seed=42))
    >>> bootstrapped.sequence_length == msa.sequence_length
    True
    >>> bootstrapped.num_taxa == msa.num_taxa
    True
    
    Notes
    -----
    Bootstrap sampling is done by sampling column indices with replacement
    using numpy's random number generator for efficiency. The function uses
    vectorized operations to create bootstrapped sequences from the coded
    array representation. Uses MSA.from_coded_array() to avoid encoding/decoding
    overhead for maximum efficiency.
    """
    if length is None:
        length = msa.sequence_length
    
    if length <= 0:
        raise PhyloZooValueError(f"Bootstrap length must be positive, got {length}")
    
    if n_bootstrap <= 0:
        raise PhyloZooValueError(f"n_bootstrap must be positive, got {n_bootstrap}")
    
    # Set up random number generator
    rng = np.random.default_rng(seed)
    
    # Get coded array (shape: num_taxa, sequence_length)
    coded_array = msa.coded_array
    
    # Efficiency: Create MSA objects directly from coded arrays without encoding/decoding
    # We use __new__ to bypass __init__ and manually set attributes
    for _ in range(n_bootstrap):
        # Sample column indices with replacement
        column_indices = rng.integers(0, msa.sequence_length, size=length)
        
        # Extract bootstrapped columns using advanced indexing
        # Shape: (num_taxa, length)
        bootstrapped_array = coded_array[:, column_indices]
        
        # Efficiency optimization: Create MSA directly from coded array
        # This avoids encoding/decoding cycles
        bootstrapped_msa = MSA.from_coded_array(
            coded_array=bootstrapped_array,
            taxa_order=msa.taxa_order,
        )
        
        yield bootstrapped_msa


def bootstrap_per_gene(
    msa: MSA,
    gene_lengths: list[int],
    n_bootstrap: int,
    seed: int | None = None,
) -> Iterator[MSA]:
    """
    Generate bootstrap replicates of an MSA with per-gene bootstrapping.
    
    This function bootstraps each gene separately, maintaining the structure
    of gene boundaries. For each gene, columns are sampled with replacement
    only from within that gene's boundaries.
    
    Parameters
    ----------
    msa : MSA
        The multiple sequence alignment to bootstrap.
    gene_lengths : list[int]
        List of gene lengths. The sum of gene lengths must equal the
        sequence length of the MSA.
    n_bootstrap : int
        Number of bootstrap replicates to generate. Must be positive.
    seed : int | None, optional
        Random seed for reproducibility. If None, uses random state.
        By default None.
    
    Yields
    ------
    MSA
        A bootstrapped MSA replicate with per-gene bootstrapping.
    
    Raises
    ------
    PhyloZooValueError
        If gene_lengths is empty, contains non-positive values, or if
        the sum of gene lengths does not equal the MSA sequence length.
    
    Examples
    --------
    >>> sequences = {
    ...     "taxon1": "ACGTACGT",
    ...     "taxon2": "ACGTACGT",
    ...     "taxon3": "ACGTACGT"
    ... }
    >>> msa = MSA(sequences)
    >>> gene_lengths = [4, 4]  # Two genes of length 4 each
    >>> bootstrapped = next(bootstrap_per_gene(msa, gene_lengths, n_bootstrap=1, seed=42))
    >>> bootstrapped.sequence_length == msa.sequence_length
    True
    
    Notes
    -----
    Per-gene bootstrapping preserves gene boundaries by sampling columns
    only within each gene's region. This is useful when genes have different
    evolutionary histories or when gene boundaries are biologically meaningful.
    Uses MSA.from_coded_array() to avoid encoding/decoding overhead for maximum
    efficiency.
    """
    if not gene_lengths:
        raise PhyloZooValueError("gene_lengths cannot be empty")
    
    if any(length <= 0 for length in gene_lengths):
        raise PhyloZooValueError(
            f"All gene lengths must be positive, got {gene_lengths}"
        )
    
    total_length = sum(gene_lengths)
    if total_length != msa.sequence_length:
        raise PhyloZooValueError(
            f"Sum of gene lengths ({total_length}) must equal MSA sequence length "
            f"({msa.sequence_length})"
        )
    
    if n_bootstrap <= 0:
        raise PhyloZooValueError(f"n_bootstrap must be positive, got {n_bootstrap}")
    
    # Set up random number generator
    rng = np.random.default_rng(seed)
    
    # Get coded array (shape: num_taxa, sequence_length)
    coded_array = msa.coded_array
    
    # Build gene boundaries
    gene_boundaries: list[tuple[int, int]] = []
    pos = 0
    for gene_len in gene_lengths:
        start = pos
        end = pos + gene_len
        gene_boundaries.append((start, end))
        pos = end
    
    # Efficiency: Create MSA objects directly from coded arrays without encoding/decoding
    for _ in range(n_bootstrap):
        # Sample columns per gene
        column_indices_list: list[np.ndarray] = []
        
        for start, end in gene_boundaries:
            gene_len = end - start
            # Sample indices within this gene's boundaries
            gene_indices = rng.integers(start, end, size=gene_len)
            column_indices_list.append(gene_indices)
        
        # Concatenate all gene indices
        column_indices = np.concatenate(column_indices_list)
        
        # Extract bootstrapped columns using advanced indexing
        # Shape: (num_taxa, sequence_length)
        bootstrapped_array = coded_array[:, column_indices]
        
        # Efficiency optimization: Create MSA directly from coded array
        # This avoids encoding/decoding cycles
        bootstrapped_msa = MSA.from_coded_array(
            coded_array=bootstrapped_array,
            taxa_order=msa.taxa_order,
        )
        
        yield bootstrapped_msa
