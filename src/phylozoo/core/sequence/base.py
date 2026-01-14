"""
Multiple Sequence Alignment (MSA) base module.

This module provides the MSA class for working with multiple sequence alignments.
Sequences are stored internally as numpy arrays for efficient computation.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ...utils.exceptions import PhyloZooValueError
from ...utils.io import IOMixin

# Default nucleotide encoding
DEFAULT_NUCLEOTIDE_ENCODING = {
    "A": 0,
    "a": 0,
    "C": 1,
    "c": 1,
    "G": 2,
    "g": 2,
    "T": 3,
    "t": 3,
    "-": -1,  # Gap
    "N": -1,  # Unknown
    "n": -1,
}


class MSA(IOMixin):
    """
    Immutable Multiple Sequence Alignment.
    
    A Multiple Sequence Alignment (MSA) represents a collection of aligned sequences
    where all sequences have the same length. Each sequence is associated with a
    taxon (sequence identifier).
    
    Sequences are stored internally as numpy arrays (int8) for efficient computation.
    This class is immutable - once created, the alignment structure cannot be
    modified. All sequences are stored internally in a canonical order.
    
    Parameters
    ----------
    sequences : dict[str, str]
        Dictionary mapping taxon names (sequence identifiers) to sequence strings.
        All sequences must have the same length.
    
    Attributes
    ----------
    taxa : frozenset[str]
        Frozen set of all taxon names (read-only).
    taxa_order : tuple[str, ...]
        Tuple of taxon names in canonical order (read-only).
    sequence_length : int
        Length of all sequences in the alignment (read-only).
    num_taxa : int
        Number of taxa in the alignment (read-only).
    
    Raises
    ------
    ValueError
        If sequences dictionary is empty, or if sequences have different lengths.
    TypeError
        If sequences is not a dictionary.
    
    Examples
    --------
    >>> sequences = {
    ...     "taxon1": "ACGTACGT",
    ...     "taxon2": "ACGTACGT",
    ...     "taxon3": "ACGTACGT"
    ... }
    >>> msa = MSA(sequences)
    >>> len(msa)
    3
    >>> msa.sequence_length
    8
    >>> "taxon1" in msa
    True
    >>> msa.get_sequence("taxon1")
    'ACGTACGT'
    >>> msa.get_sequence("nonexistent")
    None
    
    Notes
    -----
    The internal representation uses:
    - Numpy array (`_coded_array`) of shape (num_taxa, sequence_length) with dtype int8
    - Precomputed lookup table for efficient character encoding
    - Reverse lookup table for decoding back to strings when needed
    - Taxa order for consistent indexing
    """
    
    # I/O format configuration
    _default_format = 'fasta'
    _supported_formats = ['fasta', 'nexus']
    
    __slots__ = (
        "_coded_array",
        "_taxa",
        "_taxa_order",
        "_sequence_length",
        "_num_taxa",
        "_character_encoding",
        "_lookup_table",
        "_reverse_lookup",
        "_taxa_to_index",
    )
    
    def __init__(
        self,
        sequences: dict[str, str],
    ) -> None:
        """
        Initialize an MSA object.
        
        Parameters
        ----------
        sequences : dict[str, str]
            Dictionary mapping taxon names to sequence strings.
        """
        # Input validation
        if not isinstance(sequences, dict):
            raise TypeError("sequences must be a dictionary")
        
        if len(sequences) == 0:
            raise PhyloZooValueError("sequences dictionary cannot be empty")
        
        # Validate all values are strings
        for taxon, seq in sequences.items():
            if not isinstance(seq, str):
                raise TypeError(f"Sequence for taxon '{taxon}' must be a string")
            if not isinstance(taxon, str):
                raise TypeError("All taxon names must be strings")
        
        # Check sequence lengths
        lengths = [len(seq) for seq in sequences.values()]
        if len(set(lengths)) != 1:
            raise PhyloZooValueError(
                f"All sequences must have the same length. "
                f"Found lengths: {set(lengths)}"
            )
        
        # Store sequences in canonical order (sorted by taxon name) for fast iteration
        self._taxa_order: tuple[str, ...] = tuple(sorted(sequences.keys()))
        self._taxa: frozenset[str] = frozenset(self._taxa_order)
        self._sequence_length: int = lengths[0]
        self._num_taxa: int = len(self._taxa_order)
        
        # Build taxa to index mapping for O(1) lookup
        self._taxa_to_index: dict[str, int] = {
            taxon: idx for idx, taxon in enumerate(self._taxa_order)
        }
        
        # Character encoding (nucleotide only for now)
        self._character_encoding: dict[str, int] = DEFAULT_NUCLEOTIDE_ENCODING.copy()
        
        # Build lookup table for efficient encoding (256 entries for ASCII)
        self._lookup_table: np.ndarray = np.full(256, -1, dtype=np.int8)
        for char, code in self._character_encoding.items():
            if len(char) == 1:
                self._lookup_table[ord(char)] = code
        
        # Build reverse lookup for decoding (prefer uppercase)
        self._reverse_lookup: dict[int, str] = {}
        for char, code in self._character_encoding.items():
            if code >= 0 and code not in self._reverse_lookup:
                self._reverse_lookup[code] = char.upper()
        self._reverse_lookup[-1] = '-'  # Gap/unknown
        
        # Convert sequences to numpy array immediately (primary storage)
        self._coded_array: np.ndarray = self._build_coded_array(sequences)
    
    def _build_coded_array(self, sequences: dict[str, str]) -> np.ndarray:
        """
        Build the coded array representation of the MSA.
        
        Parameters
        ----------
        sequences : dict[str, str]
            Dictionary mapping taxon names to sequence strings.
        
        Returns
        -------
        np.ndarray
            Array of shape (num_taxa, sequence_length) with dtype int8.
        """
        coded_sequences = []
        for taxon in self._taxa_order:
            seq = sequences[taxon]
            # Encode sequence using lookup table (vectorized)
            coded = self._lookup_table[
                np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
            ]
            coded_sequences.append(coded)
        
        # Stack into single array: shape (num_taxa, sequence_length)
        return np.vstack(coded_sequences)
    
    @property
    def taxa(self) -> frozenset[str]:
        """
        Get the set of all taxon names (read-only).
        
        Returns
        -------
        frozenset[str]
            Frozen set of all taxon names.
        """
        return self._taxa
    
    @property
    def taxa_order(self) -> tuple[str, ...]:
        """
        Get the taxon names in canonical order (read-only).
        
        Returns
        -------
        tuple[str, ...]
            Tuple of taxon names in sorted order.
        """
        return self._taxa_order
    
    @property
    def sequence_length(self) -> int:
        """
        Get the length of sequences in the alignment (read-only).
        
        Returns
        -------
        int
            Length of all sequences.
        """
        return self._sequence_length
    
    @property
    def num_taxa(self) -> int:
        """
        Get the number of taxa in the alignment (read-only).
        
        Returns
        -------
        int
            Number of taxa.
        """
        return self._num_taxa
    
    def __len__(self) -> int:
        """
        Return the number of taxa in the alignment.
        
        Returns
        -------
        int
            Number of taxa.
        """
        return self._num_taxa
    
    def __contains__(self, taxon: str) -> bool:
        """
        Check if a taxon is in the alignment.
        
        Parameters
        ----------
        taxon : str
            Taxon name to check.
        
        Returns
        -------
        bool
            True if taxon is in the alignment, False otherwise.
        
        Examples
        --------
        >>> sequences = {"taxon1": "ACGT", "taxon2": "ACGT"}
        >>> msa = MSA(sequences)
        >>> "taxon1" in msa
        True
        >>> "taxon3" in msa
        False
        """
        return taxon in self._taxa
    
    def __repr__(self) -> str:
        """
        Return string representation of the MSA.
        
        Returns
        -------
        str
            String representation showing number of taxa and sequence length.
        """
        return (
            f"MSA(num_taxa={self._num_taxa}, "
            f"sequence_length={self._sequence_length})"
        )
    
    def __str__(self) -> str:
        """
        Return detailed string representation of the MSA.
        
        Returns
        -------
        str
            Detailed string representation with sequence previews.
        """
        preview_len = min(40, self._sequence_length)
        preview = "..."
        if self._sequence_length <= 40:
            preview = ""
        
        lines = [f"MSA with {self._num_taxa} taxa, sequence length {self._sequence_length}:"]
        for taxon in self._taxa_order:
            seq = self.get_sequence(taxon)
            if seq:
                lines.append(f"  {taxon:<30} {seq[:preview_len]}{preview}")
        
        return "\n".join(lines)
    
    def get_sequence(self, taxon: str) -> str | None:
        """
        Get the sequence for a given taxon (lazy decoding from numpy array).
        
        Parameters
        ----------
        taxon : str
            Taxon name.
        
        Returns
        -------
        str | None
            Sequence string if taxon exists, None otherwise.
        
        Examples
        --------
        >>> sequences = {"taxon1": "ACGT", "taxon2": "ACGT"}
        >>> msa = MSA(sequences)
        >>> msa.get_sequence("taxon1")
        'ACGT'
        >>> msa.get_sequence("taxon3")
        None
        """
        if taxon not in self._taxa_to_index:
            return None
        
        idx = self._taxa_to_index[taxon]
        coded_seq = self._coded_array[idx, :]
        
        # Decode using reverse lookup
        return ''.join(self._reverse_lookup.get(int(code), 'N') for code in coded_seq)
    
    def get_sequences(self) -> dict[str, str]:
        """
        Get a copy of all sequences as a dictionary.
        
        Returns
        -------
        dict[str, str]
            Dictionary mapping taxon names to sequences.
        
        Notes
        -----
        Returns a copy to maintain immutability. The returned dictionary
        can be modified without affecting the MSA.
        """
        result = {}
        for taxon in self._taxa_order:
            result[taxon] = self.get_sequence(taxon) or ''
        return result
    
    @property
    def coded_array(self) -> np.ndarray:
        """
        Get the coded array representation of the MSA.
        
        Returns a numpy array of shape (num_taxa, sequence_length) where each
        character is encoded according to the character encoding dictionary.
        Characters not in the encoding are mapped to -1.
        
        Returns
        -------
        np.ndarray
            Array of shape (num_taxa, sequence_length) with dtype int8.
            Each element is the encoded value of the corresponding character.
        
        Notes
        -----
        This returns the internal array directly (not a copy for efficiency).
        The array is read-only from the user's perspective due to immutability.
        """
        return self._coded_array
    
    def copy(self) -> MSA:
        """
        Create a copy of the MSA.
        
        Returns
        -------
        MSA
            A new MSA instance with the same sequences.
        
        Notes
        -----
        Returns a new instance with the same data. The copy shares no mutable
        state with the original.
        """
        return MSA(sequences=self.get_sequences())

