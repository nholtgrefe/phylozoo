"""
Tests for the MSA distances module.
"""

import pytest
import numpy as np

from phylozoo.core.sequence.base import MSA
from phylozoo.core.sequence.distances import hamming_distances


class TestHammingDistances:
    """Test cases for Hamming distance computation."""
    
    def test_identical_sequences(self) -> None:
        """Test that identical sequences have zero distance."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGTACGT",
            "taxon3": "ACGTACGT",
        }
        msa = MSA(sequences)
        dm = hamming_distances(msa)
        
        assert dm.get_distance("taxon1", "taxon2") == 0.0
        assert dm.get_distance("taxon1", "taxon3") == 0.0
        assert dm.get_distance("taxon2", "taxon3") == 0.0
    
    def test_different_sequences(self) -> None:
        """Test Hamming distance between different sequences."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGTTTAA",
        }
        msa = MSA(sequences)
        dm = hamming_distances(msa)
        
        # 4 differences out of 8 positions
        assert dm.get_distance("taxon1", "taxon2") == 0.5
    
    def test_symmetric_distance_matrix(self) -> None:
        """Test that distance matrix is symmetric."""
        sequences = {
            "taxon1": "ACGT",
            "taxon2": "TGCA",
            "taxon3": "AAAA",
        }
        msa = MSA(sequences)
        dm = hamming_distances(msa)
        
        assert dm.get_distance("taxon1", "taxon2") == dm.get_distance("taxon2", "taxon1")
        assert dm.get_distance("taxon1", "taxon3") == dm.get_distance("taxon3", "taxon1")
        assert dm.get_distance("taxon2", "taxon3") == dm.get_distance("taxon3", "taxon2")
    
    def test_zero_diagonal(self) -> None:
        """Test that distance from a taxon to itself is zero."""
        sequences = {
            "taxon1": "ACGT",
            "taxon2": "TGCA",
        }
        msa = MSA(sequences)
        dm = hamming_distances(msa)
        
        assert dm.get_distance("taxon1", "taxon1") == 0.0
        assert dm.get_distance("taxon2", "taxon2") == 0.0
    
    def test_gaps_excluded(self) -> None:
        """Test that gaps are excluded from distance computation."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGT----",
        }
        msa = MSA(sequences)
        dm = hamming_distances(msa)
        
        # First 4 positions match, last 4 are gaps (excluded)
        # So distance should be 0.0
        assert dm.get_distance("taxon1", "taxon2") == 0.0
    
    def test_gaps_in_both_sequences(self) -> None:
        """Test that positions with gaps in both sequences are excluded."""
        sequences = {
            "taxon1": "AC-GT",
            "taxon2": "AC-GT",
        }
        msa = MSA(sequences)
        dm = hamming_distances(msa)
        
        # Both have gaps at position 2, which is excluded
        # Remaining positions match, so distance is 0.0
        assert dm.get_distance("taxon1", "taxon2") == 0.0
    
    def test_unknown_characters_excluded(self) -> None:
        """Test that unknown characters (N) are excluded from distance computation."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGTNNNN",
        }
        msa = MSA(sequences)
        dm = hamming_distances(msa)
        
        # First 4 positions match, last 4 are unknown (excluded)
        # So distance should be 0.0
        assert dm.get_distance("taxon1", "taxon2") == 0.0
    
    def test_all_gaps_returns_zero(self) -> None:
        """Test that sequences with only gaps return zero distance."""
        sequences = {
            "taxon1": "------",
            "taxon2": "------",
        }
        msa = MSA(sequences)
        dm = hamming_distances(msa)
        
        # No valid positions to compare
        assert dm.get_distance("taxon1", "taxon2") == 0.0
    
    def test_normalized_distance(self) -> None:
        """Test that distances are normalized to [0, 1] range."""
        sequences = {
            "taxon1": "AAAA",
            "taxon2": "TTTT",
        }
        msa = MSA(sequences)
        dm = hamming_distances(msa)
        
        # All positions differ, so distance should be 1.0
        assert dm.get_distance("taxon1", "taxon2") == 1.0
    
    def test_partial_differences(self) -> None:
        """Test Hamming distance with partial differences."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGTACGA",
        }
        msa = MSA(sequences)
        dm = hamming_distances(msa)
        
        # 1 difference out of 8 positions
        assert dm.get_distance("taxon1", "taxon2") == 0.125

