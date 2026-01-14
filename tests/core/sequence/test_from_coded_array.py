"""
Tests for MSA.from_coded_array class method.
"""

from __future__ import annotations

import pytest
import numpy as np

from phylozoo.core.sequence import MSA
from phylozoo.core.sequence.base import DEFAULT_NUCLEOTIDE_ENCODING
from phylozoo.utils.exceptions import PhyloZooValueError, PhyloZooTypeError


class TestFromCodedArray:
    """Tests for MSA.from_coded_array class method."""
    
    def test_from_coded_array_basic(self) -> None:
        """Test basic from_coded_array functionality."""
        coded = np.array([[0, 1, 2, 3], [0, 1, 2, 3]], dtype=np.int8)
        taxa_order = ("taxon1", "taxon2")
        msa = MSA.from_coded_array(coded, taxa_order)
        
        assert msa.sequence_length == 4
        assert msa.num_taxa == 2
        assert msa.taxa == {"taxon1", "taxon2"}
        assert msa.taxa_order == taxa_order
    
    def test_from_coded_array_roundtrip(self) -> None:
        """Test that from_coded_array produces same result as normal init."""
        sequences = {
            "taxon1": "ACGT",
            "taxon2": "ACGT"
        }
        msa1 = MSA(sequences)
        
        # Create from coded array
        msa2 = MSA.from_coded_array(msa1.coded_array, msa1.taxa_order)
        
        assert msa2.sequence_length == msa1.sequence_length
        assert msa2.num_taxa == msa1.num_taxa
        assert msa2.taxa == msa1.taxa
        for taxon in msa1.taxa:
            assert msa2.get_sequence(taxon) == msa1.get_sequence(taxon)
    
    def test_from_coded_array_invalid_type(self) -> None:
        """Test that non-numpy array raises TypeError."""
        taxa_order = ("taxon1", "taxon2")
        
        with pytest.raises(PhyloZooTypeError, match="must be a numpy array"):
            MSA.from_coded_array([[0, 1], [0, 1]], taxa_order)
    
    def test_from_coded_array_invalid_taxa_order_type(self) -> None:
        """Test that non-tuple taxa_order raises TypeError."""
        coded = np.array([[0, 1], [0, 1]], dtype=np.int8)
        
        with pytest.raises(PhyloZooTypeError, match="must be a tuple"):
            MSA.from_coded_array(coded, ["taxon1", "taxon2"])
    
    def test_from_coded_array_empty_array(self) -> None:
        """Test that empty array raises ValueError."""
        coded = np.array([], dtype=np.int8).reshape(0, 0)
        taxa_order = ()
        
        with pytest.raises(PhyloZooValueError, match="cannot be empty"):
            MSA.from_coded_array(coded, taxa_order)
    
    def test_from_coded_array_wrong_dimensions(self) -> None:
        """Test that wrong number of dimensions raises ValueError."""
        coded = np.array([0, 1, 2, 3], dtype=np.int8)  # 1D array
        taxa_order = ("taxon1",)
        
        with pytest.raises(PhyloZooValueError, match="2-dimensional"):
            MSA.from_coded_array(coded, taxa_order)
    
    def test_from_coded_array_wrong_dtype(self) -> None:
        """Test that wrong dtype raises ValueError."""
        coded = np.array([[0, 1], [0, 1]], dtype=np.int32)
        taxa_order = ("taxon1", "taxon2")
        
        with pytest.raises(PhyloZooValueError, match="dtype int8"):
            MSA.from_coded_array(coded, taxa_order)
    
    def test_from_coded_array_length_mismatch(self) -> None:
        """Test that taxa_order length mismatch raises ValueError."""
        coded = np.array([[0, 1], [0, 1]], dtype=np.int8)
        taxa_order = ("taxon1",)  # Wrong length
        
        with pytest.raises(PhyloZooValueError, match="must match"):
            MSA.from_coded_array(coded, taxa_order)
    
    def test_from_coded_array_empty_taxa_order(self) -> None:
        """Test that empty taxa_order raises ValueError."""
        coded = np.array([[]], dtype=np.int8).reshape(1, 0)
        taxa_order = ()
        
        with pytest.raises(PhyloZooValueError, match="cannot be empty"):
            MSA.from_coded_array(coded, taxa_order)
    
    def test_from_coded_array_non_string_taxon(self) -> None:
        """Test that non-string taxon raises TypeError."""
        coded = np.array([[0, 1]], dtype=np.int8)
        taxa_order = (123,)  # Not a string
        
        with pytest.raises(PhyloZooTypeError, match="must be a string"):
            MSA.from_coded_array(coded, taxa_order)
    
    def test_from_coded_array_immutability(self) -> None:
        """Test that modifying original array doesn't affect MSA."""
        coded = np.array([[0, 1, 2, 3]], dtype=np.int8)
        taxa_order = ("taxon1",)
        msa = MSA.from_coded_array(coded, taxa_order)
        
        # Modify original array
        coded[0, 0] = 99
        
        # MSA should be unaffected
        assert msa.coded_array[0, 0] == 0  # Still original value
    
    def test_from_coded_array_uses_default_encoding(self) -> None:
        """Test that DEFAULT_NUCLEOTIDE_ENCODING is used."""
        coded = np.array([[0, 1, 2, 3]], dtype=np.int8)
        taxa_order = ("taxon1",)
        msa = MSA.from_coded_array(coded, taxa_order)
        
        # Check that encoding matches default
        assert msa._character_encoding == DEFAULT_NUCLEOTIDE_ENCODING
        assert msa.get_sequence("taxon1") == "ACGT"
