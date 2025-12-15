"""
Tests for the MSA module.
"""

import pytest
import numpy as np

from phylozoo.core.sequence.msa import MSA, DEFAULT_NUCLEOTIDE_ENCODING


class TestMSACreation:
    """Test cases for MSA creation and basic properties."""
    
    def test_basic_creation(self) -> None:
        """Test creating a valid MSA."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGTACGT",
            "taxon3": "ACGTACGT",
        }
        msa = MSA(sequences)
        assert len(msa) == 3
        assert msa.sequence_length == 8
        assert msa.num_taxa == 3
        assert len(msa.taxa) == 3
    
    def test_single_taxon(self) -> None:
        """Test MSA with a single taxon."""
        sequences = {"taxon1": "ACGT"}
        msa = MSA(sequences)
        assert len(msa) == 1
        assert msa.sequence_length == 4
    
    def test_empty_sequences_raises_error(self) -> None:
        """Test that empty sequences dictionary raises ValueError."""
        with pytest.raises(ValueError, match="cannot be empty"):
            MSA({})
    
    def test_non_dict_raises_error(self) -> None:
        """Test that non-dictionary input raises TypeError."""
        with pytest.raises(TypeError, match="must be a dictionary"):
            MSA("not a dict")  # type: ignore
    
    def test_different_lengths_raises_error(self) -> None:
        """Test that sequences with different lengths raise ValueError."""
        sequences = {
            "taxon1": "ACGT",
            "taxon2": "ACGTACGT",  # Different length
        }
        with pytest.raises(ValueError, match="same length"):
            MSA(sequences)
    
    def test_non_string_sequence_raises_error(self) -> None:
        """Test that non-string sequences raise TypeError."""
        sequences = {
            "taxon1": "ACGT",
            "taxon2": 123,  # type: ignore
        }
        with pytest.raises(TypeError, match="must be a string"):
            MSA(sequences)
    
    def test_non_string_taxon_raises_error(self) -> None:
        """Test that non-string taxon names raise TypeError."""
        sequences = {
            "taxon1": "ACGT",
            123: "ACGT",  # type: ignore
        }
        with pytest.raises(TypeError, match="must be strings"):
            MSA(sequences)
    
    def test_canonical_ordering(self) -> None:
        """Test that taxa are stored in canonical (sorted) order."""
        sequences = {
            "zebra": "ACGT",
            "alpha": "ACGT",
            "beta": "ACGT",
        }
        msa = MSA(sequences)
        assert msa.taxa_order == ("alpha", "beta", "zebra")
        assert msa.taxa == frozenset({"alpha", "beta", "zebra"})




class TestMSAAccess:
    """Test cases for accessing MSA data."""
    
    def test_get_sequence(self) -> None:
        """Test getting a sequence by taxon name."""
        sequences = {"taxon1": "ACGT", "taxon2": "TGCA"}
        msa = MSA(sequences)
        
        assert msa.get_sequence("taxon1") == "ACGT"
        assert msa.get_sequence("taxon2") == "TGCA"
        assert msa.get_sequence("nonexistent") is None
    
    def test_get_sequences(self) -> None:
        """Test getting all sequences as a dictionary."""
        sequences = {"taxon1": "ACGT", "taxon2": "TGCA"}
        msa = MSA(sequences)
        
        all_sequences = msa.get_sequences()
        assert all_sequences == {"taxon1": "ACGT", "taxon2": "TGCA"}
        
        # Verify it's a copy (modifying doesn't affect original)
        all_sequences["taxon1"] = "XXXX"
        assert msa.get_sequence("taxon1") == "ACGT"
    
    def test_contains(self) -> None:
        """Test checking if a taxon is in the MSA."""
        sequences = {"taxon1": "ACGT", "taxon2": "TGCA"}
        msa = MSA(sequences)
        
        assert "taxon1" in msa
        assert "taxon2" in msa
        assert "taxon3" not in msa
    
    def test_len(self) -> None:
        """Test getting the number of taxa."""
        sequences = {"taxon1": "ACGT", "taxon2": "TGCA", "taxon3": "AAAA"}
        msa = MSA(sequences)
        assert len(msa) == 3


class TestMSAProperties:
    """Test cases for MSA properties."""
    
    def test_taxa_property(self) -> None:
        """Test taxa property returns frozenset."""
        sequences = {"taxon1": "ACGT", "taxon2": "TGCA"}
        msa = MSA(sequences)
        
        taxa = msa.taxa
        assert isinstance(taxa, frozenset)
        assert taxa == frozenset({"taxon1", "taxon2"})
    
    def test_taxa_order_property(self) -> None:
        """Test taxa_order property returns tuple in sorted order."""
        sequences = {"zebra": "ACGT", "alpha": "ACGT", "beta": "ACGT"}
        msa = MSA(sequences)
        
        taxa_order = msa.taxa_order
        assert isinstance(taxa_order, tuple)
        assert taxa_order == ("alpha", "beta", "zebra")
    
    def test_sequence_length_property(self) -> None:
        """Test sequence_length property."""
        sequences = {"taxon1": "ACGTACGT", "taxon2": "ACGTACGT"}
        msa = MSA(sequences)
        assert msa.sequence_length == 8


class TestMSACharacterEncoding:
    """Test cases for character encoding."""
    
    def test_default_encoding(self) -> None:
        """Test that default nucleotide encoding is used."""
        sequences = {"taxon1": "ACGT", "taxon2": "ACGT"}
        msa = MSA(sequences)
        
        coded = msa._coded_array()
        assert coded.shape == (2, 4)
        assert coded.dtype == np.int8
        
        # Check encoding: A=0, C=1, G=2, T=3
        assert coded[0, 0] == 0  # A
        assert coded[0, 1] == 1  # C
        assert coded[0, 2] == 2  # G
        assert coded[0, 3] == 3  # T
    
    def test_custom_encoding(self) -> None:
        """Test custom character encoding."""
        custom_encoding = {"A": 0, "B": 1, "C": 2, "-": -1}
        sequences = {"taxon1": "ABC", "taxon2": "ABC"}
        msa = MSA(sequences, character_encoding=custom_encoding)
        
        coded = msa._coded_array()
        assert coded.shape == (2, 3)
        assert coded[0, 0] == 0  # A
        assert coded[0, 1] == 1  # B
        assert coded[0, 2] == 2  # C
    
    def test_unknown_characters_encoded_as_minus_one(self) -> None:
        """Test that unknown characters are encoded as -1."""
        sequences = {"taxon1": "AXGT", "taxon2": "AXGT"}  # X is not in encoding
        msa = MSA(sequences)
        
        coded = msa._coded_array()
        assert coded[0, 1] == -1  # X should be -1
    
    def test_gaps_encoded_as_minus_one(self) -> None:
        """Test that gaps are encoded as -1."""
        sequences = {"taxon1": "AC-GT", "taxon2": "AC-GT"}
        msa = MSA(sequences)
        
        coded = msa._coded_array()
        assert coded[0, 2] == -1  # Gap should be -1


class TestMSARepresentation:
    """Test cases for MSA string representations."""
    
    def test_repr(self) -> None:
        """Test __repr__ method."""
        sequences = {"taxon1": "ACGT", "taxon2": "TGCA"}
        msa = MSA(sequences)
        
        repr_str = repr(msa)
        assert "MSA" in repr_str
        assert "num_taxa=2" in repr_str
        assert "sequence_length=4" in repr_str
    
    def test_str(self) -> None:
        """Test __str__ method."""
        sequences = {"taxon1": "ACGT", "taxon2": "TGCA"}
        msa = MSA(sequences)
        
        str_repr = str(msa)
        assert "MSA" in str_repr
        assert "2 taxa" in str_repr
        assert "sequence length 4" in str_repr
        assert "taxon1" in str_repr
        assert "taxon2" in str_repr




class TestMSACopy:
    """Test cases for MSA copying."""
    
    def test_copy(self) -> None:
        """Test creating a copy of an MSA."""
        sequences = {"taxon1": "ACGT", "taxon2": "TGCA"}
        msa1 = MSA(sequences)
        msa2 = msa1.copy()
        
        # Verify they are different objects
        assert msa1 is not msa2
        
        # Verify they have the same data
        assert msa1.get_sequence("taxon1") == msa2.get_sequence("taxon1")
        assert msa1.get_sequence("taxon2") == msa2.get_sequence("taxon2")
        assert msa1.taxa == msa2.taxa
        assert msa1.sequence_length == msa2.sequence_length
        assert msa1.num_taxa == msa2.num_taxa


class TestMSACodedArray:
    """Test cases for coded array computation."""
    
    def test_coded_array_shape(self) -> None:
        """Test that coded array has correct shape."""
        sequences = {"taxon1": "ACGT", "taxon2": "TGCA", "taxon3": "AAAA"}
        msa = MSA(sequences)
        
        coded = msa._coded_array()
        assert coded.shape == (3, 4)  # (num_taxa, sequence_length)
        assert coded.dtype == np.int8
    
    def test_coded_array_values(self) -> None:
        """Test that coded array has correct values."""
        sequences = {"taxon1": "ACGT"}
        msa = MSA(sequences)
        
        coded = msa._coded_array()
        # A=0, C=1, G=2, T=3
        assert coded[0, 0] == 0
        assert coded[0, 1] == 1
        assert coded[0, 2] == 2
        assert coded[0, 3] == 3
    
    def test_coded_array_case_insensitive(self) -> None:
        """Test that coded array handles case-insensitive encoding."""
        sequences = {"taxon1": "acgt", "taxon2": "ACGT"}
        msa = MSA(sequences)
        
        coded = msa._coded_array()
        # Both should be encoded the same
        assert np.array_equal(coded[0, :], coded[1, :])
    
    def test_coded_array_gaps(self) -> None:
        """Test that gaps are encoded correctly."""
        sequences = {"taxon1": "AC-GT", "taxon2": "AC-GT"}
        msa = MSA(sequences)
        
        coded = msa._coded_array()
        assert coded[0, 2] == -1  # Gap
        assert coded[0, 0] == 0   # A
        assert coded[0, 1] == 1   # C
        assert coded[0, 3] == 2   # G
        assert coded[0, 4] == 3   # T

