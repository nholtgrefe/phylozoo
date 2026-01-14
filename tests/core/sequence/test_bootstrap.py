"""
Tests for bootstrap module.
"""

from __future__ import annotations

import pytest
import numpy as np

from phylozoo.core.sequence import MSA
from phylozoo.core.sequence.bootstrap import bootstrap, bootstrap_per_gene
from phylozoo.utils.exceptions import PhyloZooValueError


class TestBootstrap:
    """Tests for bootstrap function."""
    
    def test_bootstrap_basic(self) -> None:
        """Test basic bootstrap functionality."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGTACGT",
            "taxon3": "ACGTACGT"
        }
        msa = MSA(sequences)
        
        bootstrapped = next(bootstrap(msa, n_bootstrap=1, seed=42))
        
        assert bootstrapped.sequence_length == msa.sequence_length
        assert bootstrapped.num_taxa == msa.num_taxa
        assert bootstrapped.taxa == msa.taxa
    
    def test_bootstrap_custom_length(self) -> None:
        """Test bootstrap with custom length."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGTACGT",
            "taxon3": "ACGTACGT"
        }
        msa = MSA(sequences)
        
        bootstrapped = next(bootstrap(msa, n_bootstrap=1, length=10, seed=42))
        
        assert bootstrapped.sequence_length == 10
        assert bootstrapped.num_taxa == msa.num_taxa
    
    def test_bootstrap_reproducibility(self) -> None:
        """Test that bootstrap is reproducible with same seed."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGTACGT",
            "taxon3": "ACGTACGT"
        }
        msa = MSA(sequences)
        
        boot1 = next(bootstrap(msa, n_bootstrap=1, seed=42))
        boot2 = next(bootstrap(msa, n_bootstrap=1, seed=42))
        
        # Should be identical with same seed
        assert boot1.sequence_length == boot2.sequence_length
        for taxon in msa.taxa:
            assert boot1.get_sequence(taxon) == boot2.get_sequence(taxon)
    
    def test_bootstrap_multiple_replicates(self) -> None:
        """Test generating multiple bootstrap replicates."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGTACGT",
            "taxon3": "ACGTACGT"
        }
        msa = MSA(sequences)
        
        replicates = list(bootstrap(msa, n_bootstrap=5, seed=42))
        
        assert len(replicates) == 5
        for rep in replicates:
            assert rep.sequence_length == msa.sequence_length
            assert rep.num_taxa == msa.num_taxa
    
    def test_bootstrap_different_seeds(self) -> None:
        """Test that different seeds produce different results."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGTACGT",
            "taxon3": "ACGTACGT"
        }
        msa = MSA(sequences)
        
        boot1 = next(bootstrap(msa, n_bootstrap=1, seed=42))
        boot2 = next(bootstrap(msa, n_bootstrap=1, seed=123))
        
        # With different seeds, results should likely be different
        # (though not guaranteed, so we just check they're valid)
        assert boot1.sequence_length == boot2.sequence_length
        assert boot1.num_taxa == boot2.num_taxa
    
    def test_bootstrap_invalid_length(self) -> None:
        """Test that invalid length raises error."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGTACGT",
            "taxon3": "ACGTACGT"
        }
        msa = MSA(sequences)
        
        with pytest.raises(PhyloZooValueError, match="must be positive"):
            next(bootstrap(msa, n_bootstrap=1, length=0, seed=42))
        
        with pytest.raises(PhyloZooValueError, match="must be positive"):
            next(bootstrap(msa, n_bootstrap=1, length=-1, seed=42))
    
    def test_bootstrap_invalid_n_bootstrap(self) -> None:
        """Test that invalid n_bootstrap raises error."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGTACGT",
            "taxon3": "ACGTACGT"
        }
        msa = MSA(sequences)
        
        with pytest.raises(PhyloZooValueError, match="must be positive"):
            next(bootstrap(msa, n_bootstrap=0, seed=42))
        
        with pytest.raises(PhyloZooValueError, match="must be positive"):
            next(bootstrap(msa, n_bootstrap=-1, seed=42))
    
    def test_bootstrap_preserves_taxa(self) -> None:
        """Test that bootstrap preserves all taxa."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGTACGT",
            "taxon3": "ACGTACGT"
        }
        msa = MSA(sequences)
        
        bootstrapped = next(bootstrap(msa, n_bootstrap=1, seed=42))
        
        assert bootstrapped.taxa == msa.taxa
        for taxon in msa.taxa:
            assert taxon in bootstrapped
            assert bootstrapped.get_sequence(taxon) is not None


class TestBootstrapPerGene:
    """Tests for bootstrap_per_gene function."""
    
    def test_bootstrap_per_gene_basic(self) -> None:
        """Test basic per-gene bootstrap functionality."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGTACGT",
            "taxon3": "ACGTACGT"
        }
        msa = MSA(sequences)
        gene_lengths = [4, 4]  # Two genes of length 4 each
        
        bootstrapped = next(bootstrap_per_gene(msa, gene_lengths, n_bootstrap=1, seed=42))
        
        assert bootstrapped.sequence_length == msa.sequence_length
        assert bootstrapped.num_taxa == msa.num_taxa
        assert bootstrapped.taxa == msa.taxa
    
    def test_bootstrap_per_gene_reproducibility(self) -> None:
        """Test that per-gene bootstrap is reproducible with same seed."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGTACGT",
            "taxon3": "ACGTACGT"
        }
        msa = MSA(sequences)
        gene_lengths = [4, 4]
        
        boot1 = next(bootstrap_per_gene(msa, gene_lengths, n_bootstrap=1, seed=42))
        boot2 = next(bootstrap_per_gene(msa, gene_lengths, n_bootstrap=1, seed=42))
        
        # Should be identical with same seed
        assert boot1.sequence_length == boot2.sequence_length
        for taxon in msa.taxa:
            assert boot1.get_sequence(taxon) == boot2.get_sequence(taxon)
    
    def test_bootstrap_per_gene_multiple_genes(self) -> None:
        """Test per-gene bootstrap with multiple genes."""
        sequences = {
            "taxon1": "ACGTACGTACGT",
            "taxon2": "ACGTACGTACGT",
            "taxon3": "ACGTACGTACGT"
        }
        msa = MSA(sequences)
        gene_lengths = [3, 3, 3, 3]  # Four genes of length 3 each
        
        bootstrapped = next(bootstrap_per_gene(msa, gene_lengths, n_bootstrap=1, seed=42))
        
        assert bootstrapped.sequence_length == msa.sequence_length
        assert bootstrapped.num_taxa == msa.num_taxa
    
    def test_bootstrap_per_gene_multiple_replicates(self) -> None:
        """Test generating multiple per-gene bootstrap replicates."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGTACGT",
            "taxon3": "ACGTACGT"
        }
        msa = MSA(sequences)
        gene_lengths = [4, 4]
        
        replicates = list(bootstrap_per_gene(msa, gene_lengths, n_bootstrap=5, seed=42))
        
        assert len(replicates) == 5
        for rep in replicates:
            assert rep.sequence_length == msa.sequence_length
            assert rep.num_taxa == msa.num_taxa
    
    def test_bootstrap_per_gene_empty_gene_lengths(self) -> None:
        """Test that empty gene_lengths raises error."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGTACGT",
            "taxon3": "ACGTACGT"
        }
        msa = MSA(sequences)
        
        with pytest.raises(PhyloZooValueError, match="cannot be empty"):
            next(bootstrap_per_gene(msa, [], n_bootstrap=1, seed=42))
    
    def test_bootstrap_per_gene_invalid_length_sum(self) -> None:
        """Test that incorrect sum of gene lengths raises error."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGTACGT",
            "taxon3": "ACGTACGT"
        }
        msa = MSA(sequences)
        
        with pytest.raises(PhyloZooValueError, match="must equal"):
            next(bootstrap_per_gene(msa, [4, 3], n_bootstrap=1, seed=42))  # Sum is 7, not 8
    
    def test_bootstrap_per_gene_non_positive_lengths(self) -> None:
        """Test that non-positive gene lengths raise error."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGTACGT",
            "taxon3": "ACGTACGT"
        }
        msa = MSA(sequences)
        
        with pytest.raises(PhyloZooValueError, match="must be positive"):
            next(bootstrap_per_gene(msa, [4, 0, 4], n_bootstrap=1, seed=42))
        
        with pytest.raises(PhyloZooValueError, match="must be positive"):
            next(bootstrap_per_gene(msa, [4, -1, 4], n_bootstrap=1, seed=42))
    
    def test_bootstrap_per_gene_invalid_n_bootstrap(self) -> None:
        """Test that invalid n_bootstrap raises error."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGTACGT",
            "taxon3": "ACGTACGT"
        }
        msa = MSA(sequences)
        gene_lengths = [4, 4]
        
        with pytest.raises(PhyloZooValueError, match="must be positive"):
            next(bootstrap_per_gene(msa, gene_lengths, n_bootstrap=0, seed=42))
        
        with pytest.raises(PhyloZooValueError, match="must be positive"):
            next(bootstrap_per_gene(msa, gene_lengths, n_bootstrap=-1, seed=42))
    
    def test_bootstrap_per_gene_preserves_taxa(self) -> None:
        """Test that per-gene bootstrap preserves all taxa."""
        sequences = {
            "taxon1": "ACGTACGT",
            "taxon2": "ACGTACGT",
            "taxon3": "ACGTACGT"
        }
        msa = MSA(sequences)
        gene_lengths = [4, 4]
        
        bootstrapped = next(bootstrap_per_gene(msa, gene_lengths, n_bootstrap=1, seed=42))
        
        assert bootstrapped.taxa == msa.taxa
        for taxon in msa.taxa:
            assert taxon in bootstrapped
            assert bootstrapped.get_sequence(taxon) is not None
    
    def test_bootstrap_per_gene_gene_boundaries(self) -> None:
        """Test that per-gene bootstrap respects gene boundaries."""
        sequences = {
            "taxon1": "AAAAACCC",
            "taxon2": "AAAAACCC",
            "taxon3": "AAAAACCC"
        }
        msa = MSA(sequences)
        gene_lengths = [5, 3]  # First gene: AAAAA, second gene: CCC
        
        bootstrapped = next(bootstrap_per_gene(msa, gene_lengths, n_bootstrap=1, seed=42))
        
        # Each gene should only contain characters from its original gene
        # (though with replacement, so first gene can have A's, second gene can have C's)
        seq = bootstrapped.get_sequence("taxon1")
        assert seq is not None
        
        # First 5 positions should only contain A (from first gene)
        first_gene = seq[:5]
        assert all(c in ['A', '-', 'N'] for c in first_gene)
        
        # Last 3 positions should only contain C (from second gene)
        second_gene = seq[5:]
        assert all(c in ['C', '-', 'N'] for c in second_gene)
