"""
Tests for base diversity functions.

This test suite covers:
- diversity() function
- marginal_diversities() function
- greedy_max_diversity() function
- solve_max_diversity() function
"""

import pytest

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.panda.measure import (
    diversity,
    greedy_max_diversity,
    marginal_diversities,
    solve_max_diversity,
    all_paths,
)


class TestDiversityFunction:
    """Test cases for the diversity() function."""

    def test_diversity_basic(self) -> None:
        """Test basic diversity computation."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'branch_length': 0.5},
                {'u': 3, 'v': 2, 'branch_length': 0.3}
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        div = diversity(net, {"A", "B"}, all_paths)
        assert div == 0.8

    def test_diversity_single_taxon(self) -> None:
        """Test diversity computation for a single taxon."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': 0.5}],
            nodes=[(1, {'label': 'A'})]
        )
        
        div = diversity(net, {"A"}, all_paths)
        assert div == 0.5

    def test_diversity_empty_set(self) -> None:
        """Test diversity computation for empty set."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': 0.5}],
            nodes=[(1, {'label': 'A'})]
        )
        
        div = diversity(net, set(), all_paths)
        assert div == 0.0

    def test_diversity_invalid_taxa_raises_error(self) -> None:
        """Test that diversity raises error for invalid taxa."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': 0.5}],
            nodes=[(1, {'label': 'A'})]
        )
        
        with pytest.raises(ValueError, match="Taxa not found in network"):
            diversity(net, {"A", "C"}, all_paths)

    def test_diversity_all_invalid_taxa_raises_error(self) -> None:
        """Test that diversity raises error when all taxa are invalid."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': 0.5}],
            nodes=[(1, {'label': 'A'})]
        )
        
        with pytest.raises(ValueError, match="Taxa not found in network"):
            diversity(net, {"X", "Y"}, all_paths)


class TestMarginalDiversities:
    """Test cases for the marginal_diversities() function."""

    def test_marginal_diversities_basic(self) -> None:
        """Test basic marginal diversity computation."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'branch_length': 0.5},
                {'u': 3, 'v': 2, 'branch_length': 0.3}
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        marginals = marginal_diversities(net, {"A"}, all_paths)
        
        # Should have entries for both A and B
        assert "A" in marginals
        assert "B" in marginals
        
        # Marginal of A (in saved set) should be negative (decrease when removed)
        assert marginals["A"] < 0
        
        # Marginal of B (not in saved set) should be positive (increase when added)
        assert marginals["B"] > 0

    def test_marginal_diversities_empty_saved(self) -> None:
        """Test marginal diversities with empty saved set."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'branch_length': 0.5},
                {'u': 3, 'v': 2, 'branch_length': 0.3}
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        marginals = marginal_diversities(net, set(), all_paths)
        
        # All marginals should be positive (adding to empty set)
        assert all(val > 0 for val in marginals.values())

    def test_marginal_diversities_invalid_taxa_raises_error(self) -> None:
        """Test that marginal_diversities raises error for invalid taxa."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': 0.5}],
            nodes=[(1, {'label': 'A'})]
        )
        
        with pytest.raises(ValueError, match="Taxa not found in network"):
            marginal_diversities(net, {"A", "C"}, all_paths)


class TestGreedyMaxDiversity:
    """Test cases for the greedy_max_diversity() function."""

    def test_greedy_max_diversity_basic(self) -> None:
        """Test basic greedy max diversity."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'branch_length': 0.5},
                {'u': 3, 'v': 2, 'branch_length': 0.3}
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        value, solution = greedy_max_diversity(net, k=2, measure=all_paths)
        
        assert isinstance(value, float)
        assert isinstance(solution, set)
        assert len(solution) == 2
        assert solution == {"A", "B"}
        assert value == 0.8

    def test_greedy_max_diversity_k_equals_one(self) -> None:
        """Test greedy max diversity with k=1."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'branch_length': 0.5},
                {'u': 3, 'v': 2, 'branch_length': 0.3}
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        value, solution = greedy_max_diversity(net, k=1, measure=all_paths)
        
        assert len(solution) == 1
        # Should select the taxon with higher branch length (A)
        assert solution == {"A"}
        assert value == 0.5

    def test_greedy_max_diversity_k_equals_zero(self) -> None:
        """Test greedy max diversity with k=0."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': 0.5}],
            nodes=[(1, {'label': 'A'})]
        )
        
        value, solution = greedy_max_diversity(net, k=0, measure=all_paths)
        
        assert len(solution) == 0
        assert value == 0.0

    def test_greedy_max_diversity_invalid_k_negative(self) -> None:
        """Test that greedy_max_diversity raises error for negative k."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': 0.5}],
            nodes=[(1, {'label': 'A'})]
        )
        
        with pytest.raises(ValueError, match="k must be non-negative"):
            greedy_max_diversity(net, k=-1, measure=all_paths)

    def test_greedy_max_diversity_invalid_k_too_large(self) -> None:
        """Test that greedy_max_diversity raises error for k > number of taxa."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': 0.5}],
            nodes=[(1, {'label': 'A'})]
        )
        
        with pytest.raises(ValueError, match="k must be <="):
            greedy_max_diversity(net, k=5, measure=all_paths)


class TestSolveMaxDiversity:
    """Test cases for the solve_max_diversity() function."""

    def test_solve_max_diversity_works(self) -> None:
        """Test that solve_max_diversity works for all_paths."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'branch_length': 0.5},
                {'u': 3, 'v': 2, 'branch_length': 0.3}
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        # all_paths now implements solve_maximization
        value, solution = solve_max_diversity(net, k=1, measure=all_paths)
        assert isinstance(value, float)
        assert isinstance(solution, set)
        assert len(solution) == 1
        assert solution.issubset({"A", "B"})

    def test_solve_max_diversity_invalid_k_negative(self) -> None:
        """Test that solve_max_diversity validates k through measure."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': 0.5}],
            nodes=[(1, {'label': 'A'})]
        )
        
        # solve_maximization validates k first and raises ValueError
        with pytest.raises(ValueError, match="k must be between"):
            solve_max_diversity(net, k=-1, measure=all_paths)

