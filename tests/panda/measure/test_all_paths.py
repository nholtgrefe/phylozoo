"""
Tests for the AllPathsDiversity measure.

This test suite covers:
- compute_diversity() method
- solve_maximization() method
- Edge cases and error handling
"""

import pytest

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.panda.measure.all_paths import AllPathsDiversity, all_paths
from tests.fixtures.directed_networks import (
    DTREE_SMALL_BINARY,
    LEVEL_1_DNETWORK_SINGLE_HYBRID_BINARY,
    LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE,
    LEVEL_1_DNETWORK_PARALLEL_EDGES,
)


class TestAllPathsDiversity:
    """Test cases for AllPathsDiversity class."""

    def test_compute_diversity_basic(self) -> None:
        """Test basic diversity computation."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'branch_length': 0.5},
                {'u': 3, 'v': 2, 'branch_length': 0.3}
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        measure = AllPathsDiversity()
        div = measure.compute_diversity(net, {"A", "B"})
        assert div == 0.8

    def test_compute_diversity_single_taxon(self) -> None:
        """Test diversity computation for a single taxon."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': 0.5}],
            nodes=[(1, {'label': 'A'})]
        )
        
        measure = AllPathsDiversity()
        div = measure.compute_diversity(net, {"A"})
        assert div == 0.5

    def test_compute_diversity_empty_set(self) -> None:
        """Test diversity computation for empty set."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': 0.5}],
            nodes=[(1, {'label': 'A'})]
        )
        
        measure = AllPathsDiversity()
        div = measure.compute_diversity(net, set())
        assert div == 0.0

    def test_compute_diversity_no_branch_lengths(self) -> None:
        """Test diversity computation when branch lengths are not specified."""
        net = DirectedPhyNetwork(
            edges=[
                (3, 1),  # No branch_length
                (3, 2)   # No branch_length
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        measure = AllPathsDiversity()
        div = measure.compute_diversity(net, {"A", "B"})
        # Should default to 1.0 per edge
        assert div == 2.0

    def test_compute_diversity_shared_ancestors(self) -> None:
        """Test that shared ancestors are counted once."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 4, 'v': 1, 'branch_length': 0.5},
                {'u': 4, 'v': 2, 'branch_length': 0.3},
                {'u': 3, 'v': 4, 'branch_length': 0.2}
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        measure = AllPathsDiversity()
        # Diversity of {A, B} should include all edges
        div = measure.compute_diversity(net, {"A", "B"})
        assert div == 0.5 + 0.3 + 0.2

    def test_solve_maximization_small_tree(self) -> None:
        """Test solve_maximization on a small binary tree."""
        measure = AllPathsDiversity()
        
        # Test with k=1: should select one taxon
        value, solution = measure.solve_maximization(DTREE_SMALL_BINARY, k=1)
        assert len(solution) == 1
        assert value > 0.0
        # Solution should be one of the taxa
        assert solution.issubset({"A", "B", "C"})
        
        # Test with k=2: should select two taxa
        value2, solution2 = measure.solve_maximization(DTREE_SMALL_BINARY, k=2)
        assert len(solution2) == 2
        assert value2 > value  # More taxa should give higher diversity
        assert solution2.issubset({"A", "B", "C"})
        
        # Test with k=3: should select all taxa
        value3, solution3 = measure.solve_maximization(DTREE_SMALL_BINARY, k=3)
        assert len(solution3) == 3
        assert value3 >= value2
        assert solution3 == {"A", "B", "C"}
    
    def test_solve_maximization_level1_single_hybrid(self) -> None:
        """Test solve_maximization on a level-1 network with single hybrid."""
        measure = AllPathsDiversity()
        
        # Test with k=1
        value, solution = measure.solve_maximization(
            LEVEL_1_DNETWORK_SINGLE_HYBRID_BINARY, k=1
        )
        assert len(solution) == 1
        assert value > 0.0
        assert solution.issubset({"A", "B", "C", "D"})
        
        # Test with k=2
        value2, solution2 = measure.solve_maximization(
            LEVEL_1_DNETWORK_SINGLE_HYBRID_BINARY, k=2
        )
        assert len(solution2) == 2
        assert value2 > value
        assert solution2.issubset({"A", "B", "C", "D"})
    
    def test_solve_maximization_level1_two_hybrids(self) -> None:
        """Test solve_maximization on a level-1 network with two separate hybrids."""
        measure = AllPathsDiversity()
        
        # Test with k=1
        value, solution = measure.solve_maximization(
            LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE, k=1
        )
        assert len(solution) == 1
        assert value > 0.0
        assert solution.issubset({"A", "B", "C", "D", "E"})
        
        # Test with k=3
        value3, solution3 = measure.solve_maximization(
            LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE, k=3
        )
        assert len(solution3) == 3
        assert value3 > value
        assert solution3.issubset({"A", "B", "C", "D", "E"})
    
    def test_solve_maximization_k_zero(self) -> None:
        """Test solve_maximization with k=0."""
        measure = AllPathsDiversity()
        
        value, solution = measure.solve_maximization(DTREE_SMALL_BINARY, k=0)
        assert len(solution) == 0
        assert value == 0.0
    
    def test_solve_maximization_invalid_k(self) -> None:
        """Test solve_maximization with invalid k values."""
        measure = AllPathsDiversity()
        
        # k negative
        with pytest.raises(ValueError, match="k must be between"):
            measure.solve_maximization(DTREE_SMALL_BINARY, k=-1)
        
        # k too large
        with pytest.raises(ValueError, match="k must be between"):
            measure.solve_maximization(DTREE_SMALL_BINARY, k=10)
    
    def test_solve_maximization_parallel_edges_error(self) -> None:
        """Test that solve_maximization raises error for networks with parallel edges."""
        measure = AllPathsDiversity()
        
        with pytest.raises(ValueError, match="parallel edges"):
            measure.solve_maximization(LEVEL_1_DNETWORK_PARALLEL_EDGES, k=1)
    
    def test_solve_maximization_consistency(self) -> None:
        """Test that solve_maximization returns consistent results."""
        measure = AllPathsDiversity()
        
        # Run multiple times and check consistency
        value1, solution1 = measure.solve_maximization(DTREE_SMALL_BINARY, k=2)
        value2, solution2 = measure.solve_maximization(DTREE_SMALL_BINARY, k=2)
        
        assert value1 == value2
        assert solution1 == solution2
        
        # Check that the diversity value matches compute_diversity
        computed_div = measure.compute_diversity(DTREE_SMALL_BINARY, solution1)
        assert abs(value1 - computed_div) < 1e-10


class TestAllPathsInstance:
    """Test cases for the all_paths convenience instance."""

    def test_all_paths_is_instance(self) -> None:
        """Test that all_paths is an instance of AllPathsDiversity."""
        assert isinstance(all_paths, AllPathsDiversity)

    def test_all_paths_compute_diversity(self) -> None:
        """Test that all_paths instance works correctly."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': 0.5}],
            nodes=[(1, {'label': 'A'})]
        )
        
        div = all_paths.compute_diversity(net, {"A"})
        assert div == 0.5

