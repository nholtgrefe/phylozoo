"""
Tests for _validate_branchlength_constraints method in DirectedPhyNetwork.
"""

import pytest

from phylozoo.core.network.dnetwork.base import DirectedPhyNetwork
from phylozoo.utils.validation import no_validation


class TestValidateBranchLengthConstraints:
    """Test cases for _validate_branchlength_constraints."""

    def test_no_parallel_edges(self) -> None:
        """Test that networks without parallel edges pass validation."""
        net = DirectedPhyNetwork(
            edges=[(1, 2), (1, 3)],
            nodes=[(2, {'label': 'A'}), (3, {'label': 'B'})]
        )
        # Should not raise
        net._validate_branchlength_constraints()

    def test_parallel_edges_same_branchlength(self) -> None:
        """Test that parallel edges with same branch_length pass validation."""
        # Use internal node as target to allow parallel edges
        net = DirectedPhyNetwork(
            edges=[
                {'u': 1, 'v': 3, 'branch_length': 0.5},
                {'u': 1, 'v': 3, 'branch_length': 0.5},
                (3, 2)  # Internal node to leaf
            ],
            nodes=[(2, {'label': 'A'})]
        )
        # Should not raise
        net._validate_branchlength_constraints()

    def test_parallel_edges_no_branchlength(self) -> None:
        """Test that parallel edges without branch_length pass validation."""
        # Use internal node as target to allow parallel edges
        net = DirectedPhyNetwork(
            edges=[
                {'u': 1, 'v': 3},
                {'u': 1, 'v': 3},
                (3, 2)  # Internal node to leaf
            ],
            nodes=[(2, {'label': 'A'})]
        )
        # Should not raise
        net._validate_branchlength_constraints()

    def test_parallel_edges_mixed_branchlength_error(self) -> None:
        """Test that parallel edges with some having branch_length and others not raises error."""
        # Use internal node as target to allow parallel edges
        # Create network with validation suppressed to avoid branch_length error during init
        with no_validation():
            net = DirectedPhyNetwork(
                edges=[
                    {'u': 1, 'v': 3, 'branch_length': 0.5},
                    {'u': 1, 'v': 3},  # Missing branch_length
                    (3, 2)  # Internal node to leaf
                ],
                nodes=[(2, {'label': 'A'})]
            )
        with pytest.raises(ValueError, match="inconsistent branch_length"):
            net._validate_branchlength_constraints()

    def test_parallel_edges_different_branchlength_error(self) -> None:
        """Test that parallel edges with different branch_length values raise error."""
        # Use internal node as target to allow parallel edges
        # Create network with validation suppressed to avoid branch_length error during init
        with no_validation():
            net = DirectedPhyNetwork(
                edges=[
                    {'u': 1, 'v': 3, 'branch_length': 0.5},
                    {'u': 1, 'v': 3, 'branch_length': 0.7},  # Different value
                    (3, 2)  # Internal node to leaf
                ],
                nodes=[(2, {'label': 'A'})]
            )
        with pytest.raises(ValueError, match="different branch_length values"):
            net._validate_branchlength_constraints()

    def test_parallel_edges_multiple_pairs(self) -> None:
        """Test validation with multiple pairs of parallel edges."""
        # Use internal nodes as targets to allow parallel edges
        net = DirectedPhyNetwork(
            edges=[
                {'u': 1, 'v': 4, 'branch_length': 0.5},
                {'u': 1, 'v': 4, 'branch_length': 0.5},
                {'u': 1, 'v': 5, 'branch_length': 0.3},
                {'u': 1, 'v': 5, 'branch_length': 0.3},
                (4, 2), (5, 3)  # Internal nodes to leaves
            ],
            nodes=[(2, {'label': 'A'}), (3, {'label': 'B'})]
        )
        # Should not raise
        net._validate_branchlength_constraints()

    def test_parallel_edges_one_pair_invalid(self) -> None:
        """Test that one invalid pair causes error even if others are valid."""
        # Use internal nodes as targets to allow parallel edges
        # Create network with validation suppressed to avoid branch_length error during init
        with no_validation():
            net = DirectedPhyNetwork(
                edges=[
                    {'u': 1, 'v': 4, 'branch_length': 0.5},
                    {'u': 1, 'v': 4, 'branch_length': 0.5},
                    {'u': 1, 'v': 5, 'branch_length': 0.3},
                    {'u': 1, 'v': 5},  # Missing branch_length
                    (4, 2), (5, 3)  # Internal nodes to leaves
                ],
                nodes=[(2, {'label': 'A'}), (3, {'label': 'B'})]
            )
        with pytest.raises(ValueError, match="inconsistent branch_length"):
            net._validate_branchlength_constraints()

    def test_three_parallel_edges_all_same(self) -> None:
        """Test that three parallel edges with same branch_length pass validation."""
        # Use internal node as target to allow parallel edges
        net = DirectedPhyNetwork(
            edges=[
                {'u': 1, 'v': 3, 'branch_length': 0.5},
                {'u': 1, 'v': 3, 'branch_length': 0.5},
                {'u': 1, 'v': 3, 'branch_length': 0.5},
                (3, 2)  # Internal node to leaf
            ],
            nodes=[(2, {'label': 'A'})]
        )
        # Should not raise
        net._validate_branchlength_constraints()

    def test_three_parallel_edges_one_different(self) -> None:
        """Test that three parallel edges with one different branch_length raise error."""
        # Use internal node as target to allow parallel edges
        # Create network with validation suppressed to avoid branch_length error during init
        with no_validation():
            net = DirectedPhyNetwork(
                edges=[
                    {'u': 1, 'v': 3, 'branch_length': 0.5},
                    {'u': 1, 'v': 3, 'branch_length': 0.5},
                    {'u': 1, 'v': 3, 'branch_length': 0.7},  # Different
                    (3, 2)  # Internal node to leaf
                ],
                nodes=[(2, {'label': 'A'})]
            )
        with pytest.raises(ValueError, match="different branch_length values"):
            net._validate_branchlength_constraints()

    def test_parallel_edges_numeric_tolerance(self) -> None:
        """Test that very small differences in branch_length are detected."""
        # Use internal node as target to allow parallel edges
        # 1e-11 is smaller than tolerance (1e-10), so should pass
        net = DirectedPhyNetwork(
            edges=[
                {'u': 1, 'v': 3, 'branch_length': 0.5},
                {'u': 1, 'v': 3, 'branch_length': 0.5 + 1e-11},  # Very small difference, within tolerance
                (3, 2)  # Internal node to leaf
            ],
            nodes=[(2, {'label': 'A'})]
        )
        # Should pass (within tolerance)
        net._validate_branchlength_constraints()

    def test_parallel_edges_numeric_tolerance_exceeded(self) -> None:
        """Test that differences larger than tolerance are detected."""
        # Use internal node as target to allow parallel edges
        # Create network with validation suppressed to avoid branch_length error during init
        with no_validation():
            net = DirectedPhyNetwork(
                edges=[
                    {'u': 1, 'v': 3, 'branch_length': 0.5},
                    {'u': 1, 'v': 3, 'branch_length': 0.5 + 1e-8},  # Larger than tolerance
                    (3, 2)  # Internal node to leaf
                ],
                nodes=[(2, {'label': 'A'})]
            )
        # abs(0.5 + 1e-8 - 0.5) = 1e-8, which is > 1e-10, so should fail
        with pytest.raises(ValueError, match="different branch_length values"):
            net._validate_branchlength_constraints()

    def test_empty_network(self) -> None:
        """Test that empty network passes validation."""
        net = DirectedPhyNetwork(edges=[], nodes=[])
        # Should not raise
        net._validate_branchlength_constraints()

    def test_single_edge(self) -> None:
        """Test that single edge (no parallel) passes validation."""
        net = DirectedPhyNetwork(
            edges=[{'u': 1, 'v': 2, 'branch_length': 0.5}],
            nodes=[(2, {'label': 'A'})]
        )
        # Should not raise
        net._validate_branchlength_constraints()

