"""
Tests for _validate_branchlength_constraints method in MixedPhyNetwork and SemiDirectedPhyNetwork.
"""

import pytest

from phylozoo.core.network.sdnetwork.base import MixedPhyNetwork
from phylozoo.core.network.sdnetwork.sd_phynetwork import SemiDirectedPhyNetwork
from phylozoo.utils.validation import no_validation


class TestValidateBranchLengthConstraintsMixed:
    """Test cases for _validate_branchlength_constraints in MixedPhyNetwork."""

    def test_no_parallel_edges(self) -> None:
        """Test that networks without parallel edges pass validation."""
        # Internal node must have degree >= 3
        net = MixedPhyNetwork(
            undirected_edges=[(1, 2), (1, 3), (1, 4)],
            nodes=[(2, {'label': 'A'}), (3, {'label': 'B'}), (4, {'label': 'C'})]
        )
        # Should not raise
        net._validate_branchlength_constraints()

    def test_parallel_undirected_edges_same_branchlength(self) -> None:
        """Test that parallel undirected edges with same branch_length pass validation."""
        # Internal nodes must have degree >= 3 (parallel edges don't increase degree)
        # Node 1: edges to 2 (parallel, counts as 1), 3, 4 = degree 3
        # Node 2: edges to 1 (parallel, counts as 1), 5, 6 = degree 3
        net = MixedPhyNetwork(
            undirected_edges=[
                {'u': 1, 'v': 2, 'branch_length': 0.5},
                {'u': 1, 'v': 2, 'branch_length': 0.5},
                (1, 3), (1, 4), (2, 5), (2, 6)  # Additional edges to satisfy degree constraints
            ],
            nodes=[(3, {'label': 'A'}), (4, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
        )
        # Should not raise
        net._validate_branchlength_constraints()

    def test_parallel_directed_edges_same_branchlength(self) -> None:
        """Test that parallel directed edges with same branch_length pass validation."""
        # Internal node must have degree >= 3, and hybrid node needs proper structure
        net = MixedPhyNetwork(
            directed_edges=[
                {'u': 1, 'v': 3, 'gamma': 0.5, 'branch_length': 0.5},
                {'u': 2, 'v': 3, 'gamma': 0.5, 'branch_length': 0.5}
            ],
            undirected_edges=[(3, 4)],
            nodes=[(4, {'label': 'A'})]
        )
        # Should not raise
        net._validate_branchlength_constraints()

    def test_parallel_undirected_edges_no_branchlength(self) -> None:
        """Test that parallel undirected edges without branch_length pass validation."""
        # Internal nodes must have degree >= 3 (parallel edges don't increase degree)
        # Node 1: edges to 2 (parallel, counts as 1), 3, 4 = degree 3
        # Node 2: edges to 1 (parallel, counts as 1), 5, 6 = degree 3
        net = MixedPhyNetwork(
            undirected_edges=[
                {'u': 1, 'v': 2},
                {'u': 1, 'v': 2},
                (1, 3), (1, 4), (2, 5), (2, 6)  # Additional edges to satisfy degree constraints
            ],
            nodes=[(3, {'label': 'A'}), (4, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
        )
        # Should not raise
        net._validate_branchlength_constraints()

    def test_parallel_undirected_edges_mixed_branchlength_error(self) -> None:
        """Test that parallel undirected edges with some having branch_length and others not raises error."""
        # Internal nodes must have degree >= 3 (parallel edges don't increase degree)
        # Node 1: edges to 2 (parallel, counts as 1), 3, 4 = degree 3
        # Node 2: edges to 1 (parallel, counts as 1), 5, 6 = degree 3
        # Create network with validation suppressed to avoid branch_length error during init
        with no_validation():
            net = MixedPhyNetwork(
                undirected_edges=[
                    {'u': 1, 'v': 2, 'branch_length': 0.5},
                    {'u': 1, 'v': 2},  # Missing branch_length
                    (1, 3), (1, 4), (2, 5), (2, 6)  # Additional edges to satisfy degree constraints
                ],
                nodes=[(3, {'label': 'A'}), (4, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
            )
        with pytest.raises(ValueError, match="inconsistent branch_length"):
            net._validate_branchlength_constraints()

    def test_parallel_directed_edges_mixed_branchlength_error(self) -> None:
        """Test that parallel directed edges with some having branch_length and others not raises error."""
        # Parallel directed edges from same source to same hybrid (with explicit keys)
        # Gamma: 0.5 + 0.5 = 1.0
        # Create network with validation suppressed to avoid branch_length error during init
        with no_validation():
            net = MixedPhyNetwork(
                directed_edges=[
                    {'u': 5, 'v': 4, 'key': 0, 'gamma': 0.5, 'branch_length': 0.5},
                    {'u': 5, 'v': 4, 'key': 1, 'gamma': 0.5},  # Missing branch_length
                    {'u': 6, 'v': 4, 'gamma': 0.0}  # Another source to hybrid (gamma 0.0 means no gamma)
                ],
                undirected_edges=[(4, 1), (5, 2), (5, 3), (6, 7), (6, 8)],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'}), (8, {'label': 'E'})]
            )
        with pytest.raises(ValueError, match="inconsistent branch_length"):
            net._validate_branchlength_constraints()

    def test_parallel_undirected_edges_different_branchlength_error(self) -> None:
        """Test that parallel undirected edges with different branch_length values raise error."""
        # Internal nodes must have degree >= 3
        # Create network with validation suppressed to avoid branch_length error during init
        with no_validation():
            net = MixedPhyNetwork(
                undirected_edges=[
                    {'u': 1, 'v': 2, 'branch_length': 0.5},
                    {'u': 1, 'v': 2, 'branch_length': 0.7},  # Different value
                    (1, 3), (1, 6), (2, 4), (2, 5), (2, 7)  # Additional edges to satisfy degree constraints
                ],
                nodes=[(3, {'label': 'A'}), (4, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'}), (7, {'label': 'E'})]
            )
        with pytest.raises(ValueError, match="different branch_length values"):
            net._validate_branchlength_constraints()

    def test_parallel_directed_edges_different_branchlength_error(self) -> None:
        """Test that parallel directed edges with different branch_length values raise error."""
        # Parallel directed edges from same source to same hybrid (with explicit keys)
        # Gamma: 0.33 + 0.33 + 0.34 = 1.0
        # Create network with validation suppressed to avoid branch_length error during init
        with no_validation():
            net = MixedPhyNetwork(
                directed_edges=[
                    {'u': 5, 'v': 4, 'key': 0, 'gamma': 0.33, 'branch_length': 0.5},
                    {'u': 5, 'v': 4, 'key': 1, 'gamma': 0.33, 'branch_length': 0.7},  # Different value
                    {'u': 6, 'v': 4, 'gamma': 0.34}  # Another source to hybrid
                ],
                undirected_edges=[(4, 1), (5, 2), (5, 3), (6, 7), (6, 8)],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'}), (8, {'label': 'E'})]
            )
        with pytest.raises(ValueError, match="different branch_length values"):
            net._validate_branchlength_constraints()

    def test_both_directed_and_undirected_parallel(self) -> None:
        """Test validation with both directed and undirected parallel edges."""
        # Need proper structure with hybrid node and connected network
        # Hybrid node 4: indegree 2 (from 5 twice), so total_degree must be 3
        # So it can have 1 undirected edge
        # For parallel undirected edges, use a different internal node
        # Connect everything properly
        net = MixedPhyNetwork(
            directed_edges=[
                {'u': 5, 'v': 4, 'key': 0, 'gamma': 0.5, 'branch_length': 0.5},
                {'u': 5, 'v': 4, 'key': 1, 'gamma': 0.5, 'branch_length': 0.5}
            ],
            undirected_edges=[
                {'u': 6, 'v': 7, 'branch_length': 0.3},
                {'u': 6, 'v': 7, 'branch_length': 0.3},
                (4, 1), (5, 2), (5, 3), (5, 6), (6, 8), (6, 9), (7, 10), (7, 11)  # Connect 5 and 6
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (8, {'label': 'D'}), (9, {'label': 'E'}), (10, {'label': 'F'}), (11, {'label': 'G'})]
        )
        # Should not raise
        net._validate_branchlength_constraints()

    def test_empty_network(self) -> None:
        """Test that empty network passes validation."""
        net = MixedPhyNetwork(undirected_edges=[], directed_edges=[], nodes=[])
        # Should not raise
        net._validate_branchlength_constraints()


class TestValidateBranchLengthConstraintsSemiDirected:
    """Test cases for _validate_branchlength_constraints in SemiDirectedPhyNetwork."""

    def test_parallel_undirected_edges_same_branchlength(self) -> None:
        """Test that parallel undirected edges with same branch_length pass validation."""
        # Internal node must have degree >= 3, and node 2 needs degree >= 3 too
        net = SemiDirectedPhyNetwork(
            undirected_edges=[
                {'u': 1, 'v': 2, 'branch_length': 0.5},
                {'u': 1, 'v': 2, 'branch_length': 0.5},
                (1, 3), (2, 4), (2, 5)  # Additional edges to satisfy degree constraints
            ],
            nodes=[(3, {'label': 'A'}), (4, {'label': 'B'}), (5, {'label': 'C'})]
        )
        # Should not raise
        net._validate_branchlength_constraints()

    def test_parallel_directed_edges_same_branchlength(self) -> None:
        """Test that parallel directed edges with same branch_length pass validation."""
        # Parallel directed edges from same source to same hybrid (with explicit keys)
        # Need single source component - all tree nodes must be connected via undirected edges
        # Gamma: 0.33 + 0.33 + 0.34 = 1.0
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {'u': 5, 'v': 4, 'key': 0, 'gamma': 0.33, 'branch_length': 0.5},
                {'u': 5, 'v': 4, 'key': 1, 'gamma': 0.33, 'branch_length': 0.5},
                {'u': 6, 'v': 4, 'gamma': 0.34}  # Another source to hybrid
            ],
            undirected_edges=[(4, 1), (5, 2), (5, 3), (5, 6), (6, 7), (6, 8)],  # Connect 5 and 6
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'}), (8, {'label': 'E'})]
        )
        # Should not raise
        net._validate_branchlength_constraints()

