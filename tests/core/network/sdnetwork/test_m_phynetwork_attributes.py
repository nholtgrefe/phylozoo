"""
Comprehensive tests for MixedPhyNetwork edge attributes.

This module tests all aspects of edge attribute handling including:
- get_edge_attribute() method
- get_branch_length() method
- get_bootstrap() method
- get_gamma() method
- Bootstrap validation
- Gamma validation
- Custom attributes
- Parallel edges with attributes
"""

import math
import warnings

import pytest

from phylozoo.core.network.sdnetwork import MixedPhyNetwork
from tests.core.network.sdnetwork.conftest import expect_mixed_network_warning


class TestGetEdgeAttribute:
    """Test cases for get_edge_attribute() method."""

    def test_get_edge_attribute_undirected_existing(self) -> None:
        """Test getting existing edge attribute from undirected edge."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[{'u': 3, 'v': 1, 'branch_length': 0.5}, (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.get_edge_attribute(3, 1, attr='branch_length') == 0.5

    def test_get_edge_attribute_directed_existing(self) -> None:
        """Test getting existing edge attribute from directed edge (hybrid edge)."""
        # Hybrid node 4 has directed edges from 3 and 5, undirected edges to leaves
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[{'u': 3, 'v': 4, 'branch_length': 0.5}, (5, 4)],
            undirected_edges=[(4, 1), (3, 2), (3, 6), (5, 7), (5, 8)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (6, {'label': 'C'}), (7, {'label': 'D'}), (8, {'label': 'E'})]
            )
        assert net.get_edge_attribute(3, 4, attr='branch_length') == 0.5

    def test_get_edge_attribute_missing(self) -> None:
        """Test getting missing edge attribute."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.get_edge_attribute(3, 1, attr='branch_length') is None

    def test_get_edge_attribute_nonexistent_edge(self) -> None:
        """Test getting attribute from non-existent edge."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.get_edge_attribute(3, 999, attr='branch_length') is None

    def test_get_edge_attribute_custom(self) -> None:
        """Test getting custom edge attribute."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[{'u': 3, 'v': 1, 'custom_attr': 'value'}, (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.get_edge_attribute(3, 1, attr='custom_attr') == 'value'

    def test_get_edge_attribute_multiple_attributes(self) -> None:
        """Test getting one attribute when edge has multiple."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[{'u': 3, 'v': 1, 'branch_length': 0.5, 'bootstrap': 0.95, 'custom': 'x'}, (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.get_edge_attribute(3, 1, attr='branch_length') == 0.5
        assert net.get_edge_attribute(3, 1, attr='bootstrap') == 0.95
        assert net.get_edge_attribute(3, 1, attr='custom') == 'x'

    def test_get_edge_attribute_parallel_edges_with_key(self) -> None:
        """Test getting attribute from parallel edge with key."""
        # Parallel undirected edges between internal nodes 3 and 4
        # Both nodes need degree >= 3
        # Use no_validation to allow different branch_length values for testing
        from phylozoo.utils.validation import no_validation
        with no_validation():
            net = MixedPhyNetwork(
            undirected_edges=[
            {'u': 3, 'v': 4, 'key': 0, 'branch_length': 0.5},
            {'u': 3, 'v': 4, 'key': 1, 'branch_length': 0.7},
            (3, 1), (3, 2),  # Additional edges from 3
            (4, 5), (4, 6)   # Additional edges from 4
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
            )
        assert net.get_edge_attribute(3, 4, key=0, attr='branch_length') == 0.5
        assert net.get_edge_attribute(3, 4, key=1, attr='branch_length') == 0.7

    def test_get_edge_attribute_parallel_edges_without_key(self) -> None:
        """Test that parallel edges require key."""
        # Parallel edges between internal nodes 3 and 4
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[
            (3, 4, 0), (3, 4, 1),  # Parallel edges
            (3, 1), (3, 2),  # Additional edges from 3
            (4, 5), (4, 6)   # Additional edges from 4
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
            )
        with pytest.raises(ValueError, match="Multiple parallel"):
            net.get_edge_attribute(3, 4, attr='branch_length')

    def test_get_edge_attribute_directed_parameter(self) -> None:
        """Test get_edge_attribute with directed parameter."""
        # Hybrid node 5 with directed edges from 3 and 6
        # Node 2 needs degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[{'u': 3, 'v': 5, 'branch_length': 0.5}, (6, 5)],
            undirected_edges=[{'u': 2, 'v': 3, 'branch_length': 0.3}, (2, 4), (2, 11), (5, 1), (3, 7), (3, 8), (6, 9), (6, 10)],
            nodes=[(1, {'label': 'A'}), (4, {'label': 'B'}), (7, {'label': 'C'}), (8, {'label': 'D'}), (9, {'label': 'E'}), (10, {'label': 'F'}), (11, {'label': 'G'})]
            )
        # For directed edge (hybrid edge), should work
        assert net.get_edge_attribute(3, 5, attr='branch_length') == 0.5
        # For undirected edge
        assert net.get_edge_attribute(2, 3, attr='branch_length') == 0.3


class TestGetBranchLength:
    """Test cases for get_branch_length() method."""

    def test_get_branch_length_existing(self) -> None:
        """Test getting existing branch length."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[{'u': 3, 'v': 1, 'branch_length': 0.5}, (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.get_branch_length(3, 1) == 0.5

    def test_get_branch_length_missing(self) -> None:
        """Test getting missing branch length."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.get_branch_length(3, 1) is None

    def test_get_branch_length_parallel_edges(self) -> None:
        """Test getting branch length from parallel edges."""
        # Parallel edges between internal nodes
        # Use no_validation to allow different branch_length values for testing
        from phylozoo.utils.validation import no_validation
        with no_validation():
            net = MixedPhyNetwork(
            undirected_edges=[
            {'u': 3, 'v': 4, 'key': 0, 'branch_length': 0.5},
            {'u': 3, 'v': 4, 'key': 1, 'branch_length': 0.7},
            (3, 1), (3, 2),  # Additional edges from 3
            (4, 5), (4, 6)   # Additional edges from 4
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
            )
        assert net.get_branch_length(3, 4, key=0) == 0.5
        assert net.get_branch_length(3, 4, key=1) == 0.7

    def test_get_branch_length_directed_edge(self) -> None:
        """Test getting branch length from directed edge (hybrid edge)."""
        # Node 3 needs degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[{'u': 3, 'v': 4, 'branch_length': 0.5}, (5, 4)],
            undirected_edges=[(4, 1), (3, 2), (3, 7), (5, 6), (5, 8)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (6, {'label': 'C'}), (7, {'label': 'D'}), (8, {'label': 'E'})]
            )
        assert net.get_branch_length(3, 4) == 0.5


class TestGetBootstrap:
    """Test cases for get_bootstrap() method."""

    def test_get_bootstrap_existing(self) -> None:
        """Test getting existing bootstrap value."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[{'u': 3, 'v': 1, 'bootstrap': 0.95}, (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.get_bootstrap(3, 1) == 0.95

    def test_get_bootstrap_missing(self) -> None:
        """Test getting missing bootstrap value."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.get_bootstrap(3, 1) is None

    def test_get_bootstrap_parallel_edges(self) -> None:
        """Test getting bootstrap from parallel edges."""
        # Parallel edges between internal nodes
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[
            {'u': 3, 'v': 4, 'key': 0, 'bootstrap': 0.95},
            {'u': 3, 'v': 4, 'key': 1, 'bootstrap': 0.87},
            (3, 1), (3, 2),  # Additional edges from 3
            (4, 5), (4, 6)   # Additional edges from 4
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
            )
        assert net.get_bootstrap(3, 4, key=0) == 0.95
        assert net.get_bootstrap(3, 4, key=1) == 0.87

    def test_get_bootstrap_directed_edge(self) -> None:
        """Test getting bootstrap from directed edge (hybrid edge)."""
        # Node 3 needs degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[{'u': 3, 'v': 4, 'bootstrap': 0.95}, (5, 4)],
            undirected_edges=[(4, 1), (3, 2), (3, 7), (5, 6), (5, 8)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (6, {'label': 'C'}), (7, {'label': 'D'}), (8, {'label': 'E'})]
            )
        assert net.get_bootstrap(3, 4) == 0.95


class TestGetGamma:
    """Test cases for get_gamma() method."""

    def test_get_gamma_existing(self) -> None:
        """Test getting existing gamma value from hybrid edge."""
        # Hybrid node 4: indegree 2, total_degree must be 3 (only 1 outgoing)
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[
            {'u': 5, 'v': 4, 'gamma': 0.6},
            {'u': 6, 'v': 4, 'gamma': 0.4}
            ],
            undirected_edges=[(4, 1), (5, 7), (5, 8), (6, 9), (6, 10)],
            nodes=[(1, {'label': 'A'}), (7, {'label': 'B'}), (8, {'label': 'C'}), (9, {'label': 'D'}), (10, {'label': 'E'})]
            )
        assert net.get_gamma(5, 4) == 0.6
        assert net.get_gamma(6, 4) == 0.4

    def test_get_gamma_missing(self) -> None:
        """Test getting missing gamma value."""
        # Hybrid node 4: indegree 2, total_degree must be 3
        # Node 5 needs degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 1), (5, 7), (5, 8), (6, 9), (6, 10)],
            nodes=[(1, {'label': 'A'}), (7, {'label': 'B'}), (8, {'label': 'C'}), (9, {'label': 'D'}), (10, {'label': 'E'})]
            )
        assert net.get_gamma(5, 4) is None
        assert net.get_gamma(6, 4) is None

    def test_get_gamma_parallel_edges(self) -> None:
        """Test getting gamma from parallel hybrid edges."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[
            {'u': 5, 'v': 4, 'key': 0, 'gamma': 0.4},
            {'u': 5, 'v': 4, 'key': 1, 'gamma': 0.2},
            {'u': 6, 'v': 4, 'gamma': 0.4}
            ],
            undirected_edges=[(4, 1), (5, 8), (5, 9), (6, 10), (6, 11)],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        assert net.get_gamma(5, 4, key=0) == 0.4
        assert net.get_gamma(5, 4, key=1) == 0.2
        assert net.get_gamma(6, 4) == 0.4

    def test_get_gamma_undirected_edge_raises_error(self) -> None:
        """Test that getting gamma from undirected edge raises error."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        with pytest.raises(ValueError, match="cannot be set on undirected edges"):
            net.get_gamma(3, 1)

    def test_get_gamma_non_hybrid_directed_edge_raises_error(self) -> None:
        """Test that getting gamma from non-hybrid directed edge raises error."""
        # Test that get_gamma raises error when called on an undirected edge
        # Hybrid node 4: indegree 2, total_degree must be 3 (only 1 outgoing)
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],  # Node 4 is hybrid
            undirected_edges=[(4, 1), (3, 2), (3, 7), (3, 5), (5, 8), (5, 9), (6, 10), (6, 11)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (7, {'label': 'C'}), (8, {'label': 'D'}), (9, {'label': 'E'}), (10, {'label': 'F'}), (11, {'label': 'G'})]
            )
        # (3, 2) is undirected, so get_gamma should raise error
        with pytest.raises(ValueError, match="cannot be set on undirected edges"):
            net.get_gamma(3, 2)


class TestBootstrapValidation:
    """Test cases for bootstrap value validation."""

    def test_bootstrap_in_valid_range(self) -> None:
        """Test that bootstrap values in [0.0, 1.0] are valid."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
                undirected_edges=[
                    {'u': 3, 'v': 1, 'bootstrap': 0.0},
                    {'u': 3, 'v': 2, 'bootstrap': 0.5},
                    {'u': 3, 'v': 4, 'bootstrap': 1.0}
                ],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        with expect_mixed_network_warning():
            net.validate()

    def test_bootstrap_below_zero_raises_error(self) -> None:
        """Test that bootstrap < 0.0 raises error."""
        with expect_mixed_network_warning():
            with pytest.raises(ValueError, match="Bootstrap value.*must be in"):
                MixedPhyNetwork(
                    undirected_edges=[{'u': 3, 'v': 1, 'bootstrap': -0.1}, (3, 2), (3, 4)],
                    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
                )

    def test_bootstrap_above_one_raises_error(self) -> None:
        """Test that bootstrap > 1.0 raises error."""
        with expect_mixed_network_warning():
            with pytest.raises(ValueError, match="Bootstrap value.*must be in"):
                MixedPhyNetwork(
                    undirected_edges=[{'u': 3, 'v': 1, 'bootstrap': 1.1}, (3, 2), (3, 4)],
                    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
                )

    def test_bootstrap_nan_raises_error(self) -> None:
        """Test that NaN bootstrap raises error."""
        with expect_mixed_network_warning():
            with pytest.raises(ValueError, match="Bootstrap value.*must be in"):
                MixedPhyNetwork(
                    undirected_edges=[{'u': 3, 'v': 1, 'bootstrap': float('nan')}, (3, 2), (3, 4)],
                    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
                )

    def test_bootstrap_on_directed_edge(self) -> None:
        """Test bootstrap on directed edge (hybrid edge)."""
        # Node 3 needs degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[{'u': 3, 'v': 4, 'bootstrap': 0.95}, (5, 4)],
            undirected_edges=[(4, 1), (3, 2), (3, 7), (5, 6), (5, 8)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (6, {'label': 'C'}), (7, {'label': 'D'}), (8, {'label': 'E'})]
            )
        assert net.get_bootstrap(3, 4) == 0.95


class TestGammaValidation:
    """Test cases for gamma value validation."""

    def test_gamma_on_hybrid_edge_valid(self) -> None:
        """Test that gamma on hybrid edge is valid."""
        # Nodes 5 and 6 need degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
                directed_edges=[
                    {'u': 5, 'v': 4, 'gamma': 0.6},
                    {'u': 6, 'v': 4, 'gamma': 0.4}
                ],
                undirected_edges=[(4, 1), (5, 8), (5, 10), (6, 9), (6, 11)],
                nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        with expect_mixed_network_warning():
            net.validate()

    def test_gamma_sum_to_one(self) -> None:
        """Test that gamma values sum to 1.0."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
                directed_edges=[
                    {'u': 5, 'v': 4, 'gamma': 0.6},
                    {'u': 6, 'v': 4, 'gamma': 0.4}
                ],
                undirected_edges=[(4, 1), (5, 8), (5, 10), (6, 9), (6, 11)],
                nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        with expect_mixed_network_warning():
            net.validate()

    def test_gamma_sum_not_one_raises_error(self) -> None:
        """Test that gamma values not summing to 1.0 raise error."""
        with expect_mixed_network_warning():
            with pytest.raises(ValueError, match="sum to"):
                MixedPhyNetwork(
                    directed_edges=[
                        {'u': 5, 'v': 4, 'gamma': 0.6},
                        {'u': 6, 'v': 4, 'gamma': 0.3}  # Sum = 0.9
                    ],
                    undirected_edges=[(4, 1), (5, 8), (5, 10), (6, 9), (6, 11)],
                nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
                )

    def test_gamma_on_undirected_edge_raises_error(self) -> None:
        """Test that gamma on undirected edge raises error."""
        with expect_mixed_network_warning():
            with pytest.raises(ValueError, match="not a directed edge|cannot have gamma"):
                MixedPhyNetwork(
                    undirected_edges=[{'u': 3, 'v': 1, 'gamma': 0.5}, (3, 2), (3, 4)],
                    nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
                )

    def test_gamma_partial_raises_error(self) -> None:
        """Test that partial gamma (some edges have it, others don't) raises error."""
        with expect_mixed_network_warning():
            with pytest.raises(ValueError, match="some edges with gamma"):
                MixedPhyNetwork(
                    directed_edges=[
                        {'u': 5, 'v': 4, 'gamma': 0.6},
                        (6, 4)  # Missing gamma
                    ],
                    undirected_edges=[(4, 1), (5, 8), (5, 10), (6, 9), (6, 11)],
                    nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
                )

    def test_gamma_parallel_edges_sum(self) -> None:
        """Test that gamma on parallel edges sums correctly."""
        # Nodes 5 and 6 need degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
                directed_edges=[
                    {'u': 5, 'v': 4, 'key': 0, 'gamma': 0.3},
                    {'u': 5, 'v': 4, 'key': 1, 'gamma': 0.3},
                    {'u': 6, 'v': 4, 'gamma': 0.4}
                ],
                undirected_edges=[(4, 1), (5, 7), (5, 8), (6, 9), (6, 10)],
                nodes=[(1, {'label': 'A'}), (7, {'label': 'B'}), (8, {'label': 'C'}), (9, {'label': 'D'}), (10, {'label': 'E'})]
            )
        # Sum should be 0.3 + 0.3 + 0.4 = 1.0
        with expect_mixed_network_warning():
            net.validate()

    def test_gamma_in_valid_range(self) -> None:
        """Test that gamma values in [0.0, 1.0] are valid."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
                directed_edges=[
                    {'u': 5, 'v': 4, 'gamma': 0.0},
                    {'u': 6, 'v': 4, 'gamma': 1.0}
                ],
                undirected_edges=[(4, 1), (5, 8), (5, 10), (6, 9), (6, 11)],
                nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        with expect_mixed_network_warning():
            net.validate()

    def test_gamma_below_zero_raises_error(self) -> None:
        """Test that gamma < 0.0 raises error."""
        with expect_mixed_network_warning():
            with pytest.raises(ValueError, match="must be in"):
                MixedPhyNetwork(
                    directed_edges=[
                        {'u': 5, 'v': 4, 'gamma': -0.1},
                        {'u': 6, 'v': 4, 'gamma': 1.1}
                    ],
                    undirected_edges=[(4, 1), (5, 8), (5, 10), (6, 9), (6, 11)],
                nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
                )

    def test_gamma_above_one_raises_error(self) -> None:
        """Test that gamma > 1.0 raises error."""
        with expect_mixed_network_warning():
            with pytest.raises(ValueError, match="must be in"):
                MixedPhyNetwork(
                    directed_edges=[
                        {'u': 5, 'v': 4, 'gamma': 1.1},
                        {'u': 6, 'v': 4, 'gamma': -0.1}
                    ],
                    undirected_edges=[(4, 1), (5, 8), (5, 10), (6, 9), (6, 11)],
                nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
                )


class TestCustomAttributes:
    """Test cases for custom edge attributes."""

    def test_custom_attribute_preserved(self) -> None:
        """Test that custom attributes are preserved."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[{'u': 3, 'v': 1, 'custom': 'value', 'number': 42}, (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.get_edge_attribute(3, 1, attr='custom') == 'value'
        assert net.get_edge_attribute(3, 1, attr='number') == 42

    def test_multiple_custom_attributes(self) -> None:
        """Test multiple custom attributes on same edge."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[{
            'u': 3, 'v': 1,
            'branch_length': 0.5,
            'bootstrap': 0.95,
            'custom1': 'a',
            'custom2': 123,
            'custom3': [1, 2, 3]
            }, (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.get_edge_attribute(3, 1, attr='custom1') == 'a'
        assert net.get_edge_attribute(3, 1, attr='custom2') == 123
        assert net.get_edge_attribute(3, 1, attr='custom3') == [1, 2, 3]

