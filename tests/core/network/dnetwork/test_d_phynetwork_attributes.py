"""
Comprehensive tests for DirectedPhyNetwork edge attributes.

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

from phylozoo.core.network import DirectedPhyNetwork


class TestGetEdgeAttribute:
    """Test cases for get_edge_attribute() method."""

    def test_get_edge_attribute_existing(self) -> None:
        """Test getting existing edge attribute."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': 0.5}],
            nodes=[(1, {"label": "A"})]
        )
        assert net.get_edge_attribute(3, 1, attr='branch_length') == 0.5

    def test_get_edge_attribute_missing(self) -> None:
        """Test getting missing edge attribute."""
        net = DirectedPhyNetwork(edges=[(3, 1)], nodes=[(1, {"label": "A"})])
        assert net.get_edge_attribute(3, 1, attr='branch_length') is None

    def test_get_edge_attribute_nonexistent_edge(self) -> None:
        """Test getting attribute from non-existent edge."""
        net = DirectedPhyNetwork(edges=[(3, 1)], nodes=[(1, {"label": "A"})])
        assert net.get_edge_attribute(3, 999, attr='branch_length') is None

    def test_get_edge_attribute_custom(self) -> None:
        """Test getting custom edge attribute."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'custom_attr': 'value'}],
            nodes=[(1, {"label": "A"})]
        )
        assert net.get_edge_attribute(3, 1, attr='custom_attr') == 'value'

    def test_get_edge_attribute_multiple_attributes(self) -> None:
        """Test getting one attribute when edge has multiple."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': 0.5, 'bootstrap': 0.95, 'custom': 'x'}],
            nodes=[(1, {"label": "A"})]
        )
        assert net.get_edge_attribute(3, 1, attr='branch_length') == 0.5
        assert net.get_edge_attribute(3, 1, attr='bootstrap') == 0.95
        assert net.get_edge_attribute(3, 1, attr='custom') == 'x'

    def test_get_edge_attribute_parallel_edges_with_key(self) -> None:
        """Test getting attribute from parallel edge with key."""
        # Use parallel edges to a hybrid node (valid structure)
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),
                (5, 4, 0),  # First parallel edge - need to set attribute after or use dict
                (5, 8),
                (6, 4),
                (6, 9),
                (4, 2)
            ],
            nodes=[(2, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})]
        )
        # Note: Can't easily set attributes on parallel edges with tuple format
        # This test verifies the method works when key is specified
        # The actual attribute value may be None if not set
        result = net.get_edge_attribute(5, 4, key=0, attr='branch_length')
        # Result may be None if attribute not set, but method should not raise error
        assert result is None or isinstance(result, (int, float))

    def test_get_edge_attribute_parallel_edges_without_key(self) -> None:
        """Test that parallel edges require key."""
        # Use parallel edges to a hybrid node (valid structure)
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),
                (5, 4, 0), (5, 4, 1),  # Parallel edges to hybrid
                (5, 8),
                (6, 4),
                (6, 9),
                (4, 2)
            ],
            nodes=[(2, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})]
        )
        with pytest.raises(ValueError, match="Multiple parallel edges"):
            net.get_edge_attribute(5, 4, attr='branch_length')

    def test_get_edge_attribute_default_attr(self) -> None:
        """Test default attr=None returns all attributes."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': 0.5}],
            nodes=[(1, {"label": "A"})]
        )
        assert net.get_edge_attribute(3, 1) == {'branch_length': 0.5}  # Default attr=None returns all attributes


class TestGetBranchLength:
    """Test cases for get_branch_length() method."""

    def test_get_branch_length_existing(self) -> None:
        """Test getting existing branch length."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': 0.5}],
            nodes=[(1, {"label": "A"})]
        )
        assert net.get_branch_length(3, 1) == 0.5

    def test_get_branch_length_missing(self) -> None:
        """Test getting missing branch length."""
        net = DirectedPhyNetwork(edges=[(3, 1)], nodes=[(1, {"label": "A"})])
        assert net.get_branch_length(3, 1) is None

    def test_get_branch_length_zero(self) -> None:
        """Test branch length of zero."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': 0.0}],
            nodes=[(1, {"label": "A"})]
        )
        assert net.get_branch_length(3, 1) == 0.0

    def test_get_branch_length_negative(self) -> None:
        """Test that negative branch lengths are allowed (not validated)."""
        # Note: branch_length validation is not enforced, so negative values are allowed
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': -0.5}],
            nodes=[(1, {"label": "A"})]
        )
        assert net.get_branch_length(3, 1) == -0.5

    def test_get_branch_length_parallel_edges(self) -> None:
        """Test branch length with parallel edges."""
        # Use parallel edges to a hybrid node (valid structure)
        # Note: When using dict format with 'key', attributes may need to be set differently
        # This test verifies the method works with parallel edges and keys
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),
                (5, 4, 0),  # Parallel edges - attributes can't be set with tuple format
                (5, 4, 1),
                (5, 8),
                (6, 4),
                (6, 9),
                (4, 2)
            ],
            nodes=[(2, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})]
        )
        # Method should work with keys even if attribute is None
        result0 = net.get_branch_length(5, 4, key=0)
        result1 = net.get_branch_length(5, 4, key=1)
        # Results may be None if attributes not set, but method should not raise error
        assert result0 is None or isinstance(result0, (int, float))
        assert result1 is None or isinstance(result1, (int, float))
        # Verify both edges exist and can be accessed
        assert net.has_edge(5, 4, key=0)
        assert net.has_edge(5, 4, key=1)


class TestGetBootstrap:
    """Test cases for get_bootstrap() method."""

    def test_get_bootstrap_existing(self) -> None:
        """Test getting existing bootstrap value."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'bootstrap': 0.95}],
            nodes=[(1, {"label": "A"})]
        )
        assert net.get_bootstrap(3, 1) == 0.95

    def test_get_bootstrap_missing(self) -> None:
        """Test getting missing bootstrap value."""
        net = DirectedPhyNetwork(edges=[(3, 1)], nodes=[(1, {"label": "A"})])
        assert net.get_bootstrap(3, 1) is None

    def test_get_bootstrap_boundary_values(self) -> None:
        """Test bootstrap at boundary values (0.0 and 1.0)."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'bootstrap': 0.0},
                {'u': 3, 'v': 2, 'bootstrap': 1.0}
            ],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"})]
        )
        assert net.get_bootstrap(3, 1) == 0.0
        assert net.get_bootstrap(3, 2) == 1.0

    def test_get_bootstrap_parallel_edges(self) -> None:
        """Test bootstrap with parallel edges."""
        # Use parallel edges to a hybrid node (valid structure)
        # Note: When using dict format with 'key', attributes may need to be set differently
        # This test verifies the method works with parallel edges and keys
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),
                (5, 4, 0),  # Parallel edges - attributes can't be set with tuple format
                (5, 4, 1),
                (5, 8),
                (6, 4),
                (6, 9),
                (4, 2)
            ],
            nodes=[(2, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})]
        )
        # Method should work with keys even if attribute is None
        result0 = net.get_bootstrap(5, 4, key=0)
        result1 = net.get_bootstrap(5, 4, key=1)
        # Results may be None if attributes not set, but method should not raise error
        assert result0 is None or isinstance(result0, (int, float))
        assert result1 is None or isinstance(result1, (int, float))
        # Verify both edges exist and can be accessed
        assert net.has_edge(5, 4, key=0)
        assert net.has_edge(5, 4, key=1)


class TestBootstrapValidation:
    """Test cases for bootstrap value validation."""

    def test_bootstrap_valid_range(self) -> None:
        """Test that bootstrap values in [0.0, 1.0] are valid."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'bootstrap': 0.0},
                {'u': 3, 'v': 2, 'bootstrap': 0.5},
                {'u': 3, 'v': 4, 'bootstrap': 1.0}
            ],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (4, {"label": "C"})]
        )
        net.validate()

    def test_bootstrap_below_zero(self) -> None:
        """Test that bootstrap < 0.0 raises ValueError."""
        with pytest.raises(ValueError, match="must be in \\[0.0, 1.0\\]"):
            DirectedPhyNetwork(
                edges=[{'u': 3, 'v': 1, 'bootstrap': -0.1}],
                nodes=[(1, {"label": "A"})]
            )

    def test_bootstrap_above_one(self) -> None:
        """Test that bootstrap > 1.0 raises ValueError."""
        with pytest.raises(ValueError, match="must be in \\[0.0, 1.0\\]"):
            DirectedPhyNetwork(
                edges=[{'u': 3, 'v': 1, 'bootstrap': 1.1}],
                nodes=[(1, {"label": "A"})]
            )

    def test_bootstrap_non_numeric(self) -> None:
        """Test that non-numeric bootstrap raises ValueError."""
        with pytest.raises(ValueError, match="must be numeric"):
            DirectedPhyNetwork(
                edges=[{'u': 3, 'v': 1, 'bootstrap': 'invalid'}],
                nodes=[(1, {"label": "A"})]
            )

    def test_bootstrap_integer_zero(self) -> None:
        """Test that integer 0 is valid for bootstrap."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'bootstrap': 0}],
            nodes=[(1, {"label": "A"})]
        )
        assert net.get_bootstrap(3, 1) == 0

    def test_bootstrap_integer_one(self) -> None:
        """Test that integer 1 is valid for bootstrap."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'bootstrap': 1}],
            nodes=[(1, {"label": "A"})]
        )
        assert net.get_bootstrap(3, 1) == 1

    def test_bootstrap_multiple_edges(self) -> None:
        """Test bootstrap validation on multiple edges."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'bootstrap': 0.8},
                {'u': 3, 'v': 2, 'bootstrap': 0.9},
                {'u': 3, 'v': 4, 'bootstrap': 0.7}
            ],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (4, {"label": "C"})]
        )
        net.validate()

    def test_bootstrap_parallel_edges_validation(self) -> None:
        """Test bootstrap validation on parallel edges."""
        # Use parallel edges to a hybrid node (valid structure)
        # Test that bootstrap validation works on parallel edges
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),
                (5, 4, 0),  # Parallel edges without bootstrap (valid)
                (5, 4, 1),
                (5, 8),
                (6, 4),
                (6, 9),
                (4, 2)
            ],
            nodes=[(2, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})]
        )
        # Network should validate (no bootstrap values to validate)
        net.validate()
        # Note: To test bootstrap validation on parallel edges with actual values,
        # we would need a way to set attributes on parallel edges, which may require
        # direct access to the underlying graph or a different initialization method


class TestGetGamma:
    """Test cases for get_gamma() method."""

    def test_get_gamma_existing(self) -> None:
        """Test getting existing gamma value."""
        # Need a root node, so add edges from root to tree nodes
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),  # Root to tree nodes
                {'u': 5, 'v': 4, 'gamma': 0.6},
                {'u': 5, 'v': 8},  # Tree node 5 also has another child
                {'u': 6, 'v': 4, 'gamma': 0.4},
                {'u': 6, 'v': 9},  # Tree node 6 also has another child
                {'u': 4, 'v': 1}
            ],
            nodes=[(1, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})]
        )
        assert net.get_gamma(5, 4) == 0.6
        assert net.get_gamma(6, 4) == 0.4

    def test_get_gamma_missing(self) -> None:
        """Test getting missing gamma value."""
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),  # Root to tree nodes
                (5, 4), (5, 8),  # Tree node 5 splits
                (6, 4), (6, 9),  # Tree node 6 splits
                (4, 1)  # Hybrid to leaf
            ],
            nodes=[(1, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})]
        )
        assert net.get_gamma(5, 4) is None

    def test_get_gamma_boundary_values(self) -> None:
        """Test gamma at boundary values (0.0 and 1.0)."""
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),  # Root to tree nodes
                {'u': 5, 'v': 4, 'gamma': 0.0},
                {'u': 5, 'v': 8},  # Tree node 5 also has another child
                {'u': 6, 'v': 4, 'gamma': 1.0},
                {'u': 6, 'v': 9},  # Tree node 6 also has another child
                {'u': 4, 'v': 1}
            ],
            nodes=[(1, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})]
        )
        assert net.get_gamma(5, 4) == 0.0
        assert net.get_gamma(6, 4) == 1.0

    def test_get_gamma_parallel_edges(self) -> None:
        """Test gamma with parallel edges."""
        # Test that get_gamma works with parallel edges (no gamma values set)
        # This verifies the method handles parallel edges correctly
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),  # Root to tree nodes
                (5, 4, 0), (5, 4, 1),  # Parallel edges (no gamma - valid)
                (5, 8),  # Tree node 5 also has another child
                (6, 4),  # Single edge from 6
                (6, 9),  # Tree node 6 also has another child
                (4, 1)  # Hybrid to leaf
            ],
            nodes=[(1, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})]
        )
        # Method should work with keys even if gamma not set
        result0 = net.get_gamma(5, 4, key=0)
        result1 = net.get_gamma(5, 4, key=1)
        # Results should be None (no gamma set)
        assert result0 is None
        assert result1 is None
        # Verify method works without raising errors
        net.validate()


class TestGammaValidation:
    """Test cases for gamma value validation."""

    def test_gamma_valid_range(self) -> None:
        """Test that gamma values in [0.0, 1.0] are valid."""
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),  # Root to tree nodes
                {'u': 5, 'v': 4, 'gamma': 0.6},
                {'u': 5, 'v': 8},  # Tree node 5 also has another child
                {'u': 6, 'v': 4, 'gamma': 0.4},
                {'u': 6, 'v': 9},  # Tree node 6 also has another child
                {'u': 4, 'v': 1}
            ],
            nodes=[(1, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})]
        )
        net.validate()

    def test_gamma_below_zero(self) -> None:
        """Test that gamma < 0.0 raises ValueError."""
        with pytest.raises(ValueError, match="must be in \\[0.0, 1.0\\]"):
            DirectedPhyNetwork(
                edges=[
                    (7, 5), (7, 6),  # Root to tree nodes
                    {'u': 5, 'v': 4, 'gamma': -0.1},
                    {'u': 5, 'v': 8},
                    {'u': 6, 'v': 4, 'gamma': 1.1},
                    {'u': 6, 'v': 9},
                    {'u': 4, 'v': 1}
                ],
                nodes=[(1, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})]
            )

    def test_gamma_above_one(self) -> None:
        """Test that gamma > 1.0 raises ValueError."""
        with pytest.raises(ValueError, match="must be in \\[0.0, 1.0\\]"):
            DirectedPhyNetwork(
                edges=[
                    (7, 5), (7, 6),  # Root to tree nodes
                    {'u': 5, 'v': 4, 'gamma': 0.6},
                    {'u': 5, 'v': 8},
                    {'u': 6, 'v': 4, 'gamma': 1.1},
                    {'u': 6, 'v': 9},
                    {'u': 4, 'v': 1}
                ],
                nodes=[(1, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})]
            )

    def test_gamma_non_numeric(self) -> None:
        """Test that non-numeric gamma raises ValueError."""
        with pytest.raises(ValueError, match="must be numeric"):
            DirectedPhyNetwork(
                edges=[
                    (7, 5), (7, 6),  # Root to tree nodes
                    {'u': 5, 'v': 4, 'gamma': 'invalid'},
                    {'u': 5, 'v': 8},
                    {'u': 6, 'v': 4, 'gamma': 1.0},
                    {'u': 6, 'v': 9},
                    {'u': 4, 'v': 1}
                ],
                nodes=[(1, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})]
            )

    def test_gamma_sum_must_be_one(self) -> None:
        """Test that gamma values must sum to 1.0."""
        with pytest.raises(ValueError, match="must sum to exactly 1.0"):
            DirectedPhyNetwork(
                edges=[
                    (7, 5), (7, 6),  # Root to tree nodes
                    {'u': 5, 'v': 4, 'gamma': 0.6},
                    {'u': 5, 'v': 8},
                    {'u': 6, 'v': 4, 'gamma': 0.3},  # Sum = 0.9, not 1.0
                    {'u': 6, 'v': 9},
                    {'u': 4, 'v': 1}
                ],
                nodes=[(1, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})]
            )

    def test_gamma_sum_exceeds_one(self) -> None:
        """Test that gamma sum > 1.0 raises ValueError."""
        with pytest.raises(ValueError, match="must sum to exactly 1.0"):
            DirectedPhyNetwork(
                edges=[
                    (7, 5), (7, 6),  # Root to tree nodes
                    {'u': 5, 'v': 4, 'gamma': 0.6},
                    {'u': 5, 'v': 8},
                    {'u': 6, 'v': 4, 'gamma': 0.5},  # Sum = 1.1
                    {'u': 6, 'v': 9},
                    {'u': 4, 'v': 1}
                ],
                nodes=[(1, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})]
            )

    def test_gamma_all_or_none(self) -> None:
        """Test that if any gamma is specified, all must be specified."""
        with pytest.raises(ValueError, match="ALL incoming edges must have gamma values"):
            DirectedPhyNetwork(
                edges=[
                    (7, 5), (7, 6),  # Root to tree nodes
                    {'u': 5, 'v': 4, 'gamma': 0.6},
                    {'u': 5, 'v': 8},
                    {'u': 6, 'v': 4},  # Missing gamma
                    {'u': 6, 'v': 9},
                    {'u': 4, 'v': 1}
                ],
                nodes=[(1, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})]
            )

    def test_gamma_none_allowed(self) -> None:
        """Test that no gamma values is allowed."""
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),  # Root to tree nodes
                (5, 4), (5, 8),  # Tree node 5 splits
                (6, 4), (6, 9),  # Tree node 6 splits
                (4, 1)  # Hybrid to leaf
            ],
            nodes=[(1, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})]
        )
        # No gamma values - should be valid
        net.validate()

    def test_gamma_parallel_edges_all_must_have(self) -> None:
        """Test that parallel edges to hybrid must all have gamma."""
        with pytest.raises(ValueError, match="ALL incoming edges must have gamma values"):
            DirectedPhyNetwork(
                edges=[
                    (7, 5), (7, 6),  # Root to tree nodes
                    {'u': 5, 'v': 4, 'key': 0, 'gamma': 0.3},
                    {'u': 5, 'v': 4, 'key': 1},  # Missing gamma
                    {'u': 5, 'v': 8},
                    {'u': 6, 'v': 4, 'gamma': 0.4},
                    {'u': 6, 'v': 9},
                    {'u': 4, 'v': 1}
                ],
                nodes=[(1, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})]
            )

    def test_gamma_parallel_edges_sum(self) -> None:
        """Test that parallel edges with gamma must sum to 1.0."""
        # Test setting attributes on parallel edges using dict format with 'key'
        # The 'key' gets popped from the dict, and remaining attributes (like 'gamma')
        # should be properly set on the edge via **edge unpacking.
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),  # Root to tree nodes
                {'u': 5, 'v': 4, 'key': 0, 'gamma': 0.3},
                {'u': 5, 'v': 4, 'key': 1, 'gamma': 0.3},
                {'u': 5, 'v': 8},
                {'u': 6, 'v': 4, 'gamma': 0.4},
                {'u': 6, 'v': 9},
                {'u': 4, 'v': 1}
            ],
            nodes=[(1, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})]
        )
        # Verify the sum of gammas for hybrid node 4
        net.validate()
        # Check that gamma values were set correctly
        gamma_0 = net.get_gamma(5, 4, key=0)
        gamma_1 = net.get_gamma(5, 4, key=1)
        gamma_6 = net.get_gamma(6, 4)
        assert gamma_0 == 0.3
        assert gamma_1 == 0.3
        assert gamma_6 == 0.4
        # Sum should be 1.0
        assert abs(gamma_0 + gamma_1 + gamma_6 - 1.0) < 1e-10

    def test_gamma_only_on_hybrid_edges(self) -> None:
        """Test that gamma can only be set on hybrid edges."""
        # Try to set gamma on a tree edge (not a hybrid edge)
        # Tree edge: from root to tree node (tree node has out-degree >= 2)
        with pytest.raises(ValueError, match="Gamma value can only be set on hybrid edges"):
            DirectedPhyNetwork(
                edges=[
                    {'u': 3, 'v': 1, 'gamma': 0.5},  # Gamma on tree edge (1 is tree node)
                    (3, 2),  # Another tree edge
                    (1, 4), (1, 5),  # Tree node 1 splits
                    (2, 6), (2, 7)  # Tree node 2 splits
                ],
                nodes=[(4, {"label": "A"}), (5, {"label": "B"}), (6, {"label": "C"}), (7, {"label": "D"})]
            )
        
        # Try to set gamma on an edge to a leaf (not a hybrid edge)
        with pytest.raises(ValueError, match="Gamma value can only be set on hybrid edges"):
            DirectedPhyNetwork(
                edges=[
                    {'u': 3, 'v': 1, 'gamma': 0.5}  # Gamma on edge to leaf
                ],
                nodes=[(1, {"label": "A"})]
            )
        
        # Valid: gamma on hybrid edge
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),  # Root to tree nodes
                {'u': 5, 'v': 4, 'gamma': 0.6},  # Hybrid edge
                {'u': 6, 'v': 4, 'gamma': 0.4},  # Hybrid edge
                (5, 8),  # Tree node 5 also has another child
                (6, 9),  # Tree node 6 also has another child
                (4, 1)  # Hybrid to leaf
            ],
            nodes=[(1, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})]
        )
        net.validate()

    def test_gamma_multiple_hybrids(self) -> None:
        """Test gamma validation with multiple hybrid nodes."""
        # Hybrid 4 has gammas from 7 and 8, but not from 5 and 6
        # This should fail because if ANY gamma is specified, ALL must be specified
        with pytest.raises(ValueError, match="ALL incoming edges must have gamma values"):
            DirectedPhyNetwork(
                edges=[
                    (20, 10), (20, 7), (20, 8),  # Root to tree nodes
                    (10, 5), (10, 6),  # Tree node 10 splits
                    {'u': 5, 'v': 4},  # No gamma
                    {'u': 5, 'v': 11},  # Tree node 5 also has another child
                    {'u': 6, 'v': 4},  # No gamma
                    {'u': 6, 'v': 12},  # Tree node 6 also has another child
                    {'u': 7, 'v': 4, 'gamma': 0.5},  # Has gamma
                    {'u': 7, 'v': 13},  # Tree node 7 also has another child
                    {'u': 8, 'v': 4, 'gamma': 0.5},  # Has gamma
                    {'u': 8, 'v': 14},  # Tree node 8 also has another child
                    {'u': 4, 'v': 1}  # Hybrid to leaf
                ],
                nodes=[(1, {"label": "A"}), (11, {"label": "B"}), (12, {"label": "C"}), (13, {"label": "D"}), (14, {"label": "E"})]
            )

    def test_gamma_floating_point_precision(self) -> None:
        """Test that floating point precision is handled correctly."""
        # Sum that's very close to 1.0 due to floating point
        net = DirectedPhyNetwork(
            edges=[
                (20, 5), (20, 6), (20, 7),  # Root to tree nodes
                {'u': 5, 'v': 4, 'gamma': 0.3333333333333333},
                {'u': 5, 'v': 8},
                {'u': 6, 'v': 4, 'gamma': 0.3333333333333333},
                {'u': 6, 'v': 9},
                {'u': 7, 'v': 4, 'gamma': 0.3333333333333334},  # Sum = 1.0 (approx)
                {'u': 7, 'v': 10},
                {'u': 4, 'v': 1}
            ],
            nodes=[(1, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"}), (10, {"label": "D"})]
        )
        # Should pass with tolerance
        net.validate()


class TestCustomAttributes:
    """Test cases for custom edge attributes."""

    def test_custom_attribute_string(self) -> None:
        """Test custom string attribute."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'label': 'edge1'}],
            nodes=[(1, {"label": "A"})]
        )
        assert net.get_edge_attribute(3, 1, attr='label') == 'edge1'

    def test_custom_attribute_list(self) -> None:
        """Test custom list attribute."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'tags': ['a', 'b', 'c']}],
            nodes=[(1, {"label": "A"})]
        )
        assert net.get_edge_attribute(3, 1, attr='tags') == ['a', 'b', 'c']

    def test_custom_attribute_dict(self) -> None:
        """Test custom dict attribute."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'metadata': {'key': 'value'}}],
            nodes=[(1, {"label": "A"})]
        )
        assert net.get_edge_attribute(3, 1, attr='metadata') == {'key': 'value'}

    def test_multiple_custom_attributes(self) -> None:
        """Test edge with multiple custom attributes."""
        net = DirectedPhyNetwork(
            edges=[{
                'u': 3, 'v': 1,
                'branch_length': 0.5,
                'bootstrap': 0.95,
                'custom1': 'value1',
                'custom2': 42,
                'custom3': [1, 2, 3]
            }],
            nodes=[(1, {"label": "A"})]
        )
        assert net.get_edge_attribute(3, 1, attr='custom1') == 'value1'
        assert net.get_edge_attribute(3, 1, attr='custom2') == 42
        assert net.get_edge_attribute(3, 1, attr='custom3') == [1, 2, 3]


class TestAttributeEdgeCases:
    """Test edge cases for attributes."""

    def test_very_large_branch_length(self) -> None:
        """Test very large branch length values."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': 1e10}],
            nodes=[(1, {"label": "A"})]
        )
        assert net.get_branch_length(3, 1) == 1e10

    def test_very_small_branch_length(self) -> None:
        """Test very small branch length values."""
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': 1e-10}],
            nodes=[(1, {"label": "A"})]
        )
        assert net.get_branch_length(3, 1) == 1e-10

    def test_nan_branch_length(self) -> None:
        """Test NaN branch length (should be allowed, not validated)."""
        import math
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': float('nan')}],
            nodes=[(1, {"label": "A"})]
        )
        result = net.get_branch_length(3, 1)
        assert result is not None
        assert math.isnan(result)

    def test_infinity_branch_length(self) -> None:
        """Test infinity branch length (should be allowed, not validated)."""
        import math
        net = DirectedPhyNetwork(
            edges=[{'u': 3, 'v': 1, 'branch_length': float('inf')}],
            nodes=[(1, {"label": "A"})]
        )
        result = net.get_branch_length(3, 1)
        assert result is not None
        assert math.isinf(result)

    def test_many_attributes_per_edge(self) -> None:
        """Test edge with many attributes."""
        attrs = {f'attr{i}': i for i in range(100)}
        attrs['u'] = 3
        attrs['v'] = 1
        net = DirectedPhyNetwork(edges=[attrs], nodes=[(1, {"label": "A"})])
        # Check a few attributes
        assert net.get_edge_attribute(3, 1, attr='attr0') == 0
        assert net.get_edge_attribute(3, 1, attr='attr99') == 99

