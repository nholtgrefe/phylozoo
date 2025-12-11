"""
Comprehensive tests for SemiDirectedPhyNetwork validation.

Tests cover:
- Valid semi-directed networks (should pass validation)
- Invalid networks caught by validation:
  - Disconnected networks
  - Invalid degree constraints
  - Invalid semi-directed constraints
  - Invalid bootstrap/gamma values
- Edge cases
"""

import math
import pytest

from phylozoo.core.network import SemiDirectedPhyNetwork


class TestValidSemiDirectedNetworks:
    """Test that valid semi-directed networks pass validation."""

    def test_ternary_tree(self) -> None:
        """Ternary tree should be valid."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(4, 1), (4, 2), (4, 3)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})],
        )
        assert net.validate() is True

    def test_larger_tree(self) -> None:
        """Larger tree with all internal nodes having degree >= 3 should be valid."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[
                (7, 5),
                (7, 6),
                (7, 10),
                (5, 1),
                (5, 2),
                (5, 8),
                (6, 3),
                (6, 4),
                (6, 9),
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (4, {'label': 'D'}), (8, {'label': 'E'}), (9, {'label': 'F'}), (10, {'label': 'G'})],
        )
        assert net.validate() is True

    def test_single_hybrid(self) -> None:
        """Network with single hybrid node should be valid."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {"u": 5, "v": 4, "gamma": 0.6},
                {"u": 6, "v": 4, "gamma": 0.4},
            ],
            undirected_edges=[
                (7, 5),
                (7, 6),
                (7, 10),
                (5, 8),
                (6, 9),
                (4, 1),
            ],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'})],
        )
        assert net.validate() is True

    def test_multiple_hybrids(self) -> None:
        """Network with multiple hybrid nodes should be valid."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {"u": 8, "v": 5, "gamma": 0.7},
                {"u": 9, "v": 5, "gamma": 0.3},
                {"u": 8, "v": 6, "gamma": 0.4},
                {"u": 9, "v": 6, "gamma": 0.6},
            ],
            undirected_edges=[
                (10, 8),
                (10, 9),
                (10, 13),
                (8, 11),
                (9, 12),
                (5, 1),
                (6, 2),
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (11, {'label': 'C'}), (12, {'label': 'D'}), (13, {'label': 'E'})],
        )
        assert net.validate() is True

    def test_star_tree(self) -> None:
        """Star tree should be valid."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(10, i) for i in range(1, 6)],
            nodes=[(i, {'label': f"Taxon{i}"}) for i in range(1, 6)],
        )
        assert net.validate() is True


class TestInvalidConnectivity:
    """Test that disconnected networks are caught."""

    def test_disconnected_two_components(self) -> None:
        """Network with two disconnected components should be invalid."""
        with pytest.raises(ValueError, match="Network is not connected"):
            SemiDirectedPhyNetwork(
                undirected_edges=[
                    (3, 1),
                    (3, 2),
                    (3, 7),
                    # Component 2 (disconnected)
                    (6, 4),
                    (6, 5),
                    (6, 8),
                ],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (5, {'label': 'D'}), (7, {'label': 'E'}), (8, {'label': 'F'})],
            )

    def test_disconnected_with_hybrid(self) -> None:
        """Disconnected network with hybrid should be invalid."""
        with pytest.raises(ValueError, match="Network is not connected"):
            SemiDirectedPhyNetwork(
                directed_edges=[
                    {"u": 5, "v": 4, "gamma": 0.6},
                    {"u": 6, "v": 4, "gamma": 0.4},
                ],
                undirected_edges=[
                    (7, 5),
                    (7, 6),
                    (7, 10),
                    (5, 11),
                    (6, 12),
                    (4, 1),
                    # Disconnected component
                    (15, 8),
                    (15, 9),
                    (15, 13),
                ],
                nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'}), (12, {'label': 'F'}), (13, {'label': 'G'})],
            )


class TestInvalidDegreeConstraints:
    """Test that invalid degree constraints are caught."""

    def test_internal_node_degree_2(self) -> None:
        """Internal node with degree 2 should be invalid."""
        with pytest.raises(
            ValueError, match="Internal node .* has degree 2, but all internal nodes must have degree >= 3"
        ):
            SemiDirectedPhyNetwork(
                undirected_edges=[
                    (3, 1),  # Node 3 has degree 2 (invalid)
                    (3, 2),
                ],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})],
            )

    def test_leaf_with_outgoing_directed_edge(self) -> None:
        """Leaf with outgoing directed edges should be invalid."""
        with pytest.raises(ValueError, match="Internal node .* has degree 2"):
            SemiDirectedPhyNetwork(
                directed_edges=[(1, 2)],  # 1 has outgoing edge
                undirected_edges=[(3, 1), (3, 2), (3, 4)],
                nodes=[(1, {'label': 'A'}), (4, {'label': 'B'})],  # 1 is marked as taxon but has outdegree > 0
            )

    def test_invalid_indegree_constraint(self) -> None:
        """Node with invalid indegree should be caught."""
        # Create a network where a node has indegree != 0 and indegree != total_degree - 1
        # For example: node with indegree=1, total_degree=4 (1 != 0 and 1 != 3)
        with pytest.raises(
            ValueError,
            match="Node .* has indegree .* and total degree .*\\. Each node must have indegree either 0 or total_degree-1",
        ):
            SemiDirectedPhyNetwork(
                directed_edges=[
                    (5, 4),  # Node 4 has indegree=1
                ],
                undirected_edges=[
                    (4, 1),  # Node 4 has 3 undirected edges
                    (4, 2),
                    (4, 3),
                    (5, 6),
                    (5, 7),
                ],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (6, {'label': 'D'}), (7, {'label': 'E'})],
            )


class TestInvalidBootstrap:
    """Test that invalid bootstrap values are caught."""

    def test_bootstrap_below_zero(self) -> None:
        """Bootstrap value below 0 should be invalid."""
        with pytest.raises(ValueError, match="Bootstrap value .* is .*, but must be in \\[0\\.0, 1\\.0\\]"):
            SemiDirectedPhyNetwork(
                undirected_edges=[
                    {"u": 3, "v": 1, "bootstrap": -0.1},
                    (3, 2),
                    (3, 4),
                ],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})],
            )

    def test_bootstrap_above_one(self) -> None:
        """Bootstrap value above 1 should be invalid."""
        with pytest.raises(ValueError, match="Bootstrap value .* is .*, but must be in \\[0\\.0, 1\\.0\\]"):
            SemiDirectedPhyNetwork(
                undirected_edges=[
                    {"u": 3, "v": 1, "bootstrap": 1.5},
                    (3, 2),
                    (3, 4),
                ],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})],
            )

    def test_bootstrap_nan(self) -> None:
        """Bootstrap value as NaN should be invalid."""
        with pytest.raises(ValueError, match="Bootstrap value .* is nan, but must be in \\[0\\.0, 1\\.0\\]"):
            SemiDirectedPhyNetwork(
                undirected_edges=[
                    {"u": 3, "v": 1, "bootstrap": math.nan},
                    (3, 2),
                    (3, 4),
                ],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})],
            )


class TestInvalidGamma:
    """Test that invalid gamma values are caught."""

    def test_gamma_partial(self) -> None:
        """Partial gamma values (not all hybrid edges have gamma) should be invalid."""
        with pytest.raises(ValueError, match="gamma|Gamma"):
            SemiDirectedPhyNetwork(
                directed_edges=[
                    {"u": 5, "v": 4, "gamma": 0.6},
                    {"u": 6, "v": 4},  # Missing gamma
                ],
                undirected_edges=[
                    (7, 5),
                    (7, 6),
                    (7, 10),
                    (5, 8),
                    (6, 9),
                    (4, 1),
                ],
                nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'})],
            )

    def test_gamma_not_summing_to_one(self) -> None:
        """Gamma values not summing to 1 should be invalid."""
        with pytest.raises(ValueError, match="gamma|Gamma"):
            SemiDirectedPhyNetwork(
                directed_edges=[
                    {"u": 5, "v": 4, "gamma": 0.6},
                    {"u": 6, "v": 4, "gamma": 0.5},  # Sum = 1.1 (invalid)
                ],
                undirected_edges=[
                    (7, 5),
                    (7, 6),
                    (7, 10),
                    (5, 8),
                    (6, 9),
                    (4, 1),
                ],
                nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'})],
            )

    def test_gamma_below_zero(self) -> None:
        """Gamma value below 0 should be invalid."""
        with pytest.raises(ValueError, match="Gamma value .* is .*, but must be in \\[0\\.0, 1\\.0\\]"):
            SemiDirectedPhyNetwork(
                directed_edges=[
                    {"u": 5, "v": 4, "gamma": -0.1},
                    {"u": 6, "v": 4, "gamma": 1.1},
                ],
                undirected_edges=[
                    (7, 5),
                    (7, 6),
                    (7, 10),
                    (5, 8),
                    (6, 9),
                    (4, 1),
                ],
                nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'})],
            )

    def test_gamma_above_one(self) -> None:
        """Gamma value above 1 should be invalid."""
        with pytest.raises(ValueError, match="Gamma value .* is .*, but must be in \\[0\\.0, 1\\.0\\]"):
            SemiDirectedPhyNetwork(
                directed_edges=[
                    {"u": 5, "v": 4, "gamma": 0.5},
                    {"u": 6, "v": 4, "gamma": 1.5},  # > 1.0
                ],
                undirected_edges=[
                    (7, 5),
                    (7, 6),
                    (7, 10),
                    (5, 8),
                    (6, 9),
                    (4, 1),
                ],
                nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'})],
            )


class TestEdgeCases:
    """Test edge cases for validation."""

    def test_empty_network(self) -> None:
        """Empty network should be valid (special case)."""
        with pytest.warns(UserWarning, match="Empty network.*no nodes"):
            net = SemiDirectedPhyNetwork(
                directed_edges=[],
                undirected_edges=[], nodes=[],
            )
        # Empty networks skip full validation
        with pytest.warns(UserWarning, match="Empty network.*no nodes"):
            assert net.validate() is True

    def test_single_node_network(self) -> None:
        """Single-node semi-directed network is valid but warns."""
        with pytest.warns(UserWarning, match="Single-node network detected"):
            net = SemiDirectedPhyNetwork(
                directed_edges=[],
                undirected_edges=[],
                nodes=[(1, {"label": "A"})],
            )
        assert net.number_of_nodes() == 1
        assert net.number_of_edges() == 0
        assert net.leaves == {1}
        with pytest.warns(UserWarning, match="Single-node network detected"):
            assert net.validate() is True

    def test_two_nodes_undirected_edge(self) -> None:
        """Two nodes with undirected edge should be valid."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(1, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})],
        )
        assert net.validate() is True

    def test_two_nodes_directed_edge_invalid(self) -> None:
        """Two nodes with directed edge should be invalid (doesn't meet 2-node special case)."""
        with pytest.raises(ValueError):
            SemiDirectedPhyNetwork(
                directed_edges=[(1, 2)],
                undirected_edges=[],
                nodes=[(2, {'label': 'A'})],
            )

    def test_two_nodes_both_edges_invalid(self) -> None:
        """Two nodes with both directed and undirected edges should be invalid."""
        with pytest.raises(ValueError):
            SemiDirectedPhyNetwork(
                directed_edges=[(1, 2)],
                undirected_edges=[(1, 2)],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})],
            )


class TestComplexValidation:
    """Test complex validation scenarios."""

    def test_valid_with_all_attributes(self) -> None:
        """Network with branch_length, bootstrap, and gamma should validate."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {"u": 5, "v": 4, "gamma": 0.6, "branch_length": 0.1, "bootstrap": 0.9},
                {"u": 6, "v": 4, "gamma": 0.4, "branch_length": 0.2, "bootstrap": 0.85},
            ],
            undirected_edges=[
                {"u": 7, "v": 5, "branch_length": 0.3, "bootstrap": 0.95},
                {"u": 7, "v": 6, "branch_length": 0.4, "bootstrap": 0.92},
                {"u": 7, "v": 10, "branch_length": 0.45},
                {"u": 5, "v": 8, "branch_length": 0.5},
                {"u": 6, "v": 9, "branch_length": 0.6},
                {"u": 4, "v": 1, "branch_length": 0.7},
            ],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'})],
        )
        assert net.validate() is True

    def test_parallel_hybrid_edges_valid(self) -> None:
        """Network with parallel hybrid edges should validate."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {"u": 5, "v": 4, "key": 0, "gamma": 0.3},
                {"u": 5, "v": 4, "key": 1, "gamma": 0.2},
                {"u": 6, "v": 4, "gamma": 0.5},
            ],
            undirected_edges=[
                (7, 5),
                (7, 6),
                (7, 10),
                (5, 8),
                (6, 9),
                (4, 1),
            ],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'})],
        )
        assert net.validate() is True

    def test_deep_tree_valid(self) -> None:
        """Deep tree structure should validate."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[
                (10, 8),
                (10, 9),
                (10, 16),
                (8, 6),
                (8, 7),
                (6, 4),
                (6, 5),
                (4, 1),
                (4, 2),
                (5, 3),
                (5, 11),
                (7, 12),
                (7, 13),
                (9, 14),
                (9, 15),
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (11, {'label': 'D'}), (12, {'label': 'E'}), (13, {'label': 'F'}), (14, {'label': 'G'}), (15, {'label': 'H'}), (16, {'label': 'I'})],
        )
        assert net.validate() is True


class TestValidationMessages:
    """Test that validation error messages are informative."""

    def test_disconnected_error_message(self) -> None:
        """Disconnected network should have informative error."""
        with pytest.raises(ValueError) as exc_info:
            SemiDirectedPhyNetwork(
                undirected_edges=[
                    (3, 1),
                    (3, 2),
                    (3, 7),
                    (6, 4),
                    (6, 5),
                    (6, 8),
                ],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (5, {'label': 'D'}), (7, {'label': 'E'}), (8, {'label': 'F'})],
            )
        assert "not connected" in str(exc_info.value).lower()

    def test_degree_constraint_error_message(self) -> None:
        """Invalid degree should have informative error."""
        with pytest.raises(ValueError) as exc_info:
            SemiDirectedPhyNetwork(
                undirected_edges=[(3, 1), (3, 2)],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})],
            )
        assert "degree" in str(exc_info.value).lower()

    def test_gamma_error_message(self) -> None:
        """Invalid gamma should have informative error."""
        with pytest.raises(ValueError) as exc_info:
            SemiDirectedPhyNetwork(
                directed_edges=[
                    {"u": 5, "v": 4, "gamma": 0.6},
                    {"u": 6, "v": 4, "gamma": 0.5},  # Sum != 1
                ],
                undirected_edges=[
                    (7, 5),
                    (7, 6),
                    (7, 10),
                    (5, 8),
                    (6, 9),
                    (4, 1),
                ],
                nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'})],
            )
        assert "gamma" in str(exc_info.value).lower()

    def test_bootstrap_error_message(self) -> None:
        """Invalid bootstrap should have informative error."""
        with pytest.raises(ValueError) as exc_info:
            SemiDirectedPhyNetwork(
                undirected_edges=[
                    {"u": 3, "v": 1, "bootstrap": 1.5},
                    (3, 2),
                    (3, 4),
                ],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})],
            )
        assert "bootstrap" in str(exc_info.value).lower()


class TestValidationPerformance:
    """Test validation performance on larger networks."""

    def test_large_tree_validates_quickly(self) -> None:
        """Large tree should validate in reasonable time."""
        import time

        # Build a large ternary tree
        edges = []
        node_id = 1
        current_level = [0]
        nodes = []

        # Build 4 levels (81 leaves)
        for level in range(4):
            next_level = []
            for parent in current_level:
                for _ in range(3):
                    child = node_id
                    node_id += 1
                    next_level.append(child)
                    edges.append((parent, child))
            current_level = next_level

        # Mark last level as leaves
        for i, leaf in enumerate(current_level):
            nodes.append((leaf, {'label': f"Taxon{i}"}))

        start = time.time()
        net = SemiDirectedPhyNetwork(
            undirected_edges=edges,
            nodes=nodes,
        )
        elapsed = time.time() - start

        # Validation should complete in less than 1 second
        assert elapsed < 1.0
        assert net.validate() is True


class TestInvalidSemiDirectedStructures:
    """Test structures that violate semi-directed network constraints."""

    def test_multiple_source_components_disconnected(self) -> None:
        """Network with multiple source components (disconnected) should raise error."""
        # This is caught by connectivity check, but tests the error message
        with pytest.raises(ValueError, match="not connected"):
            SemiDirectedPhyNetwork(
                directed_edges=[
                    (5, 4, 0),
                    (6, 4, 1),
                    (8, 7, 0),  # Separate disconnected component
                    (9, 7, 1),
                ],
                undirected_edges=[
                    # Component 1
                    (10, 5),
                    (10, 6),
                    (10, 12),
                    (5, 14),
                    (6, 15),
                    (4, 1),
                    # Component 2 (disconnected)
                    (11, 8),
                    (11, 9),
                    (11, 13),
                    (8, 16),
                    (9, 17),
                    (7, 2),
                ],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (12, {'label': 'C'}), (13, {'label': 'D'}), (14, {'label': 'E'}), (15, {'label': 'F'}), (16, {'label': 'G'}), (17, {'label': 'H'})],
            )

    def test_source_component_only_leaves(self) -> None:
        """Source component with only leaves should be invalid."""
        # Create network where source component nodes are all leaves
        # Nodes 4, 5, 6 have no incoming directed edges, are connected via undirected edges,
        # but are all leaves (outdegree 0), so no internal node to pick as root
        # Note: This is caught earlier because nodes in source component can't be in taxa
        # if they have outgoing edges. Let's test a case where this actually happens.
        with pytest.raises(ValueError):
            # This test verifies that networks with source components containing only leaves
            # are caught. The exact error may vary based on which check catches it first.
            SemiDirectedPhyNetwork(
                directed_edges=[
                    (7, 3),
                    (8, 3),
                ],
                undirected_edges=[
                    (4, 5),  # Source component: nodes 4, 5, 6
                    (5, 6),
                    (4, 6),
                    (3, 1),
                    (7, 9),
                    (7, 10),
                    (8, 11),
                    (8, 12),
                ],
                nodes=[(1, {'label': 'A'}), (4, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'}), (9, {'label': 'E'}), (10, {'label': 'F'}), (11, {'label': 'G'}), (12, {'label': 'H'})],
            )

    def test_invalid_cross_pattern_network_creates_cycle(self) -> None:
        """
        Test invalid network with cross-pattern directed edges (from image left graph).
        
        This network has a cross-pattern of directed edges that creates a directed cycle
        when oriented, making it invalid.
        """
        with pytest.raises(ValueError):
            SemiDirectedPhyNetwork(
                directed_edges=[
                    (5, 6),  # Creates cycle: 5->6->7->8->5
                    (7, 8),
                    (6, 7),
                    (8, 5),
                ],
                undirected_edges=[
                    (19, 5),  # Source node 19
                    (19, 7),
                    (19, 20),
                    (33, 6),  # Leaf v1 (node 33)
                    (34, 8),  # Leaf v2 (node 34)
                    (3, 33),  # Connect internal node to leaf
                    (3, 13),
                    (4, 34),  # Connect internal node to leaf
                    (4, 14),
                    (9, 3),
                    (10, 4),
                    (5, 15),
                    (6, 16),
                    (7, 17),
                    (8, 18),
                ],
                nodes=[(13, {'label': 'A'}), (14, {'label': 'B'}), (15, {'label': 'C'}), (16, {'label': 'D'}), (17, {'label': 'E'}), (18, {'label': 'F'}), (20, {'label': 'G'}), (33, {'label': 'v1'}), (34, {'label': 'v2'})],
            )

    def test_invalid_dag_structure(self) -> None:
        """
        Test invalid network with DAG-like structure.
        
        This network structure may not roundtrip correctly due to the specific
        arrangement of directed and undirected edges.
        """
        with pytest.raises(ValueError, match="indegree"):
            SemiDirectedPhyNetwork(
                directed_edges=[
                    (11, 12),  # v1 -> v2 (node 12 gets indegree 1)
                    (11, 15),  # v1 -> v5
                ],
                undirected_edges=[
                    # Top section
                    (1, 2),
                    (1, 18),
                    (1, 26),
                    (2, 3),
                    (2, 4),
                    (2, 5),
                    (3, 6),
                    (3, 19),
                    (4, 7),
                    (4, 20),
                    (5, 8),
                    (5, 21),
                    (8, 9),
                    (8, 10),
                    (8, 11),  # Connect to v1
                    (9, 13),
                    (9, 22),
                    (10, 14),
                    (10, 23),
                    # Bottom section
                    (12, 13),  # v2 -> v3 (undirected) - node 12 has total_degree 3
                    (12, 24),  # But indegree 1, which violates constraint
                    (13, 16),
                    (15, 14),
                    (15, 25),
                    (14, 17),
                    (13, 14),
                ],
                nodes=[(6, {'label': 'A'}), (7, {'label': 'B'}), (16, {'label': 'C'}), (17, {'label': 'D'}), (18, {'label': 'E'}), (19, {'label': 'F'}), (20, {'label': 'G'}), (21, {'label': 'H'}), (22, {'label': 'I'}), (23, {'label': 'J'}), (24, {'label': 'K'}), (25, {'label': 'L'}), (26, {'label': 'M'})],
            )

    def test_network_with_undirected_cycle_creating_directed_cycle(self) -> None:
        """Network where undirected cycle creates directed cycle when oriented should be invalid."""
        # Create a network where an undirected cycle becomes a directed cycle when oriented
        # This should fail validation (may be caught at different stages)
        # Structure: undirected edges 4-5-6-4 form a cycle that, when oriented, creates directed cycle
        with pytest.raises(ValueError):
            SemiDirectedPhyNetwork(
                directed_edges=[
                    (15, 5),  # Hybrid edge into node 5
                    (16, 6),  # Hybrid edge into node 6
                ],
                undirected_edges=[
                    (12, 1),  # Source node 12
                    (12, 2),
                    (12, 13),
                    (1, 19),  # Leaf
                    (1, 20),  # Leaf
                    (2, 3),
                    (2, 21),
                    (3, 4),
                    (3, 22),
                    (4, 5),  # Undirected cycle: 4-5-6-4
                    (4, 6),  # Completes the cycle
                    (4, 23),
                    (5, 11),
                    (5, 14),  # Ensure node 5 has degree >= 3
                    (6, 24),  # Leaf
                    (6, 25),  # Leaf
                    (15, 17),  # Ensure node 15 has degree >= 3
                    (15, 26),  # Extra edge for node 15
                    (16, 18),  # Ensure node 16 has degree >= 3
                    (16, 27),  # Extra edge for node 16
                ],
                nodes=[(11, {'label': 'A'}), (13, {'label': 'B'}), (14, {'label': 'C'}), (17, {'label': 'D'}), (18, {'label': 'E'}), (19, {'label': 'F'}), (20, {'label': 'G'}), (21, {'label': 'H'}), (22, {'label': 'I'}), (23, {'label': 'J'}), (24, {'label': 'K'}), (25, {'label': 'L'}), (26, {'label': 'M'}), (27, {'label': 'N'})],
            )
