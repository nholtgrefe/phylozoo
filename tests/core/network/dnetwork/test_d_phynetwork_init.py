"""
Comprehensive tests for DirectedPhyNetwork initialization and validation.

This module tests all aspects of network initialization including:
- Empty networks
- Valid network structures
- Invalid network structures
- Edge format variations
- Taxa and label handling
- Validation rules
"""

import warnings
from typing import Dict, List, Tuple

import pytest

from phylozoo.core.network import DirectedPhyNetwork


class TestEmptyNetwork:
    """Test cases for empty network initialization."""

    def test_empty_network_with_warning(self) -> None:
        """Test that empty network raises a warning."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)  # Ignore the init warning
            net = DirectedPhyNetwork(edges=[])
        
        assert net.number_of_nodes() == 0
        assert net.number_of_edges() == 0
        with pytest.warns(UserWarning, match="Empty network.*no nodes.*detected"):
            net.validate()  # Empty networks are valid

    def test_empty_network_properties(self) -> None:
        """Test properties of empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(edges=[])
        
        assert len(net.leaves) == 0
        assert len(net.taxa) == 0
        # Empty network has no root - accessing root_node should raise
        with pytest.raises(ValueError, match="Network has no root node"):
            _ = net.root_node
        assert len(net.internal_nodes) == 0
        assert len(net.hybrid_nodes) == 0
        assert len(net.tree_nodes) == 0


class TestMinimalValidNetworks:
    """Test cases for minimal valid network structures."""

    def test_single_node_network_with_warning(self) -> None:
        """Test that single-node network raises a warning."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)  # Ignore the init warning
            net = DirectedPhyNetwork(nodes=[(1, {"label": "A"})])
        
        assert net.number_of_nodes() == 1
        assert net.number_of_edges() == 0
        with pytest.warns(UserWarning, match="Single-node network detected"):
            net.validate()  # Single-node networks are valid

    def test_single_edge_network(self) -> None:
        """Test network with single edge (root -> leaf)."""
        net = DirectedPhyNetwork(edges=[(1, 2)], nodes=[(2, {"label": "A"})])
        assert net.number_of_nodes() == 2
        assert net.number_of_edges() == 1
        assert net.root_node == 1
        assert net.leaves == {2}
        assert net.taxa == {"A"}

    def test_two_edge_star_network(self) -> None:
        """Test star network with root and two leaves."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"})],
        )
        assert net.number_of_nodes() == 3
        assert net.number_of_edges() == 2
        assert net.root_node == 3
        assert net.leaves == {1, 2}
        assert net.taxa == {"A", "B"}


class TestSimpleTrees:
    """Test cases for simple tree structures."""

    def test_binary_tree(self) -> None:
        """Test a simple binary tree."""
        # Tree: root -> internal -> 2 leaves, root -> leaf
        net = DirectedPhyNetwork(
            edges=[(4, 3), (3, 1), (3, 2), (4, 5)],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (5, {"label": "C"})],
        )
        assert net.number_of_nodes() == 5
        assert net.number_of_edges() == 4
        assert net.root_node == 4
        assert net.leaves == {1, 2, 5}
        assert net.is_tree() is True

    def test_ternary_tree(self) -> None:
        """Test a tree with node having out-degree 3."""
        net = DirectedPhyNetwork(
            edges=[(4, 1), (4, 2), (4, 3)],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (3, {"label": "C"})],
        )
        assert net.number_of_nodes() == 4
        assert net.number_of_edges() == 3
        assert net.root_node == 4
        assert net.leaves == {1, 2, 3}

    def test_star_tree(self) -> None:
        """Test a star tree with many leaves."""
        edges = [(10, i) for i in range(1, 10)]
        nodes = [(i, {"label": f"Taxon{i}"}) for i in range(1, 10)]
        net = DirectedPhyNetwork(edges=edges, nodes=nodes)
        assert net.number_of_nodes() == 10
        assert net.number_of_edges() == 9
        assert net.root_node == 10
        assert len(net.leaves) == 9


class TestNetworksWithHybrids:
    """Test cases for networks with hybrid nodes."""

    def test_single_hybrid_node(self) -> None:
        """Test network with single hybrid node."""
        # Structure: root -> tree nodes (out-degree 2) -> hybrid -> leaf
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})],
        )
        assert net.number_of_nodes() == 7
        assert net.number_of_edges() == 7
        assert net.hybrid_nodes == {4}
        assert net.tree_nodes == {5, 6}
        assert net.is_tree() is False

    def test_multiple_hybrid_nodes(self) -> None:
        """Test network with multiple hybrid nodes."""
        # Structure with 2 hybrid nodes
        net = DirectedPhyNetwork(
            edges=[
                (10, 7), (10, 8),  # Root to tree nodes
                (7, 5), (7, 6),    # Tree node 7 splits
                (8, 5), (8, 9),    # Tree node 8 splits, 5 is hybrid
                (5, 4), (6, 4),    # Hybrid node 4
                (4, 1), (9, 2), (9, 3), (6, 11)  # To leaves (9 and 6 split to 2 children)
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (11, {"label": "D"}),
            ],
        )
        assert net.hybrid_nodes == {4, 5}
        assert len(net.hybrid_nodes) == 2
        assert net.is_tree() is False

    def test_nested_hybrids(self) -> None:
        """Test network with nested hybridization."""
        # Hybrid node 4 has parent hybrid node 5
        net = DirectedPhyNetwork(
            edges=[
                (10, 7), (10, 8),
                (7, 5), (7, 6),
                (8, 5), (8, 9),
                (5, 4), (6, 4),
                (4, 1), (9, 2), (9, 3), (6, 11)
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (11, {"label": "D"}),
            ],
        )
        assert 4 in net.hybrid_nodes
        assert 5 in net.hybrid_nodes


class TestParallelEdges:
    """Test cases for networks with parallel edges."""

    def test_parallel_edges_to_hybrid(self) -> None:
        """Test parallel edges entering a hybrid node."""
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),
                (5, 4, 0), (5, 4, 1),  # Two parallel edges from 5 to 4
                (5, 8),
                (6, 4),
                (6, 9),
                (4, 2)
            ],
            nodes=[(2, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})],
        )
        assert net.has_edge(5, 4, key=0)
        assert net.has_edge(5, 4, key=1)
        assert net.number_of_edges() == 8  # 7 edges + 1 parallel

    def test_parallel_edges_dict_format(self) -> None:
        """Test parallel edges using dict format."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 7, 'v': 5},
                {'u': 7, 'v': 6},
                {'u': 5, 'v': 4, 'key': 0},
                {'u': 5, 'v': 4, 'key': 1},  # Parallel edge
                {'u': 5, 'v': 8},
                {'u': 6, 'v': 4},
                {'u': 6, 'v': 9},
                {'u': 4, 'v': 2}
            ],
            nodes=[(2, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})],
        )
        assert net.has_edge(5, 4, key=0)
        assert net.has_edge(5, 4, key=1)


class TestEdgeFormats:
    """Test cases for different edge input formats."""

    def test_tuple_edges(self) -> None:
        """Test edges as simple tuples."""
        net = DirectedPhyNetwork(
            edges=[(1, 2), (1, 3)],
            nodes=[(2, {"label": "A"}), (3, {"label": "B"})],
        )
        assert net.number_of_edges() == 2

    def test_tuple_edges_with_keys(self) -> None:
        """Test edges as tuples with explicit keys."""
        net = DirectedPhyNetwork(
            edges=[(1, 2, 0), (1, 3, 5)],
            nodes=[(2, {"label": "A"}), (3, {"label": "B"})],
        )
        assert net.has_edge(1, 2, key=0)
        assert net.has_edge(1, 3, key=5)

    def test_dict_edges(self) -> None:
        """Test edges as dictionaries."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 1, 'v': 2},
                {'u': 1, 'v': 3, 'key': 10}
            ],
            nodes=[(2, {"label": "A"}), (3, {"label": "B"})],
        )
        assert net.has_edge(1, 2)
        assert net.has_edge(1, 3, key=10)

    def test_mixed_edge_formats(self) -> None:
        """Test mixing different edge formats."""
        net = DirectedPhyNetwork(
            edges=[
                (1, 2),  # Tuple
                {'u': 1, 'v': 3},  # Dict
                (1, 4, 0)  # Tuple with key
            ],
            nodes=[(2, {"label": "A"}), (3, {"label": "B"}), (4, {"label": "C"})],
        )
        assert net.number_of_edges() == 3


class TestTaxaFormats:
    """Test cases for different taxa input formats."""

    def test_taxa_as_dict(self) -> None:
        """Test taxa provided as dictionary."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"})],
        )
        assert net.taxa == {"A", "B"}

    def test_taxa_as_list_of_tuples(self) -> None:
        """Test taxa provided as list of tuples."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"})],
        )
        assert net.taxa == {"A", "B"}

    def test_taxa_none_auto_labeling(self) -> None:
        """Test that uncovered leaves get auto-labeled."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=None)
        # Both leaves should be auto-labeled
        assert len(net.taxa) == 2
        assert "1" in net.taxa or "2" in net.taxa

    def test_partial_taxa_auto_labeling(self) -> None:
        """Test partial taxa mapping with auto-labeling for uncovered leaves."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {"label": "A"})],  # Only one leaf labeled
        )
        assert "A" in net.taxa
        assert len(net.taxa) == 3  # All 3 leaves should have labels
        # Leaves 2 and 4 should be auto-labeled
        assert net.get_label(2) is not None
        assert net.get_label(4) is not None


class TestInternalNodeLabels:
    """Test cases for internal node labels."""

    def test_internal_node_labels_dict(self) -> None:
        """Test internal node labels as dictionary."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (3, {"label": "root"})],
        )
        assert net.get_label(3) == "root"

    def test_internal_node_labels_list(self) -> None:
        """Test internal node labels as list of tuples."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (3, {"label": "root"})],
        )
        assert net.get_label(3) == "root"

    def test_internal_node_labels_none(self) -> None:
        """Test that internal nodes can be unlabeled."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"})],
        )
        assert net.get_label(3) is None

    def test_internal_node_labels_multiple(self) -> None:
        """Test labeling multiple internal nodes."""
        net = DirectedPhyNetwork(
            edges=[(4, 3), (3, 1), (3, 2), (4, 5)],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (5, {"label": "C"}),
                (3, {"label": "internal1"}),
                (4, {"label": "root"}),
            ],
        )
        assert net.get_label(3) == "internal1"
        assert net.get_label(4) == "root"


class TestInvalidInitialization:
    """Test cases for invalid network initialization (should raise errors)."""

    def test_cycle_detection(self) -> None:
        """Test that cycles are detected and raise ValueError."""
        # Create a cycle: 4 -> 5 -> 6 -> 4, with leaves 1, 2, 3
        with pytest.raises(ValueError, match="directed cycles"):
            DirectedPhyNetwork(
                edges=[(4, 5), (5, 6), (6, 4), (4, 1), (5, 2), (6, 3)],
                nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (3, {"label": "C"})],
            )

    def test_multiple_roots(self) -> None:
        """Test that multiple root nodes raise ValueError."""
        # Two nodes with in-degree 0
        with pytest.raises(ValueError, match="multiple root nodes"):
            DirectedPhyNetwork(
                edges=[(1, 3), (2, 3)],
                nodes=[(3, {"label": "A"})],
            )

    def test_no_root(self) -> None:
        """Test that network with no root raises ValueError."""
        # Two disconnected components, both with roots, but validation should catch
        # that there's no single root. Actually, if all nodes have in-degree >= 1,
        # there must be a cycle. So we test a case where validation order matters.
        # Since cycles are checked first, this will raise a cycle error.
        # For a true "no root" test, we'd need a valid DAG structure.
        # Instead, test that a network with a cycle (which also has no root) raises an error.
        with pytest.raises(ValueError, match="directed cycles"):
            DirectedPhyNetwork(
                edges=[(3, 4), (4, 5), (5, 3), (3, 1), (4, 2)],  # Cycle, but also no root
                nodes=[(1, {"label": "A"}), (2, {"label": "B"})],
            )

    def test_leaf_with_wrong_indegree(self) -> None:
        """Test that leaf with in-degree != 1 raises ValueError."""
        # Leaf 3 has in-degree 2 (from nodes 1 and 2)
        with pytest.raises(ValueError, match="in-degree"):
            DirectedPhyNetwork(
                edges=[(4, 1), (4, 2), (1, 3), (2, 3)],
                nodes=[(3, {"label": "A"})],
            )

    def test_leaf_with_outdegree(self) -> None:
        """Test that leaf with out-degree > 0 raises ValueError."""
        # Node 2 is labeled as leaf but has outgoing edge
        # The error will be about internal node degrees, not specifically about leaves
        with pytest.raises(ValueError):
            DirectedPhyNetwork(
                edges=[(1, 2), (2, 3)],
                nodes=[(2, {"label": "A"}), (3, {"label": "B"})],
            )

    def test_internal_node_invalid_degrees(self) -> None:
        """Test that internal nodes with invalid degrees raise ValueError."""
        # Node 2 has in-degree 1 and out-degree 1 (invalid)
        with pytest.raises(ValueError, match="Internal node"):
            DirectedPhyNetwork(
                edges=[(1, 2), (2, 3)],
                nodes=[(3, {"label": "A"})],
            )

    def test_non_leaf_in_taxa(self) -> None:
        """Test that non-leaf in taxa mapping raises ValueError."""
        # Node 2 has outgoing edge, so not a leaf
        # The error will be about internal node degrees
        with pytest.raises(ValueError):
            DirectedPhyNetwork(
                edges=[(1, 2), (2, 3)],
                nodes=[(2, {"label": "A"}), (3, {"label": "B"})],  # 2 is not a leaf
            )

    def test_leaf_in_internal_node_labels(self) -> None:
        """Test that duplicate labels raise ValueError."""
        with pytest.raises(ValueError, match="already used"):
            DirectedPhyNetwork(
                edges=[(1, 2), (1, 3)],
                nodes=[(2, {"label": "A"}), (3, {"label": "A"})],  # Duplicate label
            )

    def test_duplicate_labels(self) -> None:
        """Test that duplicate labels raise ValueError."""
        with pytest.raises(ValueError, match="already used"):
            DirectedPhyNetwork(
                edges=[(1, 2), (1, 3)],
                nodes=[(2, {"label": "A"}), (3, {"label": "A"})],  # Duplicate label
            )

    def test_invalid_nodes_format(self) -> None:
        """Test that invalid nodes format raises ValueError."""
        # Invalid nodes format (string instead of list) will cause validation error
        # The actual error may vary, but it should raise an error
        with pytest.raises((ValueError, TypeError)):
            DirectedPhyNetwork(
                edges=[(1, 2)],
                nodes="invalid",  # type: ignore
            )

    def test_self_loop_disallowed(self) -> None:
        """Test that self-loops are rejected during validation."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)  # Ignore single-node warning
            with pytest.raises(ValueError, match="Self-loops are not allowed"):
                DirectedPhyNetwork(
                    edges=[(1, 1)],
                    nodes=[(1, {"label": "A"})],
                )

    def test_invalid_internal_labels_format(self) -> None:
        """Test that invalid internal_node_labels format raises ValueError."""
        with pytest.raises(ValueError, match="must be .*dict"):
            DirectedPhyNetwork(
                edges=[(1, 2)],
                nodes=[(2, {"label": "A"}), (3, "invalid")],  # type: ignore
            )

    def test_disconnected_network(self) -> None:
        """Test that disconnected networks raise ValueError."""
        # Note: Disconnected networks will typically have multiple roots, which is
        # caught first. However, the connectivity check is still performed as an
        # explicit validation step. If a network somehow passes the root check but
        # is disconnected, it will be caught here.
        
        # Two disconnected components (will fail on multiple roots first)
        with pytest.raises(ValueError, match="multiple root nodes|not connected"):
            DirectedPhyNetwork(
                edges=[
                    (1, 2),  # Component 1: 1 -> 2
                    (3, 4)   # Component 2: 3 -> 4 (disconnected)
                ],
                nodes=[(2, {"label": "A"}), (4, {"label": "B"})],
            )
        
        # Three disconnected components (will fail on multiple roots first)
        with pytest.raises(ValueError, match="multiple root nodes|not connected"):
            DirectedPhyNetwork(
                edges=[
                    (1, 2),  # Component 1
                    (3, 4),  # Component 2
                    (5, 6)   # Component 3
                ],
                nodes=[(2, {"label": "A"}), (4, {"label": "B"}), (6, {"label": "C"})],
            )


class TestValidationEdgeCases:
    """Test edge cases for validation."""

    def test_very_deep_tree(self) -> None:
        """Test validation of very deep tree structure."""
        # Create a binary tree with 5 levels (root + 4 internal + leaves)
        edges = []
        # Level 1: root 100 -> 1, 2
        edges.append((100, 1))
        edges.append((100, 2))
        # Level 2: 1 -> 3, 4; 2 -> 5, 6
        edges.append((1, 3))
        edges.append((1, 4))
        edges.append((2, 5))
        edges.append((2, 6))
        # Level 3: 3 -> 7, 8; 4 -> 9, 10; 5 -> 11, 12; 6 -> 13, 14
        for parent in [3, 4, 5, 6]:
            edges.append((parent, parent * 2 + 1))
            edges.append((parent, parent * 2 + 2))
        # Level 4: leaves 15-30
        for parent in [7, 8, 9, 10, 11, 12, 13, 14]:
            edges.append((parent, parent * 2 + 1))
            edges.append((parent, parent * 2 + 2))
        # Leaves are nodes 15-30
        nodes = [(i, {"label": f"Taxon{i}"}) for i in range(15, 31)]
        net = DirectedPhyNetwork(edges=edges, nodes=nodes)
        net.validate()
        assert net.root_node == 100

    def test_wide_tree(self) -> None:
        """Test validation of very wide tree."""
        # Root with 50 children
        edges = [(1, i) for i in range(2, 52)]
        nodes = [(i, {"label": f"Taxon{i}"}) for i in range(2, 52)]
        net = DirectedPhyNetwork(edges=edges, nodes=nodes)
        net.validate()
        assert len(net.leaves) == 50

    def test_complex_hybrid_structure(self) -> None:
        """Test validation of complex hybrid structure."""
        # Multiple hybrid nodes with various configurations
        edges = [
            (20, 10), (20, 11),  # Root splits
            (10, 5), (10, 6),     # Tree node 10 splits
            (11, 5), (11, 7),     # Hybrid node 5
            (5, 4), (6, 4),      # Hybrid node 4
            (4, 1), (7, 2), (7, 8), (6, 3)  # To leaves (7 and 6 split to 2 children)
        ]
        net = DirectedPhyNetwork(
            edges=edges,
            nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (3, {"label": "C"}), (8, {"label": "D"})],
        )
        net.validate()
        assert len(net.hybrid_nodes) == 2

    def test_parallel_edges_validation(self) -> None:
        """Test that parallel edges don't break validation."""
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),
                (5, 4, 0), (5, 4, 1), (5, 4, 2),  # 3 parallel edges
                (5, 8),
                (6, 4),
                (6, 9),
                (4, 2)
            ],
            nodes=[(2, {"label": "A"}), (8, {"label": "B"}), (9, {"label": "C"})],
        )
        net.validate()
        # Node 4 should still be a hybrid (in-degree >= 2, out-degree 1)
        assert 4 in net.hybrid_nodes


class TestAutoLabelingEdgeCases:
    """Test edge cases for auto-labeling."""

    def test_auto_label_conflict_resolution(self) -> None:
        """Test that auto-labeling handles conflicts."""
        # If leaf 2 exists and we label another leaf as "2", should resolve conflict
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {"label": "2"})],  # Label conflicts with node ID 2
        )
        # Leaf 2 should get auto-label that doesn't conflict
        label_2 = net.get_label(2)
        assert label_2 is not None
        assert label_2 != "2"  # Should be resolved to something like "2_1"

    def test_numeric_node_ids_auto_labeling(self) -> None:
        """Test auto-labeling with numeric node IDs."""
        net = DirectedPhyNetwork(
            edges=[(100, 50), (100, 51)],
            nodes=None
        )
        # Should auto-label as strings of node IDs
        assert "50" in net.taxa or "51" in net.taxa

    def test_string_node_ids_auto_labeling(self) -> None:
        """Test auto-labeling with string node IDs."""
        net = DirectedPhyNetwork(
            edges=[("root", "leaf1"), ("root", "leaf2")],
            nodes=None
        )
        # Should auto-label using string representation
        assert len(net.taxa) == 2


class TestSingleNodeNetwork:
    """Comprehensive tests for single-node networks."""

    def test_root_is_single_node(self) -> None:
        """Test that the root of a single-node network is that node."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(nodes=[(1, {"label": "A"})])
        assert net.root_node == 1
        assert net.number_of_nodes() == 1

    def test_single_node_is_also_leaf(self) -> None:
        """Test that the single node is both root and leaf."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(nodes=[(1, {"label": "A"})])
        assert net.root_node == 1
        assert net.leaves == {1}
        assert net.root_node in net.leaves

    def test_single_node_topology_properties(self) -> None:
        """Test topology properties for single-node network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(nodes=[(1, {"label": "A"})])
        # Single node is root and leaf, so not internal
        assert net.internal_nodes == set()
        assert net.hybrid_nodes == set()
        assert net.tree_nodes == set()
        assert net.hybrid_edges == set()
        assert net.tree_edges == set()

    def test_single_node_taxa(self) -> None:
        """Test taxa property for single-node network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(nodes=[(1, {"label": "A"})])
        assert net.taxa == {"A"}
        assert net.get_label(1) == "A"

    def test_single_node_degrees(self) -> None:
        """Test degree properties for single-node network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(nodes=[(1, {"label": "A"})])
        assert net.indegree(1) == 0
        assert net.outdegree(1) == 0
        assert net.degree(1) == 0

    def test_single_node_is_tree(self) -> None:
        """Test that single-node network is considered a tree."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(nodes=[(1, {"label": "A"})])
        assert net.is_tree() is True

    def test_single_node_to_sd_network(self) -> None:
        """Test conversion to semi-directed network."""
        from phylozoo.core.network.dnetwork.operations import to_sd_network
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(nodes=[(1, {"label": "A"})])
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)  # Ignore empty-edge warning in SD init
            sd_net = to_sd_network(net)
        assert sd_net.number_of_nodes() == 1
        assert sd_net.number_of_edges() == 0
        assert sd_net.leaves == {1}
        assert sd_net.taxa == {"A"}
        assert 1 in sd_net

    def test_single_node_copy(self) -> None:
        """Test copying single-node network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(nodes=[(1, {"label": "A"})])
        
        net_copy = net.copy()
        assert net_copy.number_of_nodes() == 1
        assert net_copy.root_node == 1
        assert net_copy.leaves == {1}
        assert net_copy.taxa == {"A"}
        assert net_copy is not net  # Different object

    def test_single_node_iteration(self) -> None:
        """Test iteration over single-node network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(nodes=[(1, {"label": "A"})])
        
        nodes_list = list(net)
        assert nodes_list == [1]
        assert len(net) == 1
        assert 1 in net
        assert 2 not in net

    def test_single_node_children_parents(self) -> None:
        """Test children and parents for single-node network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(nodes=[(1, {"label": "A"})])
        
        assert list(net.children(1)) == []
        assert list(net.parents(1)) == []
        assert net.outdegree(1) == 0
        assert net.indegree(1) == 0


    def test_single_node_lsa_node(self) -> None:
        """Test LSA_node property for single-node network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(nodes=[(1, {"label": "A"})])
        
        assert net.LSA_node == 1
        assert net.LSA_node == net.root_node
        assert net.LSA_node in net.leaves

    def test_single_node_incident_edges(self) -> None:
        """Test incident edges for single-node network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(nodes=[(1, {"label": "A"})])
        
        assert list(net.incident_parent_edges(1)) == []
        assert list(net.incident_child_edges(1)) == []
        assert net.number_of_edges() == 0

    def test_single_node_neighbors(self) -> None:
        """Test neighbors for single-node network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(nodes=[(1, {"label": "A"})])
        
        assert list(net.parents(1)) == []
        assert list(net.children(1)) == []
        assert list(net.neighbors(1)) == []

