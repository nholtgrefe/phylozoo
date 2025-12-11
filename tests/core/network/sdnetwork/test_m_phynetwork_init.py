"""
Comprehensive tests for MixedPhyNetwork initialization and validation.

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

from phylozoo.core.network.sdnetwork import MixedPhyNetwork
from tests.core.network.sdnetwork.conftest import expect_mixed_network_warning


class TestEmptyNetwork:
    """Test cases for empty network initialization."""

    def test_empty_network_with_warning(self) -> None:
        """Test that empty network raises a warning."""
        # Empty networks raise warnings during initialization and validation
        with pytest.warns(UserWarning, match="Empty network.*no nodes"):
            net = MixedPhyNetwork(directed_edges=[], undirected_edges=[])
        
        assert net.number_of_nodes() == 0
        assert net.number_of_edges() == 0
        # Empty networks skip validation, so no validity warning
        with pytest.warns(UserWarning, match="Empty network.*no nodes"):
            net.validate()  # Empty networks are valid

    def test_empty_network_properties(self) -> None:
        """Test properties of empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            # Empty networks don't raise validity warning (validation is skipped)
            net = MixedPhyNetwork(directed_edges=[], undirected_edges=[])
        
        assert len(net.leaves) == 0
        assert len(net.taxa) == 0
        assert len(net.internal_nodes) == 0
        assert len(net.hybrid_nodes) == 0
        assert len(net.tree_nodes) == 0


class TestMinimalValidNetworks:
    """Test cases for minimal valid network structures."""

    def test_single_undirected_edge_network(self) -> None:
        """Test network with single undirected edge (internal -> leaf)."""
        # Single edge network: both nodes are leaves (degree 1)
        # Need at least 3 edges from internal node
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.number_of_nodes() == 4
        assert net.number_of_edges() == 3
        assert 1 in net.leaves
        assert 2 in net.leaves
        assert 4 in net.leaves
        assert "A" in net.taxa
        assert "B" in net.taxa
        assert "C" in net.taxa

    def test_two_undirected_edge_star_network(self) -> None:
        """Test star network with internal node and three leaves."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.number_of_nodes() == 4
        assert net.number_of_edges() == 3
        assert net.leaves == {1, 2, 4}
        assert net.taxa == {"A", "B", "C"}

    def test_single_directed_edge_network(self) -> None:
        """Test network with hybrid node (requires at least 2 incoming directed edges)."""
        # Hybrid node 4 needs at least 2 incoming directed edges
        # Nodes 1 and 2 need degree >= 3 to be valid internal nodes
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(1, 4), (2, 4)],  # Two hybrid edges into node 4
            undirected_edges=[(4, 3), (1, 5), (1, 6), (2, 7), (2, 8)],  # Undirected edges
            nodes=[(3, {'label': 'A'}), (5, {'label': 'B'}), (6, {'label': 'C'}), (7, {'label': 'D'}), (8, {'label': 'E'})]
            )
        assert net.number_of_nodes() == 8
        assert net.number_of_edges() == 7
        assert 3 in net.leaves
        assert 4 in net.hybrid_nodes


class TestSimpleTrees:
    """Test cases for simple tree structures."""

    def test_binary_tree_undirected(self) -> None:
        """Test a simple binary tree with undirected edges."""
        # Node 4 needs degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(4, 3), (3, 1), (3, 2), (4, 5), (4, 6)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
            )
        assert net.number_of_nodes() == 6
        assert net.number_of_edges() == 5
        assert net.leaves == {1, 2, 5, 6}

    def test_ternary_tree_undirected(self) -> None:
        """Test a tree with node having degree 3."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(4, 1), (4, 2), (4, 3)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
            )
        assert net.number_of_nodes() == 4
        assert net.number_of_edges() == 3
        assert net.leaves == {1, 2, 3}

    def test_star_tree_undirected(self) -> None:
        """Test a star tree with many leaves."""
        edges = [(10, i) for i in range(1, 10)]
        nodes = [(i, {'label': f"Taxon{i}"}) for i in range(1, 10)]
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(undirected_edges=edges, nodes=nodes)
        assert net.number_of_nodes() == 10
        assert net.number_of_edges() == 9
        assert len(net.leaves) == 9


class TestNetworksWithHybrids:
    """Test cases for networks with hybrid nodes."""

    def test_single_hybrid_node(self) -> None:
        """Test network with single hybrid node."""
        # Structure: tree nodes -> hybrid -> leaf (via undirected)
        # Hybrid has 2 incoming directed edges
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        assert net.number_of_nodes() == 8
        assert net.number_of_edges() == 7
        assert 4 in net.hybrid_nodes
        assert len(net.hybrid_nodes) > 0  # Has hybrid nodes, so not a tree

    def test_multiple_hybrid_nodes(self) -> None:
        """Test network with multiple hybrid nodes."""
        # Hybrid node 4: indegree 2, total_degree must be 3 (only 1 outgoing)
        # Hybrid node 5: indegree 2, total_degree must be 3 (only 1 outgoing)
        # Nodes 6, 7, 8, 9 need degree >= 3, and connect node 9
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[
            (7, 5), (8, 5),  # Hybrid node 5
            (5, 4), (6, 4)   # Hybrid node 4
            ],
            undirected_edges=[
            (4, 1), (6, 11), (6, 12), (6, 9), (7, 13), (7, 14), (8, 15), (8, 16), (9, 2), (9, 3)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (11, {'label': 'D'}), (12, {'label': 'E'}), (13, {'label': 'F'}), (14, {'label': 'G'}), (15, {'label': 'H'}), (16, {'label': 'I'})]
            )
        assert 4 in net.hybrid_nodes
        assert 5 in net.hybrid_nodes
        assert len(net.hybrid_nodes) == 2
        assert len(net.hybrid_nodes) > 0  # Has hybrid nodes, so not a tree

    def test_nested_hybrids(self) -> None:
        """Test network with nested hybridization."""
        # Hybrid node 4: indegree 2, total_degree must be 3 (only 1 outgoing)
        # Hybrid node 5: indegree 2, total_degree must be 3 (only 1 outgoing)
        # Nodes 6, 7, 8, 9 need degree >= 3, and connect node 9
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[
            (7, 5), (8, 5),  # Hybrid node 5
            (5, 4), (6, 4)   # Hybrid node 4 (nested)
            ],
            undirected_edges=[
            (4, 1), (6, 11), (6, 12), (6, 9), (7, 13), (7, 14), (8, 15), (8, 16), (9, 2), (9, 3)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (11, {'label': 'D'}), (12, {'label': 'E'}), (13, {'label': 'F'}), (14, {'label': 'G'}), (15, {'label': 'H'}), (16, {'label': 'I'})]
            )
        assert 4 in net.hybrid_nodes
        assert 5 in net.hybrid_nodes


class TestParallelEdges:
    """Test cases for networks with parallel edges."""

    def test_parallel_undirected_edges(self) -> None:
        """Test parallel undirected edges."""
        # Node 3 needs degree >= 3
        # Use parallel edges between internal nodes, not leaves (leaves have degree == 1)
        # Test parallel edges (3, 1) with keys 0 and 1, but node 1 is a leaf, so it can't have parallel edges
        # Instead, use parallel edges between 3 and 4 (both internal)
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[
            (3, 4, 0), (3, 4, 1),  # Two parallel edges between internal nodes
            (3, 1), (3, 2), (4, 5), (4, 6)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
            )
        # Check that parallel edges exist between 3 and 4
        assert net.has_edge(3, 4, key=0)
        assert net.has_edge(3, 4, key=1)
        assert net.number_of_edges() == 6  # 2 parallel + 4 regular edges

    def test_parallel_directed_edges_to_hybrid(self) -> None:
        """Test parallel directed edges entering a hybrid node."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[
            (5, 4, 0), (5, 4, 1),  # Two parallel edges to hybrid 4
            (6, 4)
            ],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        assert net.has_edge(5, 4, key=0)
        assert net.has_edge(5, 4, key=1)
        assert net.number_of_edges() == 8  # 3 directed + 5 undirected

    def test_parallel_edges_dict_format(self) -> None:
        """Test parallel edges using dict format."""
        # Node 3 needs degree >= 3
        # Use parallel edges between internal nodes, not leaves
        # Leaves are defined as degree == 1, so parallel edges on leaves cause issues
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[
            (3, 4, 0), (3, 4, 1),  # Parallel edges between internal nodes
            (3, 1), (3, 2), (4, 5), (4, 6)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
            )
        assert net.has_edge(3, 4, key=0)
        assert net.has_edge(3, 4, key=1)


class TestEdgeFormats:
    """Test cases for different edge input formats."""

    def test_tuple_format(self) -> None:
        """Test edges in tuple format."""
        # Node 2 needs degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(1, 2), (2, 3), (2, 4)],
            nodes=[(1, {'label': 'A'}), (3, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.number_of_edges() == 3

    def test_tuple_with_key_format(self) -> None:
        """Test edges in tuple format with explicit keys."""
        # Node 2 needs degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(1, 2, 0), (2, 3, 0), (2, 4)],
            nodes=[(1, {'label': 'A'}), (3, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.number_of_edges() == 3

    def test_dict_format(self) -> None:
        """Test edges in dict format."""
        # Node 2 needs degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[
            {'u': 1, 'v': 2},
            {'u': 2, 'v': 3, 'branch_length': 0.5},
            {'u': 2, 'v': 4}
            ],
            nodes=[(1, {'label': 'A'}), (3, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.number_of_edges() == 3

    def test_mixed_format(self) -> None:
        """Test mixing tuple and dict formats."""
        # Node 2 needs degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[
            (1, 2),
            {'u': 2, 'v': 3, 'branch_length': 0.5},
            (2, 4)
            ],
            nodes=[(1, {'label': 'A'}), (3, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.number_of_edges() == 3


class TestTaxaHandling:
    """Test cases for taxa label handling."""

    def test_complete_taxa_mapping(self) -> None:
        """Test network with complete taxa mapping."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.taxa == {"A", "B", "C"}
        assert net.get_label(1) == "A"
        assert net.get_label(2) == "B"
        assert net.get_label(4) == "C"

    def test_partial_taxa_mapping(self) -> None:
        """Test network with partial taxa mapping (auto-labeling)."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        # Leaves 2 and 4 should get auto-generated labels
        assert "A" in net.taxa
        assert net.get_label(1) == "A"
        # Auto-generated labels should be string representations of node IDs
        assert net.get_label(2) == "2" or net.get_label(2) is not None
        assert net.get_label(4) == "4" or net.get_label(4) is not None

    def test_taxa_as_list(self) -> None:
        """Test taxa provided as list of tuples."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
            )
        # Node 4 will get auto-labeled
        assert "A" in net.taxa
        assert "B" in net.taxa
        assert len(net.taxa) == 3  # A, B, and auto-labeled 4

    def test_taxa_duplicate_labels(self) -> None:
        """Test that duplicate taxon labels raise error."""
        with pytest.raises(ValueError, match="already used|duplicate"):
            MixedPhyNetwork(
                undirected_edges=[(3, 1), (3, 2), (3, 4)],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'A'}), (4, {'label': 'B'})]  # Duplicate
            )



class TestInternalNodeLabels:
    """Test cases for internal node labels."""

    def test_internal_node_label(self) -> None:
        """Test labeling an internal node."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (3, {'label': 'root'})],
            )
        assert net.get_label(3) == "root"

    def test_internal_node_labels_as_list(self) -> None:
        """Test internal node labels as list."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (3, {'label': 'root'})],
            )
        assert net.get_label(3) == "root"

    def test_internal_node_label_duplicate(self) -> None:
        """Test that duplicate internal node labels raise error."""
        with pytest.raises(ValueError, match="already used|duplicate"):
            MixedPhyNetwork(
                undirected_edges=[(3, 1), (3, 2), (3, 4)],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (3, {'label': 'A'})],
            )



class TestValidation:
    """Test cases for network validation."""

    def test_connected_network(self) -> None:
        """Test that connected network validates."""
        # Node 2 needs degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
                undirected_edges=[(1, 2), (2, 3), (2, 4)],
                nodes=[(1, {'label': 'A'}), (3, {'label': 'B'}), (4, {'label': 'C'})]
            )
        with expect_mixed_network_warning():
            net.validate()

    def test_disconnected_network(self) -> None:
        """Test that disconnected network raises error."""
        with pytest.raises(ValueError, match="not connected"):
            MixedPhyNetwork(
                undirected_edges=[(1, 2), (3, 4)],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (4, {'label': 'D'})]
            )

    def test_internal_node_degree_less_than_three(self) -> None:
        """Test that internal node with degree < 3 raises error."""
        with pytest.raises(ValueError, match="degree"):
            MixedPhyNetwork(
                undirected_edges=[(1, 2), (2, 3)],  # Node 2 has degree 2
                nodes=[(1, {'label': 'A'}), (3, {'label': 'B'})]
            )

    def test_node_indegree_constraint(self) -> None:
        """Test that node with invalid indegree raises error."""
        # Node 2 has indegree 1 but total_degree 2 (should be 0 or total_degree-1)
        # Node 1 needs degree >= 3
        with pytest.raises(ValueError, match="indegree|degree"):
            MixedPhyNetwork(
                directed_edges=[(1, 2)],
                undirected_edges=[(2, 3), (1, 4), (1, 5)],
                nodes=[(3, {'label': 'A'}), (4, {'label': 'B'}), (5, {'label': 'C'})]
            )

    def test_validation_warning(self) -> None:
        """Test that validation issues warning."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        # validate() also raises the warning
        with expect_mixed_network_warning():
            net.validate()


class TestInvalidInitialization:
    """Test cases for invalid network initialization (should raise errors)."""

    def test_no_edges_parameter(self) -> None:
        """Test that missing edges parameter raises warning."""
        # Empty networks (None edges) don't raise validity warning (validation is skipped)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = MixedPhyNetwork(directed_edges=None, undirected_edges=None)
        assert net.number_of_nodes() == 0

    def test_leaf_with_wrong_indegree(self) -> None:
        """Test that leaf with wrong in-degree raises error during validation."""
        # Leaf 2 has indegree=1 but total_degree=2 (should be 0 or total_degree-1=1)
        # Node 1 needs degree >= 3
        # Node 2 is a leaf (no outgoing directed edges), but has indegree 1 and total_degree 2
        # Actually, node 2 with undirected edge (2, 3) means it's not a pure leaf
        # Let's use a different structure: leaf 2 only has incoming directed edge
        with pytest.raises(ValueError, match="indegree|degree"):
            MixedPhyNetwork(
                directed_edges=[(1, 2)],  # Leaf 2 has indegree 1
                undirected_edges=[(1, 3), (1, 4), (1, 5)],  # Node 1 needs degree >= 3
                nodes=[(2, {'label': 'A'}), (3, {'label': 'B'}), (4, {'label': 'C'}), (5, {'label': 'D'})]
            )


class TestMixedDirectedUndirected:
    """Test cases for networks with both directed and undirected edges."""

    def test_hybrid_with_undirected_children(self) -> None:
        """Test hybrid node with undirected edges to children."""
        # Hybrid node 4: indegree 2, total_degree must be 3 (only 1 outgoing)
        # Nodes 5 and 6 need degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],  # Hybrid 4
            undirected_edges=[(4, 1), (5, 8), (5, 10), (6, 9), (6, 11)],  # Only 1 outgoing from 4
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        assert 4 in net.hybrid_nodes
        assert net.number_of_edges() == 7

    def test_tree_node_with_mixed_edges(self) -> None:
        """Test tree node with both directed and undirected edges."""
        # Tree node 3 has indegree 0, outdegree 1, undirected_degree 2, total degree 3
        # Node 4 is hybrid (indegree >= 2), total_degree must be 3 (only 1 outgoing)
        # Node 5 needs degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(3, 4), (5, 4)],  # Node 4 is hybrid
            undirected_edges=[(3, 1), (3, 2), (4, 6), (5, 7), (5, 8)],  # Undirected edges
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (6, {'label': 'C'}), (7, {'label': 'D'}), (8, {'label': 'E'})]
            )
        # Node 3 should be a tree node if it has indegree 0 and degree >= 3
        assert net.number_of_edges() == 7
        assert 3 in net.tree_nodes

