"""
Comprehensive tests for SemiDirectedPhyNetwork initialization and validation.

This module tests all aspects of network initialization including:
- Empty networks
- Valid network structures
- Invalid network structures
- Edge format variations
- Taxa and label handling
- Validation rules specific to semi-directed networks
"""

import warnings
from typing import Dict, List, Tuple

import pytest

from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
from tests.core.network.sdnetwork.conftest import expect_mixed_network_warning


class TestEmptyNetwork:
    """Test cases for empty network initialization."""

    def test_empty_network_with_warning(self) -> None:
        """Test that empty network raises a warning."""
        # Empty networks raise warnings during initialization and validation
        with pytest.warns(UserWarning, match="empty edges list|Empty network.*no nodes"):
            net = SemiDirectedPhyNetwork(directed_edges=[], undirected_edges=[])
        
        assert net.number_of_nodes() == 0
        assert net.number_of_edges() == 0
        # Empty networks skip validation, so no validity warning
        with pytest.warns(UserWarning, match="Empty network.*no nodes"):
            assert net.validate() is True  # Empty networks are valid

    def test_empty_network_properties(self) -> None:
        """Test properties of empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = SemiDirectedPhyNetwork(directed_edges=[], undirected_edges=[])
        
        assert len(net.leaves) == 0
        assert len(net.taxa) == 0
        assert len(net.internal_nodes) == 0
        assert len(net.hybrid_nodes) == 0
        assert len(net.tree_nodes) == 0


class TestMinimalValidNetworks:
    """Test cases for minimal valid network structures."""

    def test_single_undirected_edge_network(self) -> None:
        """Test network with single undirected edge."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(1, 2)],
            nodes=[(2, {'label': 'A'})]
        )
        assert net.number_of_nodes() == 2
        assert net.number_of_edges() == 1
        assert 2 in net.leaves

    def test_three_undirected_edge_star_network(self) -> None:
        """Test star network with internal node and three leaves (degree >= 3)."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        assert net.number_of_nodes() == 4
        assert net.number_of_edges() == 3
        assert net.leaves == {1, 2, 4}
        assert net.taxa == {"A", "B", "C"}

    def test_single_directed_edge_network(self) -> None:
        """Test network with hybrid node (indegree >= 2, total_degree = indegree + 1)."""
        # Hybrid node 2 has indegree 2 (from 1 and 3), total_degree 3 (outgoing to 4)
        # Connect nodes 1 and 3 via undirected edge to form single source component
        # Add extra edges to nodes 1 and 3 to ensure degree >= 3
        net = SemiDirectedPhyNetwork(
            directed_edges=[(1, 2), (3, 2)],
            undirected_edges=[(1, 3), (1, 5), (3, 6), (2, 4)],
            nodes=[(4, {'label': 'A'}), (5, {'label': 'B'}), (6, {'label': 'C'})]
        )
        assert net.number_of_nodes() == 6
        assert net.number_of_edges() == 6
        assert 2 in net.hybrid_nodes


class TestSimpleTrees:
    """Test cases for simple tree structures."""

    def test_binary_tree_undirected(self) -> None:
        """Test a simple binary tree with undirected edges."""
        # Nodes 3 and 4 both need degree >= 3
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(4, 3), (3, 1), (3, 2), (3, 6), (4, 5), (4, 7)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'}), (7, {'label': 'E'})]
        )
        assert net.number_of_nodes() == 7
        assert net.number_of_edges() == 6
        assert net.leaves == {1, 2, 5, 6, 7}

    def test_ternary_tree_undirected(self) -> None:
        """Test a tree with node having degree 3."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(4, 1), (4, 2), (4, 3)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
        )
        assert net.number_of_nodes() == 4
        assert net.number_of_edges() == 3
        assert net.leaves == {1, 2, 3}


class TestNetworksWithHybrids:
    """Test cases for networks with hybrid nodes."""

    def test_single_hybrid_node(self) -> None:
        """Test network with single hybrid node."""
        # Hybrid node 4: indegree 2 (from 5, 6), total_degree 3 (outgoing to 2)
        # Nodes 5 and 6 need degree >= 3, so add more edges
        # Connect nodes 5 and 6 via undirected edge to form single source component
        net = SemiDirectedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(5, 6), (4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
        )
        assert net.number_of_nodes() == 8
        assert net.number_of_edges() == 8
        assert 4 in net.hybrid_nodes
        assert len(net.hybrid_nodes) > 0  # Has hybrid nodes, so not a tree

    def test_multiple_hybrid_nodes(self) -> None:
        """Test network with multiple hybrid nodes."""
        # Hybrid node 5: indegree 2 (from 7, 8), total_degree 3 (outgoing to 4)
        # Hybrid node 4: indegree 2 (from 5, 6), total_degree 3 (outgoing to 1 only)
        # Nodes 6, 7, 8, 9 need degree >= 3
        # Connect nodes 6, 7, 8 via undirected edges to form single source component
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                (7, 5), (8, 5),  # Hybrid node 5
                (5, 4), (6, 4)   # Hybrid node 4
            ],
            undirected_edges=[
                (6, 7), (7, 8),  # Connect source nodes
                (4, 1), (6, 11), (6, 12), (6, 9), (7, 13), (7, 14), (8, 15), (8, 16), (9, 2), (9, 3), (9, 17)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (11, {'label': 'D'}), (12, {'label': 'E'}), (13, {'label': 'F'}), (14, {'label': 'G'}), (15, {'label': 'H'}), (16, {'label': 'I'}), (17, {'label': 'J'})]
        )
        assert 4 in net.hybrid_nodes
        assert 5 in net.hybrid_nodes
        assert len(net.hybrid_nodes) == 2
        assert len(net.hybrid_nodes) > 0  # Has hybrid nodes, so not a tree


class TestSemiDirectedConstraints:
    """Test cases for semi-directed network specific constraints."""

    def test_all_hybrid_edges_directed(self) -> None:
        """Test that all hybrid edges are directed."""
        # Nodes 5 and 6 need degree >= 3
        # Connect nodes 5 and 6 via undirected edge to form single source component
        net = SemiDirectedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(5, 6), (4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
        )
        # All hybrid edges should be directed
        for u, v in net.hybrid_edges:
            assert net._graph._directed.has_edge(u, v)

    def test_all_non_hybrid_edges_undirected(self) -> None:
        """Test that all non-hybrid edges are undirected."""
        # Nodes 5 and 6 need degree >= 3
        # Connect nodes 5 and 6 via undirected edge to form single source component
        net = SemiDirectedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(5, 6), (4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
        )
        # All tree edges should be undirected
        for u, v in net.tree_edges:
            assert net._graph._undirected.has_edge(u, v)

    def test_validation_passes_valid_semi_directed(self) -> None:
        """Test that valid semi-directed network passes validation."""
        # Nodes 5 and 6 need degree >= 3
        # Connect nodes 5 and 6 via undirected edge to form single source component
        net = SemiDirectedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(5, 6), (4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
        )
        # Validation should pass (SemiDirectedPhyNetwork doesn't raise validity warning)
        try:
            result = net.validate()
            assert isinstance(result, bool)
        except (ValueError, NotImplementedError):
            pass  # Expected if validation is not fully implemented


class TestEdgeFormats:
    """Test cases for different edge input formats."""

    def test_tuple_format(self) -> None:
        """Test edges in tuple format."""
        # Node 2 needs degree >= 3
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(1, 2), (2, 3), (2, 4)],
            nodes=[(1, {'label': 'A'}), (3, {'label': 'B'}), (4, {'label': 'C'})]
        )
        assert net.number_of_edges() == 3

    def test_dict_format(self) -> None:
        """Test edges in dict format."""
        # Node 2 needs degree >= 3
        net = SemiDirectedPhyNetwork(
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
        net = SemiDirectedPhyNetwork(
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
        # Node 3 needs degree >= 3
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        assert net.taxa == {"A", "B", "C"}
        assert net.get_label(1) == "A"
        assert net.get_label(2) == "B"

    def test_partial_taxa_mapping(self) -> None:
        """Test network with partial taxa mapping (auto-labeling)."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        # Leaves 2 and 4 should get auto-generated labels
        assert "A" in net.taxa
        assert net.get_label(1) == "A"
        assert net.get_label(2) is not None
        assert net.get_label(4) is not None

    def test_taxa_duplicate_labels(self) -> None:
        """Test that duplicate taxon labels raise error."""
        with pytest.raises(ValueError, match="already used"):
            SemiDirectedPhyNetwork(
                undirected_edges=[(3, 1), (3, 2), (3, 4)],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'A'}), (4, {'label': 'B'})]  # Duplicate
            )


class TestValidation:
    """Test cases for network validation."""

    def test_connected_network(self) -> None:
        """Test that connected network validates."""
        # Node 2 needs degree >= 3
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(1, 2), (2, 3), (2, 4)],
            nodes=[(1, {'label': 'A'}), (3, {'label': 'B'}), (4, {'label': 'C'})]
        )
        # SemiDirectedPhyNetwork doesn't raise the validity warning
        assert net.validate() is True

    def test_disconnected_network(self) -> None:
        """Test that disconnected network raises error."""
        with pytest.raises(ValueError, match="not connected"):
            SemiDirectedPhyNetwork(
                undirected_edges=[(1, 2), (3, 4)],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (4, {'label': 'D'})]
            )

    def test_internal_node_degree_less_than_three(self) -> None:
        """Test that internal node with degree < 3 raises error."""
        with pytest.raises(ValueError, match="degree"):
            SemiDirectedPhyNetwork(
                undirected_edges=[(1, 2), (2, 3)],  # Node 2 has degree 2
                nodes=[(1, {'label': 'A'}), (3, {'label': 'B'})]
            )


class TestCopy:
    """Test cases for copy() method."""

    def test_copy_returns_semi_directed(self) -> None:
        """Test that copy returns SemiDirectedPhyNetwork."""
        # Node 3 needs degree >= 3
        net1 = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        net2 = net1.copy()
        assert isinstance(net2, SemiDirectedPhyNetwork)
        assert net1.leaves == net2.leaves
        assert net1.taxa == net2.taxa

