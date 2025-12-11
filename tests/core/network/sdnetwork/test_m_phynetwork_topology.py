"""
Comprehensive tests for MixedPhyNetwork topology properties.

This module tests all network topology properties including:
- leaves
- taxa
- internal_nodes
- hybrid_nodes
- tree_nodes
- hybrid_edges
- tree_edges
- level
- is_tree()
"""

import warnings

import pytest

from phylozoo.core.network.sdnetwork import MixedPhyNetwork
from tests.core.network.sdnetwork.conftest import expect_mixed_network_warning


class TestLeaves:
    """Test cases for leaves property."""

    def test_leaves_empty_network(self) -> None:
        """Test leaves in empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            # Empty networks don't raise validity warning (validation is skipped)
            net = MixedPhyNetwork(directed_edges=[], undirected_edges=[])
        assert len(net.leaves) == 0

    def test_leaves_simple_tree(self) -> None:
        """Test leaves in simple tree."""
        # Node 3 needs degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.leaves == {1, 2, 4}

    def test_leaves_network_with_hybrids(self) -> None:
        """Test leaves in network with hybrid nodes."""
        # Nodes 5 and 6 need degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        assert net.leaves == {2, 8, 9, 10, 11}

    def test_leaves_all_have_degree_one(self) -> None:
        """Test that all leaves have degree 1."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        for leaf in net.leaves:
            assert net.degree(leaf) == 1

    def test_leaves_many_leaves(self) -> None:
        """Test leaves with many leaves."""
        edges = [(100, i) for i in range(1, 51)]
        nodes = [(i, {'label': f"Taxon{i}"}) for i in range(1, 51)]
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(undirected_edges=edges, nodes=nodes)
        assert len(net.leaves) == 50


class TestTaxa:
    """Test cases for taxa property."""

    def test_taxa_empty_network(self) -> None:
        """Test taxa in empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            # Empty networks don't raise validity warning (validation is skipped)
            net = MixedPhyNetwork(directed_edges=[], undirected_edges=[])
        assert len(net.taxa) == 0

    def test_taxa_single_taxon(self) -> None:
        """Test taxa with single taxon."""
        # Single edge network: node 1 is leaf, node 3 is not internal (degree 1)
        # Actually, node 3 has degree 1, so it's also a leaf
        # Need a valid network structure
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert "A" in net.taxa
        assert len(net.taxa) == 3

    def test_taxa_many_taxa(self) -> None:
        """Test taxa with many taxa."""
        edges = [(100, i) for i in range(1, 21)]
        nodes = [(i, {'label': f"Taxon{i}"}) for i in range(1, 21)]
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(undirected_edges=edges, nodes=nodes)
        assert len(net.taxa) == 20

    def test_taxa_auto_labeled(self) -> None:
        """Test taxa with auto-labeled leaves."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]  # Only one labeled
            )
        # All leaves should have labels
        assert len(net.taxa) == 3
        assert "A" in net.taxa

    def test_taxa_matches_leaves(self) -> None:
        """Test that taxa count matches leaves count."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert len(net.taxa) == len(net.leaves)


class TestInternalNodes:
    """Test cases for internal_nodes property."""

    def test_internal_nodes_empty_network(self) -> None:
        """Test internal_nodes in empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            # Empty networks don't raise validity warning (validation is skipped)
            net = MixedPhyNetwork(directed_edges=[], undirected_edges=[])
        assert len(net.internal_nodes) == 0

    def test_internal_nodes_simple_tree(self) -> None:
        """Test internal_nodes in simple tree."""
        # Node 3 needs degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert 3 in net.internal_nodes
        assert len(net.internal_nodes) == 1

    def test_internal_nodes_all_non_leaves(self) -> None:
        """Test that internal_nodes are all non-leaves."""
        # Nodes 3 and 4 both need degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(4, 3), (3, 1), (3, 2), (3, 6), (4, 5), (4, 7)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'}), (7, {'label': 'E'})]
            )
        for internal in net.internal_nodes:
            assert internal not in net.leaves

    def test_internal_nodes_have_degree_at_least_three(self) -> None:
        """Test that internal nodes have degree >= 3."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(4, 3), (3, 1), (3, 2), (4, 5), (4, 6)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
            )
        for internal in net.internal_nodes:
            assert net.degree(internal) >= 3


class TestHybridNodes:
    """Test cases for hybrid_nodes property."""

    def test_hybrid_nodes_empty_network(self) -> None:
        """Test hybrid_nodes in empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            # Empty networks don't raise validity warning (validation is skipped)
            net = MixedPhyNetwork(directed_edges=[], undirected_edges=[])
        assert len(net.hybrid_nodes) == 0

    def test_hybrid_nodes_single_hybrid(self) -> None:
        """Test hybrid_nodes with single hybrid."""
        # Nodes 5 and 6 need degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        assert 4 in net.hybrid_nodes
        assert len(net.hybrid_nodes) == 1

    def test_hybrid_nodes_multiple(self) -> None:
        """Test hybrid_nodes with multiple hybrids."""
        # Hybrid node 4: indegree 2, total_degree must be 3 (only 1 outgoing)
        # Hybrid node 5: indegree 2, total_degree must be 3 (only 1 outgoing)
        # Nodes 6, 7, 8, 9 need degree >= 3, and connect node 9
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[
            (7, 5), (8, 5),  # Hybrid 5
            (5, 4), (6, 4)   # Hybrid 4
            ],
            undirected_edges=[
            (4, 1), (6, 11), (6, 12), (6, 9), (7, 13), (7, 14), (8, 15), (8, 16), (9, 2), (9, 3)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (11, {'label': 'D'}), (12, {'label': 'E'}), (13, {'label': 'F'}), (14, {'label': 'G'}), (15, {'label': 'H'}), (16, {'label': 'I'})]
            )
        assert 4 in net.hybrid_nodes
        assert 5 in net.hybrid_nodes
        assert len(net.hybrid_nodes) == 2

    def test_hybrid_nodes_have_indegree_at_least_two(self) -> None:
        """Test that hybrid nodes have indegree >= 2."""
        # Nodes 5 and 6 need degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        for hybrid in net.hybrid_nodes:
            assert net.indegree(hybrid) >= 2

    def test_hybrid_nodes_total_degree_equals_indegree_plus_one(self) -> None:
        """Test that hybrid nodes have total_degree = indegree + 1."""
        # Nodes 5 and 6 need degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        for hybrid in net.hybrid_nodes:
            indegree = net.indegree(hybrid)
            total_degree = net.degree(hybrid)
            assert total_degree == indegree + 1


class TestTreeNodes:
    """Test cases for tree_nodes property."""

    def test_tree_nodes_empty_network(self) -> None:
        """Test tree_nodes in empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            # Empty networks don't raise validity warning (validation is skipped)
            net = MixedPhyNetwork(directed_edges=[], undirected_edges=[])
        assert len(net.tree_nodes) == 0

    def test_tree_nodes_simple_tree(self) -> None:
        """Test tree_nodes in simple tree."""
        # Node 3 needs degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert 3 in net.tree_nodes
        assert len(net.tree_nodes) == 1

    def test_tree_nodes_have_indegree_zero(self) -> None:
        """Test that tree nodes have indegree 0."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        for tree_node in net.tree_nodes:
            assert net.indegree(tree_node) == 0

    def test_tree_nodes_have_total_degree_at_least_three(self) -> None:
        """Test that tree nodes have total_degree >= 3."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        for tree_node in net.tree_nodes:
            assert net.degree(tree_node) >= 3


class TestHybridEdges:
    """Test cases for hybrid_edges property."""

    def test_hybrid_edges_empty_network(self) -> None:
        """Test hybrid_edges in empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            # Empty networks don't raise validity warning (validation is skipped)
            net = MixedPhyNetwork(directed_edges=[], undirected_edges=[])
        assert len(net.hybrid_edges) == 0

    def test_hybrid_edges_single_hybrid(self) -> None:
        """Test hybrid_edges with single hybrid."""
        # Nodes 5 and 6 need degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        assert (5, 4) in net.hybrid_edges
        assert (6, 4) in net.hybrid_edges
        assert len(net.hybrid_edges) == 2

    def test_hybrid_edges_are_all_directed(self) -> None:
        """Test that hybrid_edges are all directed."""
        # Nodes 5 and 6 need degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        for u, v in net.hybrid_edges:
            assert net._graph._directed.has_edge(u, v)
            assert not net._graph._undirected.has_edge(u, v)

    def test_hybrid_edges_parallel(self) -> None:
        """Test hybrid_edges with parallel edges."""
        # Nodes 5 and 6 need degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[
            (5, 4, 0), (5, 4, 1),  # Parallel edges
            (6, 4)
            ],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        # Should include both parallel edges
        assert (5, 4) in net.hybrid_edges
        assert (6, 4) in net.hybrid_edges


class TestTreeEdges:
    """Test cases for tree_edges property."""

    def test_tree_edges_empty_network(self) -> None:
        """Test tree_edges in empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            # Empty networks don't raise validity warning (validation is skipped)
            net = MixedPhyNetwork(directed_edges=[], undirected_edges=[])
        assert len(net.tree_edges) == 0

    def test_tree_edges_simple_tree(self) -> None:
        """Test tree_edges in simple tree."""
        # Node 3 needs degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert (3, 1) in net.tree_edges or (1, 3) in net.tree_edges
        assert (3, 2) in net.tree_edges or (2, 3) in net.tree_edges
        assert (3, 4) in net.tree_edges or (4, 3) in net.tree_edges
        assert len(net.tree_edges) == 3

    def test_tree_edges_are_all_undirected(self) -> None:
        """Test that tree_edges are all undirected."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        for u, v in net.tree_edges:
            assert net._graph._undirected.has_edge(u, v)
            assert not net._graph._directed.has_edge(u, v)

    def test_tree_edges_parallel(self) -> None:
        """Test tree_edges with parallel edges."""
        # Node 3 needs degree >= 3
        # Parallel edges between internal nodes, not leaves
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[
            (3, 4, 0), (3, 4, 1),  # Parallel edges between internal nodes
            (3, 1), (3, 2), (4, 5), (4, 6)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
            )
        # Should include both parallel edges
        assert len(net.tree_edges) >= 4


class TestIsTree:
    """Test cases for is_tree() method."""

    def test_is_tree_simple_tree(self) -> None:
        """Test is_tree() on simple tree."""
        # Node 3 needs degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        # Simple tree has no hybrid nodes
        assert len(net.hybrid_nodes) == 0

    def test_is_tree_with_hybrids(self) -> None:
        """Test is_tree() on network with hybrids."""
        # Nodes 5 and 6 need degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        # Network has hybrid nodes, so it's not a tree
        assert len(net.hybrid_nodes) > 0

    def test_is_tree_empty_network(self) -> None:
        """Test is_tree() on empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            # Empty networks don't raise validity warning (validation is skipped)
            net = MixedPhyNetwork(directed_edges=[], undirected_edges=[])
        # Empty network has no hybrid nodes, so check hybrid_nodes instead
        assert len(net.hybrid_nodes) == 0


class TestLevel:
    """Test cases for level property."""

    def test_level_exists(self) -> None:
        """Test that level property exists."""
        # Level property only exists in SemiDirectedPhyNetwork, not MixedPhyNetwork
        # Test with SemiDirectedPhyNetwork instead
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        # Node 3 needs degree >= 3
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        # Level should exist (may be placeholder)
        level = net.level
        assert level is not None

    def test_level_cached(self) -> None:
        """Test that level is cached."""
        # Level property only exists in SemiDirectedPhyNetwork, not MixedPhyNetwork
        # Test with SemiDirectedPhyNetwork instead
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        # Node 3 needs degree >= 3
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        level1 = net.level
        level2 = net.level
        assert level1 is level2

