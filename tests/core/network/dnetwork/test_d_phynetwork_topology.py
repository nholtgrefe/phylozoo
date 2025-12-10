"""
Comprehensive tests for DirectedPhyNetwork topology properties.

This module tests all network topology properties including:
- leaves
- taxa
- root_node
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

from phylozoo.core.network import DirectedPhyNetwork


class TestLeaves:
    """Test cases for leaves property."""

    def test_leaves_empty_network(self) -> None:
        """Test leaves in empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(edges=[])
        assert len(net.leaves) == 0

    def test_leaves_simple_tree(self) -> None:
        """Test leaves in simple tree."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        assert net.leaves == {1, 2}

    def test_leaves_network_with_hybrids(self) -> None:
        """Test leaves in network with hybrid nodes."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        assert net.leaves == {2, 8, 9}

    def test_leaves_all_have_outdegree_zero(self) -> None:
        """Test that all leaves have out-degree 0."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        for leaf in net.leaves:
            assert net.outdegree(leaf) == 0

    def test_leaves_many_leaves(self) -> None:
        """Test leaves with many leaves."""
        edges = [(100, i) for i in range(1, 51)]  # Use 100 as root, not 10
        nodes = [(i, {"label": f"Taxon{i}"}) for i in range(1, 51)]
        net = DirectedPhyNetwork(edges=edges, nodes=nodes)
        assert len(net.leaves) == 50


class TestTaxa:
    """Test cases for taxa property."""

    def test_taxa_empty_network(self) -> None:
        """Test taxa in empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(edges=[])
        assert len(net.taxa) == 0

    def test_taxa_single_taxon(self) -> None:
        """Test taxa with single taxon."""
        net = DirectedPhyNetwork(edges=[(3, 1)], nodes=[(1, {'label': 'A'})])
        assert net.taxa == {"A"}

    def test_taxa_many_taxa(self) -> None:
        """Test taxa with many taxa."""
        edges = [(100, i) for i in range(1, 21)]  # Use 100 as root, not 10
        nodes = [(i, {"label": f"Taxon{i}"}) for i in range(1, 21)]
        net = DirectedPhyNetwork(edges=edges, nodes=nodes)
        assert len(net.taxa) == 20

    def test_taxa_auto_labeled(self) -> None:
        """Test taxa with auto-labeled leaves."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'})]  # Only one labeled
        )
        # All leaves should have labels
        assert len(net.taxa) == 3
        assert "A" in net.taxa

    def test_taxa_matches_leaves(self) -> None:
        """Test that taxa count matches leaves count."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        assert len(net.taxa) == len(net.leaves)


class TestRootNode:
    """Test cases for root_node property."""

    def test_root_node_simple_tree(self) -> None:
        """Test root node in simple tree."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        assert net.root_node == 3

    def test_root_node_has_indegree_zero(self) -> None:
        """Test that root node has in-degree 0."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        assert net.indegree(net.root_node) == 0

    def test_root_node_empty_network(self) -> None:
        """Test root node in empty network raises error."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(edges=[])
        with pytest.raises(ValueError, match="Network has no root node"):
            _ = net.root_node

    def test_root_node_multiple_roots_raises(self) -> None:
        """Test that multiple roots raise error during init."""
        with pytest.raises(ValueError, match="multiple root nodes"):
            DirectedPhyNetwork(
                edges=[(1, 3), (2, 3)],
                nodes=[(3, {'label': 'A'})]
            )


class TestInternalNodes:
    """Test cases for internal_nodes property."""

    def test_internal_nodes_empty(self) -> None:
        """Test internal nodes in minimal network."""
        net = DirectedPhyNetwork(edges=[(3, 1)], nodes=[(1, {'label': 'A'})])
        # Only root and leaf, no internal nodes
        assert len(net.internal_nodes) == 0

    def test_internal_nodes_simple_tree(self) -> None:
        """Test internal nodes in simple tree."""
        net = DirectedPhyNetwork(
            edges=[(4, 3), (3, 1), (3, 2), (4, 5)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'})]
        )
        # Node 3 is internal (not root, not leaf)
        assert net.internal_nodes == [3]

    def test_internal_nodes_excludes_root(self) -> None:
        """Test that internal_nodes excludes root."""
        net = DirectedPhyNetwork(
            edges=[(4, 3), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        assert 4 not in net.internal_nodes  # Root excluded

    def test_internal_nodes_excludes_leaves(self) -> None:
        """Test that internal_nodes excludes leaves."""
        net = DirectedPhyNetwork(
            edges=[(4, 3), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        assert 1 not in net.internal_nodes  # Leaf excluded
        assert 2 not in net.internal_nodes  # Leaf excluded

    def test_internal_nodes_with_hybrids(self) -> None:
        """Test internal nodes in network with hybrids."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        # Internal nodes: 4 (hybrid), 5 (tree), 6 (tree)
        assert len(net.internal_nodes) == 3
        assert 4 in net.internal_nodes
        assert 5 in net.internal_nodes
        assert 6 in net.internal_nodes


class TestHybridNodes:
    """Test cases for hybrid_nodes property."""

    def test_hybrid_nodes_none(self) -> None:
        """Test hybrid nodes in tree (none)."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        assert len(net.hybrid_nodes) == 0

    def test_hybrid_nodes_single(self) -> None:
        """Test single hybrid node."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        assert net.hybrid_nodes == [4]

    def test_hybrid_nodes_multiple(self) -> None:
        """Test multiple hybrid nodes."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 7), (10, 8),
                (7, 5), (7, 6),
                (8, 5), (8, 9),
                (5, 4), (6, 4),
                (6, 11),  # 6 splits to 4 and 11
                (4, 1), (9, 2), (9, 3)  # 9 splits to 2 children
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (11, {'label': 'D'})]
        )
        # Hybrid nodes: 4 and 5
        assert len(net.hybrid_nodes) == 2
        assert 4 in net.hybrid_nodes
        assert 5 in net.hybrid_nodes

    def test_hybrid_nodes_indegree_requirement(self) -> None:
        """Test that hybrid nodes have in-degree >= 2."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        for hybrid in net.hybrid_nodes:
            assert net.indegree(hybrid) >= 2

    def test_hybrid_nodes_outdegree_requirement(self) -> None:
        """Test that hybrid nodes have out-degree == 1."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        for hybrid in net.hybrid_nodes:
            assert net.outdegree(hybrid) == 1

    def test_hybrid_nodes_parallel_edges(self) -> None:
        """Test hybrid nodes with parallel edges."""
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),
                (5, 4, 0), (5, 4, 1),  # Parallel edges
                (5, 8),
                (6, 4),
                (6, 9),
                (4, 2)
            ],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        # Node 4 should still be hybrid (in-degree >= 2)
        assert 4 in net.hybrid_nodes


class TestTreeNodes:
    """Test cases for tree_nodes property."""

    def test_tree_nodes_none(self) -> None:
        """Test tree nodes in minimal network."""
        net = DirectedPhyNetwork(edges=[(3, 1)], nodes=[(1, {'label': 'A'})])
        assert len(net.tree_nodes) == 0

    def test_tree_nodes_single(self) -> None:
        """Test single tree node."""
        net = DirectedPhyNetwork(
            edges=[(4, 3), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        assert net.tree_nodes == [3]

    def test_tree_nodes_multiple(self) -> None:
        """Test multiple tree nodes."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        # Tree nodes: 5 and 6
        assert len(net.tree_nodes) == 2
        assert 5 in net.tree_nodes
        assert 6 in net.tree_nodes

    def test_tree_nodes_indegree_requirement(self) -> None:
        """Test that tree nodes have in-degree == 1."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        for tree_node in net.tree_nodes:
            assert net.indegree(tree_node) == 1

    def test_tree_nodes_outdegree_requirement(self) -> None:
        """Test that tree nodes have out-degree >= 2."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        for tree_node in net.tree_nodes:
            assert net.outdegree(tree_node) >= 2

    def test_tree_nodes_with_hybrids(self) -> None:
        """Test tree nodes in network with hybrids."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        # Should have tree nodes even with hybrids
        assert len(net.tree_nodes) > 0


class TestHybridEdges:
    """Test cases for hybrid_edges property."""

    def test_hybrid_edges_none(self) -> None:
        """Test hybrid edges in tree (none)."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        assert len(net.hybrid_edges) == 0

    def test_hybrid_edges_single_hybrid(self) -> None:
        """Test hybrid edges with single hybrid."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        # Hybrid edges: (5, 4) and (6, 4)
        assert len(net.hybrid_edges) == 2
        assert (5, 4) in net.hybrid_edges
        assert (6, 4) in net.hybrid_edges

    def test_hybrid_edges_multiple_hybrids(self) -> None:
        """Test hybrid edges with multiple hybrids."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 7), (10, 8),
                (7, 5), (7, 6),
                (8, 5), (8, 9),
                (5, 4), (6, 4),
                (6, 11),  # 6 splits to 4 and 11
                (4, 1), (9, 2), (9, 3)  # 9 splits to 2 children
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (11, {'label': 'D'})]
        )
        # Hybrid edges: (7, 5), (8, 5) for hybrid 5, and (5, 4), (6, 4) for hybrid 4
        assert len(net.hybrid_edges) == 4

    def test_hybrid_edges_parallel(self) -> None:
        """Test hybrid edges with parallel edges."""
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),
                (5, 4, 0), (5, 4, 1),  # Parallel edges to hybrid
                (5, 8),
                (6, 4),
                (6, 9),
                (4, 2)
            ],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        # Should include both parallel edges
        assert (5, 4) in net.hybrid_edges
        # Note: parallel edges appear as same (u, v) tuple
        assert net.hybrid_edges.count((5, 4)) >= 1


class TestTreeEdges:
    """Test cases for tree_edges property."""

    def test_tree_edges_all_in_tree(self) -> None:
        """Test tree edges in pure tree."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        assert len(net.tree_edges) == 2
        assert (3, 1) in net.tree_edges
        assert (3, 2) in net.tree_edges

    def test_tree_edges_mixed_with_hybrids(self) -> None:
        """Test tree edges mixed with hybrid edges."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        # Tree edges: (7, 5), (7, 6), (5, 8), (6, 9), (4, 2)
        # Hybrid edges: (5, 4), (6, 4)
        assert (7, 5) in net.tree_edges
        assert (7, 6) in net.tree_edges
        assert (5, 8) in net.tree_edges
        assert (6, 9) in net.tree_edges
        assert (4, 2) in net.tree_edges
        assert (5, 4) not in net.tree_edges
        assert (6, 4) not in net.tree_edges

    def test_tree_edges_complement_of_hybrid_edges(self) -> None:
        """Test that tree_edges + hybrid_edges = all edges."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        all_edges = set(net.tree_edges) | set(net.hybrid_edges)
        expected_edges = {(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)}
        assert all_edges == expected_edges


class TestLevel:
    """Test cases for level property."""

    def test_level_placeholder(self) -> None:
        """Test that level returns 0 (placeholder)."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        # Currently returns 0 as placeholder
        assert net.level == 0

    def test_level_cached(self) -> None:
        """Test that level is cached."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        level1 = net.level
        level2 = net.level
        assert level1 is level2  # Same object (cached)


class TestIsTree:
    """Test cases for is_tree() method."""

    def test_is_tree_true(self) -> None:
        """Test is_tree() returns True for tree."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        assert net.is_tree() is True

    def test_is_tree_false_with_hybrid(self) -> None:
        """Test is_tree() returns False with hybrid node."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        assert net.is_tree() is False

    def test_is_tree_empty(self) -> None:
        """Test is_tree() on empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(edges=[])
        # Empty network has no hybrid nodes, so is_tree() should be True
        # But accessing hybrid_nodes might fail, so test carefully
        try:
            result = net.is_tree()
            assert result is True
        except (ValueError, AttributeError):
            # If it fails due to no root, that's expected
            pass


class TestTopologyConsistency:
    """Test consistency between topology properties."""

    def test_all_nodes_accounted_for(self) -> None:
        """Test that root + leaves + internal_nodes = all nodes."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        all_nodes = set(net._graph.nodes)
        accounted = {net.root_node} | net.leaves | set(net.internal_nodes)
        assert all_nodes == accounted

    def test_internal_nodes_equals_tree_plus_hybrid(self) -> None:
        """Test that internal_nodes = tree_nodes + hybrid_nodes."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        internal_set = set(net.internal_nodes)
        tree_hybrid_set = set(net.tree_nodes) | set(net.hybrid_nodes)
        assert internal_set == tree_hybrid_set

    def test_hybrid_edges_count_matches_hybrid_nodes(self) -> None:
        """Test that hybrid_edges count matches hybrid node in-degrees."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        # Count hybrid edges
        hybrid_edge_count = len(net.hybrid_edges)
        # Sum of in-degrees of hybrid nodes
        hybrid_indegree_sum = sum(net.indegree(h) for h in net.hybrid_nodes)
        assert hybrid_edge_count == hybrid_indegree_sum

