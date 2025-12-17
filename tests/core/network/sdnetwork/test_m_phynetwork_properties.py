"""
Comprehensive tests for MixedPhyNetwork cached properties.

This module tests cached property behavior including:
- Property caching (same object returned)
- Property invalidation on copy
- All cached properties
- Property consistency
"""

import warnings

import pytest

from phylozoo.core.network.sdnetwork import MixedPhyNetwork
from tests.core.network.sdnetwork.conftest import expect_mixed_network_warning


class TestCachedPropertyBehavior:
    """Test cases for cached property behavior."""

    def test_leaves_cached(self) -> None:
        """Test that leaves property is cached."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        leaves1 = net.leaves
        leaves2 = net.leaves
        assert leaves1 is leaves2  # Same object (cached)

    def test_taxa_cached(self) -> None:
        """Test that taxa property is cached."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        taxa1 = net.taxa
        taxa2 = net.taxa
        assert taxa1 is taxa2  # Same object (cached)

    def test_internal_nodes_cached(self) -> None:
        """Test that internal_nodes property is cached."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(4, 3), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        internal1 = net.internal_nodes
        internal2 = net.internal_nodes
        assert internal1 is internal2  # Same object (cached)

    def test_hybrid_nodes_cached(self) -> None:
        """Test that hybrid_nodes property is cached."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        hybrid1 = net.hybrid_nodes
        hybrid2 = net.hybrid_nodes
        assert hybrid1 is hybrid2  # Same object (cached)

    def test_tree_nodes_cached(self) -> None:
        """Test that tree_nodes property is cached."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        tree1 = net.tree_nodes
        tree2 = net.tree_nodes
        assert tree1 is tree2  # Same object (cached)

    def test_hybrid_edges_cached(self) -> None:
        """Test that hybrid_edges property is cached."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        hybrid_edges1 = net.hybrid_edges
        hybrid_edges2 = net.hybrid_edges
        assert hybrid_edges1 is hybrid_edges2  # Same object (cached)

    def test_tree_edges_cached(self) -> None:
        """Test that tree_edges property is cached."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        tree_edges1 = net.tree_edges
        tree_edges2 = net.tree_edges
        assert tree_edges1 is tree_edges2  # Same object (cached)

    def test_level_cached(self) -> None:
        """Test that level property is cached."""
        # Level property only exists in SemiDirectedPhyNetwork
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        level1 = net.level
        level2 = net.level
        assert level1 is level2  # Same object (cached)


class TestPropertyInvalidationOnCopy:
    """Test cases for property invalidation when copying."""

    def test_copy_does_not_copy_cache(self) -> None:
        """Test that copy creates new cache."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        leaves1 = net.leaves
        copy_net = net.copy()
        leaves2 = copy_net.leaves
        
        # Should have same values but different objects
        assert leaves1 == leaves2
        assert leaves1 is not leaves2  # Different objects

    def test_copy_preserves_properties(self) -> None:
        """Test that copy preserves property values."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        copy_net = net.copy()
        
        assert net.leaves == copy_net.leaves
        assert net.taxa == copy_net.taxa
        assert net.hybrid_nodes == copy_net.hybrid_nodes
        assert net.tree_nodes == copy_net.tree_nodes


class TestPropertyConsistency:
    """Test cases for property consistency."""

    def test_leaves_are_subset_of_nodes(self) -> None:
        """Test that leaves are a subset of all nodes."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        all_nodes = set(net)
        assert net.leaves.issubset(all_nodes)

    def test_internal_nodes_are_subset_of_nodes(self) -> None:
        """Test that internal nodes are a subset of all nodes."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        all_nodes = set(net)
        assert net.internal_nodes.issubset(all_nodes)

    def test_leaves_and_internal_nodes_disjoint(self) -> None:
        """Test that leaves and internal nodes are disjoint."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        leaves_set = net.leaves
        internal_set = net.internal_nodes
        assert leaves_set.isdisjoint(internal_set)

    def test_leaves_and_internal_nodes_cover_all_nodes(self) -> None:
        """Test that leaves and internal nodes cover all nodes."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        all_nodes = set(net)
        leaves_set = net.leaves
        internal_set = net.internal_nodes
        assert leaves_set | internal_set == all_nodes

    def test_hybrid_nodes_are_subset_of_internal_nodes(self) -> None:
        """Test that hybrid nodes are a subset of internal nodes."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        internal_set = net.internal_nodes
        hybrid_set = net.hybrid_nodes
        assert hybrid_set.issubset(internal_set)

    def test_tree_nodes_are_subset_of_internal_nodes(self) -> None:
        """Test that tree nodes are a subset of internal nodes."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        internal_set = net.internal_nodes
        tree_set = net.tree_nodes
        assert tree_set.issubset(internal_set)

    def test_hybrid_and_tree_nodes_disjoint(self) -> None:
        """Test that hybrid and tree nodes are disjoint."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        hybrid_set = net.hybrid_nodes
        tree_set = net.tree_nodes
        assert hybrid_set.isdisjoint(tree_set)

    def test_taxa_matches_leaves(self) -> None:
        """Test that taxa set matches number of leaves."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert len(net.taxa) == len(net.leaves)

    def test_hybrid_edges_are_directed(self) -> None:
        """Test that hybrid edges are all directed."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        for u, v, k in net.hybrid_edges:
            assert net._graph._directed.has_edge(u, v, key=k)
            assert not net._graph._undirected.has_edge(u, v)

    def test_tree_edges_are_undirected(self) -> None:
        """Test that tree edges are all undirected."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        for u, v, k in net.tree_edges:
            assert net._graph._undirected.has_edge(u, v, key=k)
            assert not net._graph._directed.has_edge(u, v)


class TestPropertyValues:
    """Test cases for property value correctness."""

    def test_leaves_have_no_outgoing_directed_edges(self) -> None:
        """Test that leaves have no outgoing directed edges."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        for leaf in net.leaves:
            assert net.outdegree(leaf) == 0

    def test_hybrid_nodes_have_indegree_at_least_two(self) -> None:
        """Test that hybrid nodes have indegree >= 2."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        for hybrid in net.hybrid_nodes:
            assert net._graph.indegree(hybrid) >= 2

    def test_hybrid_nodes_total_degree_equals_indegree_plus_one(self) -> None:
        """Test that hybrid nodes have total_degree = indegree + 1."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        for hybrid in net.hybrid_nodes:
            indegree = net._graph.indegree(hybrid)
            total_degree = net._graph.degree(hybrid)
            assert total_degree == indegree + 1

    def test_tree_nodes_have_indegree_zero(self) -> None:
        """Test that tree nodes have indegree 0."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        for tree_node in net.tree_nodes:
            assert net._graph.indegree(tree_node) == 0

    def test_tree_nodes_have_total_degree_at_least_three(self) -> None:
        """Test that tree nodes have total_degree >= 3."""
        # Node 6 needs degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        for tree_node in net.tree_nodes:
            total_degree = net._graph.degree(tree_node)
            assert total_degree >= 3

    def test_internal_nodes_have_degree_at_least_three(self) -> None:
        """Test that internal nodes have degree >= 3."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        for internal in net.internal_nodes:
            degree = net._graph.degree(internal)
            assert degree >= 3


class TestHasParallelEdges:
    """Test cases for has_parallel_edges property."""

    def test_no_parallel_edges(self) -> None:
        """Test network with no parallel edges."""
        from phylozoo.core.network.sdnetwork.classifications import has_parallel_edges
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        # Valid SemiDirectedPhyNetwork with one source component
        net = SemiDirectedPhyNetwork(
            directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.6}, {'u': 6, 'v': 4, 'gamma': 0.4}],
            undirected_edges=[
                (7, 5), (7, 6), (7, 8),  # Connect nodes 5, 6 to common source via node 7
                (4, 2), (5, 10), (6, 11)
            ],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (10, {'label': 'C'}), (11, {'label': 'D'})]
        )
        assert not has_parallel_edges(net)

    def test_with_parallel_undirected_edges(self) -> None:
        """Test network with parallel undirected edges."""
        from phylozoo.core.network.sdnetwork.classifications import has_parallel_edges
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        # Valid SemiDirectedPhyNetwork with parallel undirected edges
        # Node 7 has degree 4 (2 parallel to node 5, plus 2 others) - satisfies degree >= 3
        net = SemiDirectedPhyNetwork(
            directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.6}, {'u': 6, 'v': 4, 'gamma': 0.4}],
            undirected_edges=[
                (7, 5, 0), (7, 5, 1),  # Parallel undirected edges
                (7, 6), (7, 8),  # Additional edges to satisfy degree constraints
                (4, 2), (6, 9), (5, 10)
            ],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'})]
        )
        assert has_parallel_edges(net)

