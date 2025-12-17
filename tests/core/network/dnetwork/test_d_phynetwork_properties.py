"""
Comprehensive tests for DirectedPhyNetwork cached properties.

This module tests cached property behavior including:
- Property caching (same object returned)
- Property invalidation on copy
- All cached properties
- Property consistency
"""

import warnings

import pytest

from phylozoo.core.network import DirectedPhyNetwork


class TestCachedPropertyBehavior:
    """Test cases for cached property behavior."""

    def test_leaves_cached(self) -> None:
        """Test that leaves property is cached."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        leaves1 = net.leaves
        leaves2 = net.leaves
        assert leaves1 is leaves2  # Same object (cached)

    def test_taxa_cached(self) -> None:
        """Test that taxa property is cached."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        taxa1 = net.taxa
        taxa2 = net.taxa
        assert taxa1 is taxa2  # Same object (cached)

    def test_internal_nodes_cached(self) -> None:
        """Test that internal_nodes property is cached."""
        net = DirectedPhyNetwork(
            edges=[(4, 3), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        internal1 = net.internal_nodes
        internal2 = net.internal_nodes
        assert internal1 is internal2  # Same object (cached)

    def test_root_node_cached(self) -> None:
        """Test that root_node property is cached."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        root1 = net.root_node
        root2 = net.root_node
        assert root1 == root2  # Same value

    def test_hybrid_nodes_cached(self) -> None:
        """Test that hybrid_nodes property is cached."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        hybrid1 = net.hybrid_nodes
        hybrid2 = net.hybrid_nodes
        assert hybrid1 is hybrid2  # Same object (cached)

    def test_tree_nodes_cached(self) -> None:
        """Test that tree_nodes property is cached."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        tree1 = net.tree_nodes
        tree2 = net.tree_nodes
        assert tree1 is tree2  # Same object (cached)

    def test_hybrid_edges_cached(self) -> None:
        """Test that hybrid_edges property is cached."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        hybrid_edges1 = net.hybrid_edges
        hybrid_edges2 = net.hybrid_edges
        assert hybrid_edges1 is hybrid_edges2  # Same object (cached)

    def test_tree_edges_cached(self) -> None:
        """Test that tree_edges property is cached."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        tree_edges1 = net.tree_edges
        tree_edges2 = net.tree_edges
        assert tree_edges1 is tree_edges2  # Same object (cached)

    def test_level_cached(self) -> None:
        """Test that level property is cached."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        level1 = net.level
        level2 = net.level
        assert level1 is level2  # Same object (cached)


class TestPropertyInvalidationOnCopy:
    """Test cases for property invalidation when copying."""

    def test_copy_does_not_copy_cache(self) -> None:
        """Test that copy() doesn't copy cached properties."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        # Access property to cache it
        _ = net.leaves
        
        net_copy = net.copy()
        # Cached properties should be recomputed (different objects)
        # Note: For sets/lists, we check they're equal but may be different objects
        assert net.leaves == net_copy.leaves
        # They should be equal but may be different objects (recomputed)

    def test_copy_properties_recomputed(self) -> None:
        """Test that properties are recomputed on copy."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        net_copy = net.copy()
        
        # All properties should be equal
        assert net.leaves == net_copy.leaves
        assert net.taxa == net_copy.taxa
        assert net.root_node == net_copy.root_node
        assert net.internal_nodes == net_copy.internal_nodes
        assert net.hybrid_nodes == net_copy.hybrid_nodes
        assert net.tree_nodes == net_copy.tree_nodes
        assert net.hybrid_edges == net_copy.hybrid_edges
        assert net.tree_edges == net_copy.tree_edges
        assert net.level == net_copy.level


class TestAllCachedProperties:
    """Test that all expected properties are cached."""

    def test_all_properties_are_cached(self) -> None:
        """Test that all topology properties use @cached_property."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        
        # Access all properties twice
        leaves1, leaves2 = net.leaves, net.leaves
        taxa1, taxa2 = net.taxa, net.taxa
        internal1, internal2 = net.internal_nodes, net.internal_nodes
        root1, root2 = net.root_node, net.root_node
        hybrid1, hybrid2 = net.hybrid_nodes, net.hybrid_nodes
        tree1, tree2 = net.tree_nodes, net.tree_nodes
        hybrid_edges1, hybrid_edges2 = net.hybrid_edges, net.hybrid_edges
        tree_edges1, tree_edges2 = net.tree_edges, net.tree_edges
        level1, level2 = net.level, net.level
        
        # Check they're cached (same object for mutable types, same value for immutable)
        assert leaves1 is leaves2
        assert taxa1 is taxa2
        assert internal1 is internal2
        assert root1 == root2
        assert hybrid1 is hybrid2
        assert tree1 is tree2
        assert hybrid_edges1 is hybrid_edges2
        assert tree_edges1 is tree_edges2
        assert level1 is level2


class TestPropertyConsistency:
    """Test consistency between different properties."""

    def test_leaves_plus_internal_plus_root_equals_all_nodes(self) -> None:
        """Test that leaves + internal_nodes + root = all nodes."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        all_nodes = set(net._graph.nodes)
        accounted = {net.root_node} | net.leaves | net.internal_nodes
        assert all_nodes == accounted

    def test_internal_nodes_equals_tree_plus_hybrid(self) -> None:
        """Test that internal_nodes = tree_nodes + hybrid_nodes."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        internal_set = net.internal_nodes
        tree_hybrid_set = net.tree_nodes | net.hybrid_nodes
        assert internal_set == tree_hybrid_set

    def test_taxa_matches_leaves(self) -> None:
        """Test that taxa count matches leaves count."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        assert len(net.taxa) == len(net.leaves)

    def test_hybrid_edges_count_matches_hybrid_indegrees(self) -> None:
        """Test that hybrid_edges count matches sum of hybrid in-degrees."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        hybrid_edge_count = len(net.hybrid_edges)
        hybrid_indegree_sum = sum(net.indegree(h) for h in net.hybrid_nodes)
        assert hybrid_edge_count == hybrid_indegree_sum

    def test_tree_edges_plus_hybrid_edges_equals_all_edges(self) -> None:
        """Test that tree_edges + hybrid_edges = all edges."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        all_edges = net.tree_edges | net.hybrid_edges
        expected_edges = {(u, v, k) for u, v, k in net._graph.edges(keys=True)}
        assert all_edges == expected_edges

    def test_root_not_in_internal_nodes(self) -> None:
        """Test that root is not in internal_nodes."""
        net = DirectedPhyNetwork(
            edges=[(4, 3), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        assert net.root_node not in net.internal_nodes

    def test_leaves_not_in_internal_nodes(self) -> None:
        """Test that leaves are not in internal_nodes."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        for leaf in net.leaves:
            assert leaf not in net.internal_nodes

    def test_hybrid_nodes_in_internal_nodes(self) -> None:
        """Test that hybrid_nodes are in internal_nodes."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        for hybrid in net.hybrid_nodes:
            assert hybrid in net.internal_nodes

    def test_tree_nodes_in_internal_nodes(self) -> None:
        """Test that tree_nodes are in internal_nodes."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        for tree_node in net.tree_nodes:
            assert tree_node in net.internal_nodes


class TestHasParallelEdges:
    """Test cases for has_parallel_edges property."""

    def test_no_parallel_edges(self) -> None:
        """Test network with no parallel edges."""
        from phylozoo.core.network.dnetwork.classifications import has_parallel_edges
        # Use a simple tree structure (no hybrid nodes, no parallel edges)
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        assert not has_parallel_edges(net)

    def test_with_parallel_edges(self) -> None:
        """Test network with parallel edges."""
        from phylozoo.core.network.dnetwork.classifications import has_parallel_edges
        # Valid network with parallel edges from root to internal node
        # Root 5 has 2 parallel edges to node 3 (out-degree 2)
        # Node 3 has in-degree 2, so must have out-degree 1
        # Node 3 connects to node 4 which splits to leaves
        net = DirectedPhyNetwork(
            edges=[
                (5, 3, 0), (5, 3, 1),  # Parallel edges from root to internal node
                (3, 4),  # Internal node 3 -> 4
                (4, 1), (4, 2)  # Internal node 4 splits to leaves
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        assert has_parallel_edges(net)

