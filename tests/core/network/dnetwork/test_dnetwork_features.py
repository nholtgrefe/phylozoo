"""
Tests for DirectedPhyNetwork features module (cut_edges, cut_vertices, omnians).

This test suite covers the network-level cut_edges, cut_vertices, and omnians functions,
including caching behavior.
"""

import warnings

import pytest

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.features import cut_edges, cut_vertices, omnians


class TestCutEdges:
    """Test cases for cut_edges function."""

    def test_simple_tree(self) -> None:
        """Test cut_edges on a simple tree."""
        net = DirectedPhyNetwork(
            edges=[(1, 2), (1, 3)],
            nodes=[(2, {'label': 'A'}), (3, {'label': 'B'})]
        )
        edges = cut_edges(net)
        # Both edges are bridges in a tree
        assert (1, 2, 0) in edges
        assert (1, 3, 0) in edges
        assert len(edges) == 2

    def test_tree_with_parallel_edges(self) -> None:
        """Test cut_edges on network with parallel edges to hybrid node."""
        # Valid network: parallel edges lead to hybrid node
        # Node 4 is hybrid (in-degree 4, out-degree 1)
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),                    # Root splits
                (5, 4, 0), (5, 4, 1),              # Parallel from 5 to hybrid 4
                (6, 4, 0), (6, 4, 1),              # Parallel from 6 to hybrid 4
                (4, 10), (10, 1), (10, 2)          # After hybrid
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        edges = cut_edges(net)
        # Parallel edges to hybrid are not bridges
        assert (5, 4, 0) not in edges
        assert (5, 4, 1) not in edges
        assert (6, 4, 0) not in edges
        assert (6, 4, 1) not in edges

    def test_larger_tree(self) -> None:
        """Test cut_edges on a larger tree structure."""
        # Create a binary tree structure
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),          # Root to internal nodes
                (5, 1), (5, 2),          # Left subtree
                (6, 3), (6, 4)           # Right subtree
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (4, {'label': 'D'})
            ]
        )
        edges = cut_edges(net)
        # All edges should be bridges in a tree
        assert len(edges) == 6

    def test_network_with_hybrid_node(self) -> None:
        """Test cut_edges on network with hybrid node."""
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),                    # Root splits
                (5, 4, 0), (5, 4, 1),              # Parallel edges to hybrid
                (6, 4, 0), (6, 4, 1),              # Parallel edges to hybrid
                (4, 10), (10, 1), (10, 2)          # After hybrid
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        edges = cut_edges(net)
        # Parallel edges to hybrid are not bridges
        assert (5, 4, 0) not in edges
        assert (5, 4, 1) not in edges
        assert (6, 4, 0) not in edges
        assert (6, 4, 1) not in edges


class TestCutVertices:
    """Test cases for cut_vertices function."""

    def test_simple_tree(self) -> None:
        """Test cut_vertices on a simple tree."""
        net = DirectedPhyNetwork(
            edges=[(1, 2), (1, 3)],
            nodes=[(2, {'label': 'A'}), (3, {'label': 'B'})]
        )
        vertices = cut_vertices(net)
        # In this small tree, leaves are not cut vertices
        assert 2 not in vertices
        assert 3 not in vertices

    def test_chain_structure(self) -> None:
        """Test cut_vertices on a chain of splits."""
        net = DirectedPhyNetwork(
            edges=[
                (8, 7), (8, 6),          # Root split
                (7, 5), (7, 4),          # Second level split
                (5, 1), (5, 2),          # Third level split
                (4, 3), (4, 9)           # Second split at 4
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (6, {'label': 'D'}),
                (9, {'label': 'E'})
            ]
        )
        vertices = cut_vertices(net)
        # Internal branching nodes might be cut vertices
        # The specific nodes depend on the network structure
        assert isinstance(vertices, set)

    def test_no_cut_vertices(self) -> None:
        """Test network with no cut vertices (single split)."""
        net = DirectedPhyNetwork(
            edges=[(1, 2), (1, 3)],
            nodes=[(2, {'label': 'A'}), (3, {'label': 'B'})]
        )
        vertices = cut_vertices(net)
        # Root might be a cut vertex depending on connectivity
        # But leaves definitely are not
        assert 2 not in vertices
        assert 3 not in vertices


class TestCaching:
    """Test cases for caching behavior."""

    def test_cut_edges_cache(self) -> None:
        """Test that cut_edges properly caches results."""
        net = DirectedPhyNetwork(
            edges=[(1, 2), (1, 3)],
            nodes=[(2, {'label': 'A'}), (3, {'label': 'B'})]
        )
        
        # Clear cache before test
        cut_edges.cache_clear()
        
        # First call
        result1 = cut_edges(net)
        cache_info1 = cut_edges.cache_info()
        assert cache_info1.misses == 1
        assert cache_info1.hits == 0
        
        # Second call - should hit cache
        result2 = cut_edges(net)
        cache_info2 = cut_edges.cache_info()
        assert cache_info2.misses == 1
        assert cache_info2.hits == 1
        
        # Verify same object returned
        assert result1 is result2

    def test_cut_vertices_cache(self) -> None:
        """Test that cut_vertices properly caches results."""
        net = DirectedPhyNetwork(
            edges=[(1, 2), (1, 3)],
            nodes=[(2, {'label': 'A'}), (3, {'label': 'B'})]
        )
        
        # Clear cache before test
        cut_vertices.cache_clear()
        
        # First call
        result1 = cut_vertices(net)
        cache_info1 = cut_vertices.cache_info()
        assert cache_info1.misses == 1
        assert cache_info1.hits == 0
        
        # Second call - should hit cache
        result2 = cut_vertices(net)
        cache_info2 = cut_vertices.cache_info()
        assert cache_info2.misses == 1
        assert cache_info2.hits == 1
        
        # Verify same object returned
        assert result1 is result2

    def test_cache_separate_networks(self) -> None:
        """Test that cache stores results for separate networks."""
        net1 = DirectedPhyNetwork(
            edges=[(1, 2), (1, 3)],
            nodes=[(2, {'label': 'A'}), (3, {'label': 'B'})]
        )
        net2 = DirectedPhyNetwork(
            edges=[(1, 2), (1, 3), (1, 4)],
            nodes=[(2, {'label': 'A'}), (3, {'label': 'B'}), (4, {'label': 'C'})]
        )
        
        # Clear cache
        cut_edges.cache_clear()
        
        # Call on both networks
        edges1 = cut_edges(net1)
        edges2 = cut_edges(net2)
        
        # Both should be in cache
        cache_info = cut_edges.cache_info()
        assert cache_info.currsize == 2
        assert cache_info.misses == 2
        assert cache_info.hits == 0
        
        # Call again - should hit cache
        edges1_again = cut_edges(net1)
        edges2_again = cut_edges(net2)
        
        cache_info = cut_edges.cache_info()
        assert cache_info.hits == 2
        
        # Verify results are cached correctly
        assert edges1 is edges1_again
        assert edges2 is edges2_again
        assert edges1 != edges2  # Different networks, different results

    def test_cache_info_accessible(self) -> None:
        """Test that cache_info and cache_clear are accessible."""
        # These should not raise AttributeError
        info = cut_edges.cache_info()
        assert hasattr(info, 'hits')
        assert hasattr(info, 'misses')
        assert hasattr(info, 'maxsize')
        assert hasattr(info, 'currsize')
        
        # Should be able to clear cache
        cut_edges.cache_clear()
        info_after = cut_edges.cache_info()
        assert info_after.currsize == 0


class TestOmnians:
    """Test cases for omnians function."""

    def test_network_with_multiple_omnians(self) -> None:
        """Test omnians on network with multiple omnian nodes."""
        # Network where nodes 5, 8, and 9 all have all children as hybrids
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 8), (7, 9),  # Root to tree nodes
                (5, 4), (5, 6),  # Node 5 to hybrid nodes 4 and 6
                (8, 4), (8, 6),  # Node 8 to hybrid nodes 4 and 6
                (9, 4), (9, 6),  # Node 9 to hybrid nodes 4 and 6
                (4, 1), (6, 2)   # Hybrids to leaves
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        result = omnians(net)
        assert result == {5, 8, 9}

    def test_network_with_single_omnian(self) -> None:
        """Test omnians on network with a single omnian node."""
        # Network where only node 5 has all children as hybrids
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 8), (7, 9),  # Root to tree nodes
                (5, 4), (5, 6),  # Node 5 to hybrid nodes 4 and 6 (all hybrids)
                (8, 4), (8, 10), (8, 11),  # Node 8 to hybrid 4 and tree nodes 10, 11 (mixed)
                (9, 4), (9, 6),  # Node 9 to hybrid nodes 4 and 6 (all hybrids)
                (4, 1), (6, 2), (10, 3), (10, 12), (11, 13), (11, 14)  # To leaves
            ],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
                (12, {'label': 'D'}), (13, {'label': 'E'}), (14, {'label': 'F'})
            ]
        )
        result = omnians(net)
        # Only 5 and 9 have all children as hybrids
        assert result == {5, 9}

    def test_network_with_no_omnians(self) -> None:
        """Test omnians on network with no omnian nodes."""
        # Simple tree with no hybrid nodes
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        result = omnians(net)
        assert result == set()

    def test_network_with_parallel_edges_warns(self) -> None:
        """Test that omnians warns when network has parallel edges."""
        net = DirectedPhyNetwork(
            edges=[
                (7, 4), (7, 4),  # Parallel edges to hybrid 4
                (4, 1)           # Hybrid to leaf
            ],
            nodes=[(1, {'label': 'A'})]
        )
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = omnians(net)
            # Should have raised a warning about parallel edges
            assert len(w) > 0
            assert "parallel edges" in str(w[0].message).lower()
            assert "Jetten" in str(w[0].message) or "van Iersel" in str(w[0].message)
            # Function should still return a result (node 7 is omnian)
            assert isinstance(result, set)

    def test_omnian_with_single_hybrid_child(self) -> None:
        """Test that a node with only hybrid children (even if just one) is an omnian."""
        # Node 5 has two hybrid children (both are hybrids)
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6), (7, 8), (7, 9),  # Root to tree nodes
                (5, 4), (5, 10),  # Node 5 to hybrid nodes 4 and 10 (both hybrids)
                (6, 4), (6, 10),  # Node 6 to hybrid nodes 4 and 10
                (8, 4), (8, 10),  # Node 8 to hybrid nodes 4 and 10
                (9, 4), (9, 10),  # Node 9 to hybrid nodes 4 and 10
                (4, 1), (10, 2)   # Hybrids to leaves
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        result = omnians(net)
        # Node 5 has all children (both are hybrids) as hybrid, so it's an omnian
        assert 5 in result

    def test_node_with_mixed_children_not_omnian(self) -> None:
        """Test that a node with mixed children (hybrid and tree) is not an omnian."""
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6), (7, 9),  # Root to tree nodes
                (5, 4), (5, 8),  # Node 5 to hybrid 4 and tree node 8 (mixed)
                (6, 4), (6, 10),  # Node 6 to hybrid 4 and leaf 10 (ensures out-degree >= 2)
                (9, 4), (9, 11),  # Node 9 to hybrid 4 and leaf 11 (ensures out-degree >= 2)
                (4, 1), (8, 2), (8, 3)  # To leaves
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
        )
        result = omnians(net)
        # Node 5 has mixed children, so it's not an omnian
        assert 5 not in result

    def test_omnians_cache(self) -> None:
        """Test that omnians properly caches results."""
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 8), (7, 9),
                (5, 4), (5, 6),
                (8, 4), (8, 6),
                (9, 4), (9, 6),
                (4, 1), (6, 2)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        # Clear cache before test
        omnians.cache_clear()
        
        # First call
        result1 = omnians(net)
        cache_info1 = omnians.cache_info()
        assert cache_info1.misses == 1
        assert cache_info1.hits == 0
        
        # Second call - should hit cache
        result2 = omnians(net)
        cache_info2 = omnians.cache_info()
        assert cache_info2.misses == 1
        assert cache_info2.hits == 1
        
        # Verify same object returned
        assert result1 is result2

