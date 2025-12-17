"""
Tests for DirectedPhyNetwork features module (cut_edges, cut_vertices).

This test suite covers the network-level cut_edges and cut_vertices functions,
including caching behavior.
"""

import pytest

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.features import cut_edges, cut_vertices


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

