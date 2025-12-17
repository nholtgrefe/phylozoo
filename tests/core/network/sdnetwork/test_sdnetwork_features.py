"""
Tests for SemiDirectedPhyNetwork and MixedPhyNetwork features (cut_edges, cut_vertices).

This test suite covers the network-level cut_edges and cut_vertices functions,
including caching behavior.
"""

import pytest

from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork, MixedPhyNetwork
from phylozoo.core.network.sdnetwork.features import cut_edges, cut_vertices


class TestCutEdgesSemiDirected:
    """Test cases for cut_edges on SemiDirectedPhyNetwork."""

    def test_simple_network(self) -> None:
        """Test cut_edges on a simple semi-directed network."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.6}, {'u': 6, 'v': 4, 'gamma': 0.4}],
            undirected_edges=[(7, 5), (7, 6), (7, 8), (4, 2), (5, 10), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (10, {'label': 'C'}), (11, {'label': 'D'})]
        )
        edges = cut_edges(net)
        # Should have some bridges
        assert len(edges) > 0
        # All edges should be 3-tuples
        for edge in edges:
            assert len(edge) == 3

    def test_with_parallel_undirected(self) -> None:
        """Test cut_edges with parallel undirected edges."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.6}, {'u': 6, 'v': 4, 'gamma': 0.4}],
            undirected_edges=[
                (7, 5, 0), (7, 5, 1),  # Parallel
                (7, 6), (7, 8),
                (4, 2), (6, 9), (5, 10)
            ],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'})]
        )
        edges = cut_edges(net)
        # Parallel edges should not be bridges
        parallel_count = sum(1 for u, v, k in edges if (u, v) == (7, 5) or (u, v) == (5, 7))
        assert parallel_count == 0  # No parallel edges should be bridges

    def test_larger_network(self) -> None:
        """Test cut_edges on a larger semi-directed network."""
        # Use same basic structure but ensure all internal nodes have degree >= 3
        net = SemiDirectedPhyNetwork(
            directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.6}, {'u': 6, 'v': 4, 'gamma': 0.4}],
            undirected_edges=[
                (7, 5), (7, 6), (7, 8), (7, 9),  # Node 7 has degree 4
                (4, 2), (5, 10), (6, 11)
            ],
            nodes=[
                (2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}),
                (10, {'label': 'D'}), (11, {'label': 'E'})
            ]
        )
        edges = cut_edges(net)
        assert len(edges) > 0


class TestCutVerticesSemiDirected:
    """Test cases for cut_vertices on SemiDirectedPhyNetwork."""

    def test_simple_network(self) -> None:
        """Test cut_vertices on a simple semi-directed network."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.6}, {'u': 6, 'v': 4, 'gamma': 0.4}],
            undirected_edges=[(7, 5), (7, 6), (7, 8), (4, 2), (5, 10), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (10, {'label': 'C'}), (11, {'label': 'D'})]
        )
        vertices = cut_vertices(net)
        # Should have some articulation points
        assert isinstance(vertices, set)

    def test_star_structure(self) -> None:
        """Test cut_vertices on star-like structure."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.6}, {'u': 6, 'v': 4, 'gamma': 0.4}],
            undirected_edges=[
                (7, 5), (7, 6), (7, 8), (7, 9),  # Star center at 7
                (4, 2), (5, 10), (6, 11)
            ],
            nodes=[
                (2, {'label': 'A'}), (8, {'label': 'B'}),
                (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})
            ]
        )
        vertices = cut_vertices(net)
        # Node 7 should be a cut vertex (star center)
        assert 7 in vertices


class TestCutEdgesMixed:
    """Test cases for cut_edges on MixedPhyNetwork."""

    def test_mixed_network(self) -> None:
        """Test cut_edges on a mixed network."""
        net = MixedPhyNetwork(
            directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.6}, {'u': 6, 'v': 4, 'gamma': 0.4}],
            undirected_edges=[(4, 2), (5, 8), (6, 9), (5, 10), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
        )
        edges = cut_edges(net)
        # Should return 3-tuples
        for edge in edges:
            assert len(edge) == 3


class TestCachingSemiDirected:
    """Test cases for caching behavior on semi-directed networks."""

    def test_cut_edges_cache(self) -> None:
        """Test that cut_edges properly caches results."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.6}, {'u': 6, 'v': 4, 'gamma': 0.4}],
            undirected_edges=[(7, 5), (7, 6), (7, 8), (4, 2), (5, 10), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (10, {'label': 'C'}), (11, {'label': 'D'})]
        )
        
        # Clear cache
        cut_edges.cache_clear()
        
        # First call
        result1 = cut_edges(net)
        cache_info1 = cut_edges.cache_info()
        assert cache_info1.misses == 1
        assert cache_info1.hits == 0
        
        # Second call - should hit cache
        result2 = cut_edges(net)
        cache_info2 = cut_edges.cache_info()
        assert cache_info2.hits == 1
        
        # Verify same object returned
        assert result1 is result2

    def test_cut_vertices_cache(self) -> None:
        """Test that cut_vertices properly caches results."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.6}, {'u': 6, 'v': 4, 'gamma': 0.4}],
            undirected_edges=[(7, 5), (7, 6), (7, 8), (4, 2), (5, 10), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (10, {'label': 'C'}), (11, {'label': 'D'})]
        )
        
        # Clear cache
        cut_vertices.cache_clear()
        
        # First call
        result1 = cut_vertices(net)
        cache_info1 = cut_vertices.cache_info()
        assert cache_info1.misses == 1
        
        # Second call - should hit cache
        result2 = cut_vertices(net)
        cache_info2 = cut_vertices.cache_info()
        assert cache_info2.hits == 1
        
        # Verify same object
        assert result1 is result2

    def test_cache_multiple_networks(self) -> None:
        """Test caching with multiple network instances."""
        net1 = SemiDirectedPhyNetwork(
            directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.6}, {'u': 6, 'v': 4, 'gamma': 0.4}],
            undirected_edges=[(7, 5), (7, 6), (7, 8), (4, 2), (5, 10), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (10, {'label': 'C'}), (11, {'label': 'D'})]
        )
        net2 = SemiDirectedPhyNetwork(
            directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.5}, {'u': 6, 'v': 4, 'gamma': 0.5}],
            undirected_edges=[(7, 5), (7, 6), (7, 8), (4, 2), (5, 10), (6, 11)],  # Same structure, different gammas
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (10, {'label': 'C'}), (11, {'label': 'D'})]
        )
        
        # Clear cache
        cut_edges.cache_clear()
        
        # Call on both
        edges1 = cut_edges(net1)
        edges2 = cut_edges(net2)
        
        # Both in cache
        cache_info = cut_edges.cache_info()
        assert cache_info.currsize == 2
        assert cache_info.misses == 2
        
        # Call again
        edges1_again = cut_edges(net1)
        edges2_again = cut_edges(net2)
        
        # Should have hits
        cache_info = cut_edges.cache_info()
        assert cache_info.hits == 2
        
        # Verify cached correctly
        assert edges1 is edges1_again
        assert edges2 is edges2_again

    def test_cache_lru_behavior(self) -> None:
        """Test that cache has proper LRU behavior."""
        # Cache has maxsize=128, so this shouldn't overflow
        networks = []
        for i in range(10):
            net = SemiDirectedPhyNetwork(
                directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.6}, {'u': 6, 'v': 4, 'gamma': 0.4}],
                undirected_edges=[(7, 5), (7, 6), (7, 8), (4, 2), (5, 10), (6, 11)],
                nodes=[(2, {'label': str(i)}), (8, {'label': 'B'}), (10, {'label': 'C'}), (11, {'label': 'D'})]
            )
            networks.append(net)
        
        # Clear and populate cache
        cut_edges.cache_clear()
        for net in networks:
            cut_edges(net)
        
        # All should be in cache
        cache_info = cut_edges.cache_info()
        assert cache_info.currsize == 10
        assert cache_info.misses == 10
        
        # Access first network again
        cut_edges(networks[0])
        cache_info = cut_edges.cache_info()
        assert cache_info.hits == 1

