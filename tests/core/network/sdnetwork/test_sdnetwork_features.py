"""
Tests for SemiDirectedPhyNetwork and MixedPhyNetwork features (cut_edges, cut_vertices).

This test suite covers the network-level cut_edges and cut_vertices functions,
including caching behavior.
"""

import pytest

from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork, MixedPhyNetwork
from phylozoo.core.network.sdnetwork.features import cut_edges, cut_vertices, root_locations


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


class TestRootLocations:
    """Test cases for root_locations function."""

    def test_simple_undirected_network(self) -> None:
        """Test root_locations on a simple network with only undirected edges."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        node_locs, undir_locs, dir_locs = root_locations(net)
        
        # Node 3 should be a valid root location (non-leaf in source component)
        assert 3 in node_locs
        assert len(node_locs) == 1
        
        # All undirected edges should be valid root locations
        assert len(undir_locs) == 3
        # Edges are normalized, so check for both orderings
        edge_tuples = {(u, v, k) for u, v, k in undir_locs}
        assert (1, 3, 0) in edge_tuples or (3, 1, 0) in edge_tuples
        assert (2, 3, 0) in edge_tuples or (3, 2, 0) in edge_tuples
        assert (3, 4, 0) in edge_tuples or (4, 3, 0) in edge_tuples
        
        # No directed edges exiting source component
        assert len(dir_locs) == 0

    def test_network_with_directed_edges_exiting(self) -> None:
        """Test root_locations with directed edges exiting the source component."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[
                (4, 1),  # Hybrid node 4: indegree 2, total degree 3 (2+1)
                (5, 7), (5, 8), (5, 9),  # Node 5 has degree 4
                (6, 7), (6, 10), (6, 11),  # Node 6 has degree 4
                (7, 12)  # Node 7 has degree 3
            ],
            nodes=[
                (1, {'label': 'A'}),
                (8, {'label': 'B'}), (9, {'label': 'C'}),
                (10, {'label': 'D'}), (11, {'label': 'E'}),
                (12, {'label': 'F'})
            ]
        )
        node_locs, undir_locs, dir_locs = root_locations(net)
        
        # Source component should contain nodes 5, 6, 7, 8, 9, 10, 11, 12 (connected via undirected edges)
        # Node 7, 8, 9, 10, 11, or 12 should be a valid root location
        assert len(node_locs) > 0
        
        # Undirected edges in source component
        assert len(undir_locs) >= 7  # At least (5, 7), (5, 8), (5, 9), (6, 7), (6, 10), (6, 11), (7, 12)
        
        # Directed edges exiting source component: (5, 4) and (6, 4)
        assert len(dir_locs) == 2
        edge_tuples = {(u, v, k) for u, v, k in dir_locs}
        assert (5, 4, 0) in edge_tuples
        assert (6, 4, 0) in edge_tuples

    def test_network_with_hybrid_node(self) -> None:
        """Test root_locations on a network with a hybrid node."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {'u': 5, 'v': 4, 'gamma': 0.6},
                {'u': 6, 'v': 4, 'gamma': 0.4}
            ],
            undirected_edges=[
                (4, 1),  # Hybrid node 4: indegree 2, total degree 3 (2+1)
                (5, 7), (5, 8), (5, 9),  # Node 5 has degree 4
                (6, 7), (6, 10), (6, 11),  # Node 6 has degree 4
                (7, 12)  # Node 7 has degree 3
            ],
            nodes=[
                (1, {'label': 'A'}),
                (8, {'label': 'B'}), (9, {'label': 'C'}),
                (10, {'label': 'D'}), (11, {'label': 'E'}),
                (12, {'label': 'F'})
            ]
        )
        node_locs, undir_locs, dir_locs = root_locations(net)
        
        # Source component should contain nodes 5, 6, 7, 8, 9, 10, 11, 12 (connected via undirected)
        # Node 7, 8, 9, 10, 11, or 12 should be a valid root location
        assert len(node_locs) > 0
        
        # Undirected edges in source component
        assert len(undir_locs) >= 7  # At least (5, 7), (5, 8), (5, 9), (6, 7), (6, 10), (6, 11), (7, 12)
        
        # Directed edges exiting source component: (5, 4) and (6, 4)
        assert len(dir_locs) == 2
        edge_tuples = {(u, v, k) for u, v, k in dir_locs}
        assert (5, 4, 0) in edge_tuples
        assert (6, 4, 0) in edge_tuples

    def test_only_node_locations(self) -> None:
        """Test root_locations when source component has nodes with outgoing directed edges."""
        # Similar structure to test_network_with_directed_edges_exiting but simpler
        # Source component: nodes 5, 6, 7, 8 (connected via undirected edges)
        # Nodes 5 and 6 have outgoing directed edges to hybrid node 4
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {'u': 5, 'v': 4, 'gamma': 0.6},
                {'u': 6, 'v': 4, 'gamma': 0.4}  # Hybrid node 4: indegree 2, total degree 3
            ],
            undirected_edges=[
                (4, 1),  # Hybrid to leaf
                (5, 7), (5, 8), (5, 9),  # Node 5 has degree 4
                (6, 7), (6, 10), (6, 11),  # Node 6 has degree 4
                (7, 12)  # Node 7 has degree 3
            ],
            nodes=[
                (1, {'label': 'A'}),
                (8, {'label': 'B'}), (9, {'label': 'C'}),
                (10, {'label': 'D'}), (11, {'label': 'E'}),
                (12, {'label': 'F'})
            ]
        )
        node_locs, undir_locs, dir_locs = root_locations(net)
        
        # Nodes 7, 8, 9, 10, 11, 12 should be valid root locations (in source component)
        assert len(node_locs) >= 1  # At least one internal node in source component
        
        # Undirected edges in source component
        assert len(undir_locs) >= 7  # At least (5, 7), (5, 8), (5, 9), (6, 7), (6, 10), (6, 11), (7, 12)
        
        # Directed edges exiting source component: (5, 4) and (6, 4)
        assert len(dir_locs) == 2
        edge_tuples = {(u, v, k) for u, v, k in dir_locs}
        assert (5, 4, 0) in edge_tuples
        assert (6, 4, 0) in edge_tuples

    def test_multiple_source_components_error(self) -> None:
        """Test that root_locations raises error for multiple source components."""
        # Create a network with disconnected components (should fail validation)
        # But if we bypass validation, we can test the error
        from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
        from phylozoo.core.network.sdnetwork.base import MixedPhyNetwork
        from phylozoo.utils.validation import no_validation
        
        graph = MixedMultiGraph(
            undirected_edges=[(1, 2), (3, 4)]
        )
        
        # This should have two source components
        with no_validation(classes=["MixedPhyNetwork"]):
            net = MixedPhyNetwork(
                undirected_edges=[(1, 2), (3, 4)],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (4, {'label': 'D'})]
            )
        
        # root_locations should raise PhyloZooNetworkStructureError for multiple source components
        from phylozoo.utils.exceptions import PhyloZooNetworkStructureError
        with pytest.raises(PhyloZooNetworkStructureError, match="exactly one source component"):
            root_locations(net)

    def test_empty_source_component_error(self) -> None:
        """Test that root_locations raises error for empty source component."""
        # This is a theoretical edge case - an empty network would fail validation
        # But we can test the error handling
        from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
        from phylozoo.utils.validation import no_validation
        
        # Create a network with no nodes (would fail validation normally)
        with no_validation(classes=["MixedPhyNetwork"]):
            from phylozoo.core.network.sdnetwork.base import MixedPhyNetwork
            net = MixedPhyNetwork()
        
        # root_locations should raise PhyloZooNetworkStructureError for no source components
        from phylozoo.utils.exceptions import PhyloZooNetworkStructureError
        with pytest.raises(PhyloZooNetworkStructureError, match="exactly one source component"):
            root_locations(net)

    def test_return_type(self) -> None:
        """Test that root_locations returns the correct tuple structure."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        result = root_locations(net)
        
        # Should return a tuple of three lists
        assert isinstance(result, tuple)
        assert len(result) == 3
        
        node_locs, undir_locs, dir_locs = result
        
        # All should be lists
        assert isinstance(node_locs, list)
        assert isinstance(undir_locs, list)
        assert isinstance(dir_locs, list)
        
        # Node locations should contain nodes
        for node in node_locs:
            assert isinstance(node, (int, str))  # Node ID type
        
        # Edge locations should be tuples of (u, v, key)
        for edge in undir_locs + dir_locs:
            assert isinstance(edge, tuple)
            assert len(edge) == 3
            u, v, key = edge
            assert isinstance(key, int)

    def test_parallel_undirected_edges(self) -> None:
        """Test root_locations with parallel undirected edges."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[
                (3, 1, 0), (3, 1, 1),  # Parallel edges (3, 1)
                (3, 2), (3, 4), (3, 5),  # Node 3 has degree 5
                (1, 6), (1, 7)  # Node 1 has degree 4 (2 parallel + 2 other)
            ],
            nodes=[
                (2, {'label': 'B'}), (4, {'label': 'C'}), (5, {'label': 'D'}),
                (6, {'label': 'E'}), (7, {'label': 'F'})
            ]
        )
        node_locs, undir_locs, dir_locs = root_locations(net)
        
        # Node 3 should be a valid root location
        assert 3 in node_locs
        
        # Both parallel edges should be in undir_locs
        assert len(undir_locs) == 7  # (3, 1, 0), (3, 1, 1), (3, 2, 0), (3, 4, 0), (3, 5, 0), (1, 6, 0), (1, 7, 0)
        
        # Check that parallel edges are included
        edge_tuples = {(u, v, k) for u, v, k in undir_locs}
        # Edges might be normalized, so check both orderings
        assert (1, 3, 0) in edge_tuples or (3, 1, 0) in edge_tuples
        assert (1, 3, 1) in edge_tuples or (3, 1, 1) in edge_tuples

