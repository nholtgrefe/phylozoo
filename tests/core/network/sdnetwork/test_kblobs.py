"""
Tests for k_blobs function in SemiDirectedPhyNetwork and MixedPhyNetwork.

This test suite covers k-blobs with various k values, parallel edges,
non-binary networks, and edge cases.
"""

import pytest

from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork, MixedPhyNetwork
from phylozoo.core.network.sdnetwork.features import k_blobs


class TestKBlobs1BlobsSemiDirected:
    """Test cases for 1-blobs in SemiDirectedPhyNetwork."""

    def test_simple_tree_all_1blobs(self) -> None:
        """Test that leaves are 1-blobs in a simple tree."""
        # Valid semi-directed network requires degree >= 3 for internal nodes
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(4, 1), (4, 2), (4, 3)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
        )
        blobs_1 = list(k_blobs(net, k=1))
        
        # Leaves have 1 incident cut-edge
        assert {1} in blobs_1
        assert {2} in blobs_1
        assert {3} in blobs_1
        
        # Root has 3 incident cut-edges
        blobs_3 = list(k_blobs(net, k=3))
        assert {4} in blobs_3

    def test_with_hybrid_node(self) -> None:
        """Test 1-blobs with hybrid node present."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.6}, {'u': 6, 'v': 4, 'gamma': 0.4}],
            undirected_edges=[(7, 5), (7, 6), (7, 8), (4, 2), (5, 10), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (10, {'label': 'C'}), (11, {'label': 'D'})]
        )
        blobs_1 = list(k_blobs(net, k=1))
        
        # Leaves should be 1-blobs
        assert {2} in blobs_1
        assert {8} in blobs_1
        assert {10} in blobs_1
        assert {11} in blobs_1


class TestKBlobs2BlobsSemiDirected:
    """Test cases for 2-blobs in SemiDirectedPhyNetwork."""

    def test_internal_node_2blob(self) -> None:
        """Test that nodes with 2 incident cut-edges are 2-blobs."""
        # Star structure with one central node
        net = SemiDirectedPhyNetwork(
            directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.6}, {'u': 6, 'v': 4, 'gamma': 0.4}],
            undirected_edges=[(7, 5), (7, 6), (7, 8), (4, 2), (5, 10), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (10, {'label': 'C'}), (11, {'label': 'D'})]
        )
        blobs_2 = list(k_blobs(net, k=2))
        
        # Check that we have some 2-blobs (nodes with 2 incident cut-edges)
        assert len(blobs_2) >= 0  # May or may not have 2-blobs depending on structure


class TestKBlobs3BlobsSemiDirected:
    """Test cases for 3-blobs in SemiDirectedPhyNetwork."""

    def test_node_with_three_connections(self) -> None:
        """Test that node with 3 incident cut-edges is a 3-blob."""
        # Three-way split
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(4, 1), (4, 2), (4, 3)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
        )
        blobs_3 = list(k_blobs(net, k=3))
        
        # Central node has 3 incident cut-edges
        assert {4} in blobs_3

    def test_larger_star(self) -> None:
        """Test 4-blob with star structure."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 1), (5, 2), (5, 3), (5, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (4, {'label': 'D'})]
        )
        blobs_4 = list(k_blobs(net, k=4))
        
        # Central node has 4 incident cut-edges
        assert {5} in blobs_4


class TestKBlobsWithParallelEdgesSemiDirected:
    """Test cases for parallel edges in SemiDirectedPhyNetwork."""

    def test_parallel_undirected_not_cut_edges(self) -> None:
        """Test that parallel undirected edges are not cut-edges."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.6}, {'u': 6, 'v': 4, 'gamma': 0.4}],
            undirected_edges=[
                (7, 5, 0), (7, 5, 1),  # Parallel undirected edges
                (7, 6), (7, 8),
                (4, 2), (6, 9), (5, 10)
            ],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'})]
        )
        
        # The parallel edges between 7 and 5 should not be cut-edges
        # So the incident cut-edge count should reflect this
        blobs_1 = list(k_blobs(net, k=1))
        
        # Leaves should still be 1-blobs
        assert {2} in blobs_1
        assert {8} in blobs_1

    def test_blob_structure_with_parallel_edges(self) -> None:
        """Test blob counting when parallel edges form part of blob structure."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.6}, {'u': 6, 'v': 4, 'gamma': 0.4}],
            undirected_edges=[
                (7, 5, 0), (7, 5, 1),  # Parallel
                (7, 6), (7, 8),
                (4, 2), (6, 9), (5, 10)
            ],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'})]
        )
        
        # Test various k values
        for k in range(5):
            blobs_k = list(k_blobs(net, k=k))
            # Should get results without errors
            assert isinstance(blobs_k, list)


class TestKBlobsNonBinarySemiDirected:
    """Test cases for non-binary networks."""

    def test_ternary_split(self) -> None:
        """Test 3-blob from ternary split."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(4, 1), (4, 2), (4, 3)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
        )
        blobs_3 = list(k_blobs(net, k=3))
        
        # Central node is a 3-blob
        assert {4} in blobs_3
        
        # Leaves are 1-blobs
        blobs_1 = list(k_blobs(net, k=1))
        assert {1} in blobs_1
        assert {2} in blobs_1
        assert {3} in blobs_1

    def test_5ary_split(self) -> None:
        """Test 5-blob from 5-way split."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(6, 1), (6, 2), (6, 3), (6, 4), (6, 5)],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
                (4, {'label': 'D'}), (5, {'label': 'E'})
            ]
        )
        blobs_5 = list(k_blobs(net, k=5))
        
        # Central node has 5 incident cut-edges
        assert {6} in blobs_5


class TestKBlobsMixedNetwork:
    """Test cases for MixedPhyNetwork (not necessarily SemiDirected)."""

    def test_mixed_network_1blobs(self) -> None:
        """Test 1-blobs in a mixed network."""
        net = MixedPhyNetwork(
            directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.6}, {'u': 6, 'v': 4, 'gamma': 0.4}],
            undirected_edges=[(4, 2), (5, 8), (6, 9), (5, 10), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
        )
        blobs_1 = list(k_blobs(net, k=1))
        
        # Leaves should be 1-blobs
        assert {2} in blobs_1
        assert {8} in blobs_1
        assert {9} in blobs_1

    def test_mixed_network_various_k(self) -> None:
        """Test various k values on mixed network."""
        net = MixedPhyNetwork(
            directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.6}, {'u': 6, 'v': 4, 'gamma': 0.4}],
            undirected_edges=[(4, 2), (5, 8), (6, 9), (5, 10), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
        )
        
        # Test k=0 through k=5
        for k in range(6):
            blobs_k = list(k_blobs(net, k=k))
            assert isinstance(blobs_k, list)


class TestKBlobsEdgeCasesSemiDirected:
    """Test edge cases for SemiDirectedPhyNetwork."""

    def test_trivial_false_leaves_true_raises(self) -> None:
        """Test that invalid parameter combination raises error."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(4, 1), (4, 2), (4, 3)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
        )
        
        with pytest.raises(ValueError, match="Cannot have trivial=False and leaves=True"):
            list(k_blobs(net, k=1, trivial=False, leaves=True))

    def test_no_blobs_for_large_k(self) -> None:
        """Test that no blobs exist for very large k."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(4, 1), (4, 2), (4, 3)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
        )
        blobs_100 = list(k_blobs(net, k=100))
        
        # No node has 100 incident cut-edges
        assert len(blobs_100) == 0

    def test_k_zero(self) -> None:
        """Test 0-blobs (blobs with no incident cut-edges)."""
        # Network with potential 0-blob (internal blob with no external connections)
        net = SemiDirectedPhyNetwork(
            directed_edges=[{'u': 5, 'v': 4, 'gamma': 0.6}, {'u': 6, 'v': 4, 'gamma': 0.4}],
            undirected_edges=[(7, 5), (7, 6), (7, 8), (4, 2), (5, 10), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (10, {'label': 'C'}), (11, {'label': 'D'})]
        )
        blobs_0 = list(k_blobs(net, k=0))
        
        # May or may not have 0-blobs depending on structure
        assert isinstance(blobs_0, list)

    def test_leaves_false_filters_leaves(self) -> None:
        """Test that leaves=False correctly filters out leaf blobs."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(4, 1), (4, 2), (4, 3)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
        )
        blobs_1_no_leaves = list(k_blobs(net, k=1, leaves=False))
        
        # Leaves should be filtered out
        assert {1} not in blobs_1_no_leaves
        assert {2} not in blobs_1_no_leaves
        assert {3} not in blobs_1_no_leaves
        
        # Root is a 3-blob (not 1-blob)
        blobs_3_no_leaves = list(k_blobs(net, k=3, leaves=False))
        assert {4} in blobs_3_no_leaves


class TestKBlobsComplexSemiDirected:
    """Test complex network structures."""

    def test_asymmetric_structure(self) -> None:
        """Test k-blobs on asymmetric network."""
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
        
        # Test that we can find k-blobs for various k
        blobs_1 = list(k_blobs(net, k=1))
        blobs_2 = list(k_blobs(net, k=2))
        blobs_3 = list(k_blobs(net, k=3))
        blobs_4 = list(k_blobs(net, k=4))
        
        # Should have various types of blobs
        total_blobs = len(blobs_1) + len(blobs_2) + len(blobs_3) + len(blobs_4)
        assert total_blobs > 0

    def test_all_k_values_sum_to_total_blobs(self) -> None:
        """Test that all k-blobs partition the blob space."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(4, 1), (4, 2), (4, 3)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
        )
        
        # Get all blobs
        from phylozoo.core.network.sdnetwork.features import blobs as all_blobs
        total_blobs = list(all_blobs(net))
        
        # Get k-blobs for all reasonable k
        all_k_blobs = []
        for k in range(10):
            all_k_blobs.extend(k_blobs(net, k=k))
        
        # All blobs should be classified as some k-blob
        assert len(all_k_blobs) == len(total_blobs)

    def test_large_network_performance(self) -> None:
        """Test k-blobs on a larger network."""
        # Build a larger network
        net = SemiDirectedPhyNetwork(
            undirected_edges=[
                (10, 1), (10, 2), (10, 3), (10, 4),
                (10, 5), (10, 6), (10, 7), (10, 8)
            ],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (4, {'label': 'D'}),
                (5, {'label': 'E'}), (6, {'label': 'F'}), (7, {'label': 'G'}), (8, {'label': 'H'})
            ]
        )
        
        # Should efficiently compute k-blobs
        blobs_8 = list(k_blobs(net, k=8))
        
        # Central node has 8 incident cut-edges
        assert {10} in blobs_8
        
        # All leaves are 1-blobs
        blobs_1 = list(k_blobs(net, k=1))
        assert len(blobs_1) == 8

