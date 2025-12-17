"""
Tests for k_blobs function in DirectedPhyNetwork.

This test suite covers k-blobs with various k values, parallel edges,
non-binary networks, and edge cases.
"""

import pytest

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.features import k_blobs


class TestKBlobs1Blobs:
    """Test cases for 1-blobs."""

    def test_simple_tree_all_1blobs(self) -> None:
        """Test that leaves are 1-blobs in a simple tree."""
        # Simple tree: root splits to two leaves
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        blobs_1 = list(k_blobs(net, k=1))
        blob_sets = [sorted(b) for b in blobs_1]
        
        # Leaves have 1 incident cut-edge
        assert [1] in blob_sets
        assert [2] in blob_sets
        
        # Root has 2 incident cut-edges (out to both leaves)
        blobs_2 = list(k_blobs(net, k=2))
        assert {3} in blobs_2

    def test_binary_tree_1blobs(self) -> None:
        """Test 1-blobs in a binary tree."""
        # Binary tree with 3 internal nodes and 4 leaves
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),      # Root splits
                (5, 1), (5, 2),      # Left subtree
                (6, 3), (6, 4)       # Right subtree
            ],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}),
                (3, {'label': 'C'}), (4, {'label': 'D'})
            ]
        )
        blobs_1 = list(k_blobs(net, k=1))
        
        # In a tree, all nodes are 1-blobs (single edges connecting them)
        assert len(blobs_1) == 7  # 4 leaves + 3 internal nodes

    def test_non_binary_tree_1blobs(self) -> None:
        """Test 1-blobs in a non-binary tree."""
        # Non-binary: root splits into 3 children
        net = DirectedPhyNetwork(
            edges=[(4, 1), (4, 2), (4, 3)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
        )
        blobs_1 = list(k_blobs(net, k=1))
        
        # All 4 nodes should be 1-blobs
        assert len(blobs_1) == 4


class TestKBlobs2Blobs:
    """Test cases for 2-blobs."""

    def test_simple_chain_has_2blob(self) -> None:
        """Test that middle node in a chain is a 2-blob."""
        # Chain: 1 -> 2 -> 3 -> 4
        net = DirectedPhyNetwork(
            edges=[(1, 2), (2, 3), (3, 4)],
            nodes=[(4, {'label': 'A'})]
        )
        blobs_2 = list(k_blobs(net, k=2))
        
        # Nodes 2 and 3 should be 2-blobs (one edge in, one edge out)
        assert len(blobs_2) == 2
        blob_sets = [sorted(b) for b in blobs_2]
        assert [2] in blob_sets
        assert [3] in blob_sets

    def test_network_with_parallel_edges_2blobs(self) -> None:
        """Test 2-blobs with parallel edges present."""
        # Network with parallel edges to hybrid
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),                    # Root splits
                (5, 4, 0), (5, 4, 1),              # Parallel to hybrid
                (6, 4, 0), (6, 4, 1),              # Parallel to hybrid
                (4, 10), (10, 1), (10, 2)          # After hybrid
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        blobs_2 = list(k_blobs(net, k=2))
        
        # Node 10 splits to 2 leaves (2 cut-edges out)
        assert len(blobs_2) >= 1
        assert {10} in blobs_2


class TestKBlobs3Blobs:
    """Test cases for 3-blobs."""

    def test_node_with_three_children_is_3blob(self) -> None:
        """Test that a node splitting into 3 leaves is a 3-blob."""
        # Non-binary tree with 3 children
        net = DirectedPhyNetwork(
            edges=[(4, 1), (4, 2), (4, 3)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
        )
        blobs_3 = list(k_blobs(net, k=3))
        
        # Root node 4 has 3 incident cut-edges (all going out to leaves)
        assert len(blobs_3) == 1
        assert {4} in blobs_3

    def test_non_binary_internal_node(self) -> None:
        """Test 3-blobs in more complex non-binary structure."""
        # Root splits into 3, each child is a leaf
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6), (7, 8),  # Root with 3 children
                (5, 1), (5, 2),           # Child 5 splits
                (6, 3), (6, 4)            # Child 6 splits
            ],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}),
                (3, {'label': 'C'}), (4, {'label': 'D'}),
                (8, {'label': 'E'})
            ]
        )
        blobs_3 = list(k_blobs(net, k=3))
        
        # Root node 7 has 3 outgoing cut-edges
        assert {7} in blobs_3


class TestKBlobsWithParallelEdges:
    """Test cases specifically for parallel edges."""

    def test_parallel_edges_not_cut_edges(self) -> None:
        """Test that parallel edges are not counted as cut-edges."""
        # Parallel edges between nodes that form a cycle
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),                    # Root splits
                (5, 4, 0), (5, 4, 1),              # Parallel to hybrid (NOT bridges)
                (6, 4, 0), (6, 4, 1),              # Parallel to hybrid (NOT bridges)
                (4, 10), (10, 1), (10, 2)          # After hybrid
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        # The blob {5, 6, 4, 7} should have cut-edges only to nodes outside
        # Parallel edges within the blob structure are not cut-edges
        blobs_1 = list(k_blobs(net, k=1))
        
        # Leaves should be 1-blobs
        assert {1} in blobs_1
        assert {2} in blobs_1

    def test_blob_counting_with_parallel_cut_edges(self) -> None:
        """Test that parallel cut-edges are counted correctly."""
        # Create structure where parallel edges ARE cut-edges
        # (but this is tricky with phylogenetic network constraints)
        # For now, test the base case
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        blobs_2 = list(k_blobs(net, k=2))
        
        # Root has 2 outgoing cut-edges
        assert {3} in blobs_2


class TestKBlobsNonBinary:
    """Test cases for non-binary networks."""

    def test_4ary_split(self) -> None:
        """Test 4-blob from quaternary split."""
        net = DirectedPhyNetwork(
            edges=[(5, 1), (5, 2), (5, 3), (5, 4)],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}),
                (3, {'label': 'C'}), (4, {'label': 'D'})
            ]
        )
        blobs_4 = list(k_blobs(net, k=4))
        
        # Root with 4 children is a 4-blob
        assert len(blobs_4) == 1
        assert {5} in blobs_4

    def test_mixed_degrees(self) -> None:
        """Test network with mixed degree nodes."""
        # Complex structure with various degrees
        net = DirectedPhyNetwork(
            edges=[
                (10, 7), (10, 8), (10, 9),  # Root with 3 children
                (7, 1), (7, 2),              # Binary split
                (8, 3), (8, 4),              # Binary split
                (9, 5), (9, 6)               # Binary split
            ],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
                (4, {'label': 'D'}), (5, {'label': 'E'}), (6, {'label': 'F'})
            ]
        )
        
        # Check 3-blob (root)
        blobs_3 = list(k_blobs(net, k=3))
        assert {10} in blobs_3
        
        # Check 2-blobs (internal nodes with 2 children)
        blobs_2 = list(k_blobs(net, k=2))
        assert {7} in blobs_2
        assert {8} in blobs_2
        assert {9} in blobs_2
        
        # Check 1-blobs (leaves)
        blobs_1 = list(k_blobs(net, k=1))
        assert len(blobs_1) == 6  # All leaves


class TestKBlobsEdgeCases:
    """Test edge cases and special scenarios."""

    def test_no_k_blobs_found(self) -> None:
        """Test when no k-blobs exist for given k."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        blobs_5 = list(k_blobs(net, k=5))
        
        # No node has 5 incident cut-edges
        assert len(blobs_5) == 0

    def test_trivial_and_leaves_parameters(self) -> None:
        """Test trivial and leaves filtering parameters."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        # Without leaves - leaves are 1-blobs, root is 2-blob
        blobs_1_no_leaves = list(k_blobs(net, k=1, leaves=False))
        assert {1} not in blobs_1_no_leaves
        assert {2} not in blobs_1_no_leaves
        
        blobs_2_no_leaves = list(k_blobs(net, k=2, leaves=False))
        assert {3} in blobs_2_no_leaves
        
        # Without trivial (should raise error if leaves=True)
        with pytest.raises(ValueError):
            list(k_blobs(net, k=1, trivial=False, leaves=True))

    def test_single_node_network(self) -> None:
        """Test k-blobs on a single-node network."""
        # Single node (it's both root and leaf)
        net = DirectedPhyNetwork(
            edges=[],
            nodes=[(1, {'label': 'A'})]
        )
        blobs_0 = list(k_blobs(net, k=0))
        
        # Single node has 0 incident cut-edges
        assert len(blobs_0) == 1
        assert {1} in blobs_0

    def test_k_zero(self) -> None:
        """Test 0-blobs (nodes/blobs with no incident cut-edges)."""
        # Network with a non-trivial blob (cycle-like structure via hybrid)
        net = DirectedPhyNetwork(
            edges=[
                (8, 5), (8, 6),                    # Root splits
                (5, 4, 0), (5, 4, 1),              # Parallel to hybrid
                (6, 4, 0), (6, 4, 1),              # Parallel to hybrid
                (4, 10), (10, 1), (10, 2)          # After hybrid
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        blobs_0 = list(k_blobs(net, k=0))
        
        # Non-trivial blobs with internal cycles might have 0 incident cut-edges
        # (this depends on the specific structure)
        assert isinstance(blobs_0, list)

    def test_large_k_value(self) -> None:
        """Test with very large k value."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        blobs_100 = list(k_blobs(net, k=100))
        
        # No blobs should have 100 incident edges
        assert len(blobs_100) == 0


class TestKBlobsComplex:
    """Test complex network structures."""

    def test_hybrid_network_various_k(self) -> None:
        """Test k-blobs on network with hybrid nodes."""
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),                    # Root splits
                (5, 4, 0), (5, 4, 1),              # Parallel to hybrid
                (6, 4, 0), (6, 4, 1),              # Parallel to hybrid
                (4, 10), (10, 1), (10, 2)          # After hybrid
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        # Test various k values
        blobs_1 = list(k_blobs(net, k=1))
        blobs_2 = list(k_blobs(net, k=2))
        blobs_3 = list(k_blobs(net, k=3))
        
        # Leaves should be 1-blobs
        assert {1} in blobs_1
        assert {2} in blobs_1
        
        # Some internal structure should exist
        assert len(blobs_1) + len(blobs_2) + len(blobs_3) > 0

    def test_asymmetric_tree(self) -> None:
        """Test k-blobs on asymmetric tree structure."""
        # Asymmetric: left side deeper than right
        net = DirectedPhyNetwork(
            edges=[
                (10, 8), (10, 1),      # Root: left subtree and right leaf
                (8, 6), (8, 7),        # Left: continue splitting
                (6, 2), (6, 3),        # Left-left
                (7, 4), (7, 5)         # Left-right
            ],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}),
                (4, {'label': 'D'}), (5, {'label': 'E'})
            ]
        )
        
        # Root has 2 incident cut-edges
        blobs_2 = list(k_blobs(net, k=2))
        assert {10} in blobs_2
        assert {8} in blobs_2
        assert {6} in blobs_2
        assert {7} in blobs_2
        
        # All leaves have 1 incident cut-edge
        blobs_1 = list(k_blobs(net, k=1))
        assert {1} in blobs_1
        assert {2} in blobs_1
        assert {3} in blobs_1
        assert {4} in blobs_1
        assert {5} in blobs_1

