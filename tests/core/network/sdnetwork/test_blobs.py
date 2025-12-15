"""
Tests for blob computation in SemiDirectedPhyNetwork and MixedPhyNetwork.
"""

import warnings

import pytest

from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
from phylozoo.core.network.sdnetwork.features import blobs


class TestBlobs:
    """Test blob computation for SemiDirectedPhyNetwork."""

    def test_empty_network(self) -> None:
        """Test that empty network has no blobs."""
        with pytest.warns(UserWarning, match="Empty network.*no nodes"):
            net = SemiDirectedPhyNetwork(undirected_edges=[], directed_edges=[], nodes=[])
        blobs_list = list(blobs(net))
        # Empty network should have no blobs
        assert len(blobs_list) == 0
        assert net.number_of_nodes() == 0

    def test_single_node_network(self) -> None:
        """Test that single-node network has one trivial blob."""
        with pytest.warns(UserWarning, match="Single-node network detected"):
            net = SemiDirectedPhyNetwork(nodes=[(1, {'label': 'A'})])
        blobs_list = list(blobs(net))
        # Single node should be one trivial blob
        assert len(blobs_list) == 1
        assert {1} in blobs_list
        assert all(len(b) == 1 for b in blobs_list)
        # With filtering
        assert len(list(blobs(net, trivial=False, leaves=False))) == 0
        assert len(list(blobs(net, leaves=False))) == 0  # Single node is a leaf

    def test_single_edge_network(self) -> None:
        """Test that single-edge network has two trivial blobs."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(1, 2)],
            nodes=[(2, {'label': 'A'})]
        )
        blobs_list = list(blobs(net))
        # Single edge should create two trivial blobs (one per node)
        assert len(blobs_list) == 2
        assert {1} in blobs_list
        assert {2} in blobs_list
        assert all(len(b) == 1 for b in blobs_list)
        # With filtering
        assert len(list(blobs(net, trivial=False, leaves=False))) == 0
        # In SemiDirectedPhyNetwork, both nodes are leaves (no directed edges)
        assert len(list(blobs(net, leaves=False))) == 0  # Both nodes are leaves

    def test_tree_all_trivial_blobs(self) -> None:
        """Test that trees have all trivial blobs."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (4, {'label': 'C'}),
            ],
        )
        blobs_list = list(blobs(net))
        # All nodes should be trivial blobs
        assert len(blobs_list) == net.number_of_nodes()
        assert all(len(b) == 1 for b in blobs_list)

    def test_hybrid_network_non_trivial_blob(self) -> None:
        """Test network with hybrid node creating non-trivial blob."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {'u': 6, 'v': 5, 'gamma': 0.6},
                {'u': 7, 'v': 5, 'gamma': 0.4},
            ],
            undirected_edges=[
                (5, 1),
                (6, 2), (6, 3), (6, 7),
                (7, 8), (7, 9),
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (8, {'label': 'D'}),
                (9, {'label': 'E'}),
            ],
        )
        blobs_list = list(blobs(net))
        # Should have one non-trivial blob and several trivial blobs
        non_trivial = [b for b in blobs_list if len(b) > 1]
        trivial = [b for b in blobs_list if len(b) == 1]
        assert len(non_trivial) == 1
        assert len(trivial) == 5
        # Non-trivial blob should contain the hybrid and its parents
        assert {5, 6, 7} in blobs_list

    def test_parallel_edges_two_node_blob(self) -> None:
        """Test network with parallel edges creating 2-node blob."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {'u': 6, 'v': 5, 'gamma': 0.6},
                {'u': 7, 'v': 5, 'gamma': 0.4},
            ],
            undirected_edges=[
                (5, 6), (5, 6), (5, 6),  # Parallel edges
                (5, 1),
                (6, 2), (6, 3), (6, 7),
                (7, 8), (7, 9),
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (8, {'label': 'D'}),
                (9, {'label': 'E'}),
            ],
        )
        blobs_list = list(blobs(net))
        # Should have non-trivial blob containing nodes with parallel edges
        non_trivial = [b for b in blobs_list if len(b) > 1]
        assert len(non_trivial) == 1
        # The blob should include nodes 5, 6, 7 (biconnected component)
        assert {5, 6, 7} in blobs_list

    def test_filtering_trivial_false(self) -> None:
        """Test filtering with trivial=False."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {'u': 6, 'v': 5, 'gamma': 0.6},
                {'u': 7, 'v': 5, 'gamma': 0.4},
            ],
            undirected_edges=[
                (5, 1),
                (6, 2), (6, 3), (6, 7),
                (7, 8), (7, 9),
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (8, {'label': 'D'}),
                (9, {'label': 'E'}),
            ],
        )
        # trivial=False requires leaves=False
        non_trivial = list(blobs(net, trivial=False, leaves=False))
        assert len(non_trivial) == 1
        assert all(len(b) > 1 for b in non_trivial)

    def test_filtering_leaves_false(self) -> None:
        """Test filtering with leaves=False."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {'u': 6, 'v': 5, 'gamma': 0.6},
                {'u': 7, 'v': 5, 'gamma': 0.4},
            ],
            undirected_edges=[
                (5, 1),
                (6, 2), (6, 3), (6, 7),
                (7, 8), (7, 9),
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (8, {'label': 'D'}),
                (9, {'label': 'E'}),
            ],
        )
        without_leaves = list(blobs(net, leaves=False))
        # Should exclude blobs containing only leaves
        assert len(without_leaves) == 1
        assert all(not b.issubset(net.leaves) for b in without_leaves)

    def test_error_trivial_false_leaves_true(self) -> None:
        """Test that trivial=False and leaves=True raises error."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (4, {'label': 'C'}),
            ],
        )
        with pytest.raises(ValueError, match="Cannot have trivial=False and leaves=True"):
            list(blobs(net, trivial=False, leaves=True))

    def test_all_nodes_covered(self) -> None:
        """Test that all nodes are covered by blobs."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {'u': 6, 'v': 5, 'gamma': 0.6},
                {'u': 7, 'v': 5, 'gamma': 0.4},
            ],
            undirected_edges=[
                (5, 1),
                (6, 2), (6, 3), (6, 7),
                (7, 8), (7, 9),
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (8, {'label': 'D'}),
                (9, {'label': 'E'}),
            ],
        )
        blobs_list = list(blobs(net))
        all_nodes = set()
        for blob in blobs_list:
            all_nodes.update(blob)
        assert all_nodes == set(net._graph.nodes())

