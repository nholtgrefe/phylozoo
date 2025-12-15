"""
Tests for blob computation in DirectedPhyNetwork.
"""

import warnings

import pytest

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.features import blobs


class TestBlobs:
    """Test blob computation for DirectedPhyNetwork."""

    def test_empty_network(self) -> None:
        """Test that empty network has no blobs."""
        with pytest.warns(UserWarning, match="Empty network.*no nodes"):
            net = DirectedPhyNetwork(edges=[], nodes=[])
        blobs_list = list(blobs(net))
        # Empty network should have no blobs
        assert len(blobs_list) == 0
        assert net.number_of_nodes() == 0

    def test_single_node_network(self) -> None:
        """Test that single-node network has one trivial blob."""
        with pytest.warns(UserWarning, match="Single-node network detected"):
            net = DirectedPhyNetwork(nodes=[(1, {'label': 'A'})])
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
        net = DirectedPhyNetwork(
            edges=[(1, 2)],
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
        assert len(list(blobs(net, leaves=False))) == 1  # Only root node remains

    def test_tree_all_trivial_blobs(self) -> None:
        """Test that trees have all trivial blobs."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 1), (5, 2), (6, 3), (6, 4)],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (4, {'label': 'D'}),
            ],
        )
        blobs_list = list(blobs(net))
        # All nodes should be trivial blobs
        assert len(blobs_list) == net.number_of_nodes()
        assert all(len(b) == 1 for b in blobs_list)

    def test_hybrid_network_non_trivial_blob(self) -> None:
        """Test network with hybrid node creating non-trivial blob."""
        net = DirectedPhyNetwork(
            edges=[
                (8, 5), (8, 6),
                (5, 1), (5, 2),
                (6, 3), (6, 9),
                (5, 4), (6, 4),
                (4, 7),
                (7, 10), (7, 11),
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (9, {'label': 'D'}),
                (10, {'label': 'E'}),
                (11, {'label': 'F'}),
            ],
        )
        blobs_list = list(blobs(net))
        # Should have one non-trivial blob and several trivial blobs
        non_trivial = [b for b in blobs_list if len(b) > 1]
        trivial = [b for b in blobs_list if len(b) == 1]
        assert len(non_trivial) == 1
        assert len(trivial) == 7
        # Non-trivial blob should contain the hybrid and its parents
        assert {4, 5, 6, 8} in blobs_list

    def test_parallel_edges_two_node_blob(self) -> None:
        """Test network with parallel edges creating 2-node blob."""
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),
                (5, 4, 0), (5, 4, 1), (5, 4, 2),  # Parallel edges
                (6, 4),
                (4, 8),
                (5, 1), (5, 2),
                (6, 3),
                (8, 9), (8, 10),
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (9, {'label': 'D'}),
                (10, {'label': 'E'}),
            ],
        )
        blobs_list = list(blobs(net))
        # Should have non-trivial blob containing nodes with parallel edges
        non_trivial = [b for b in blobs_list if len(b) > 1]
        assert len(non_trivial) == 1
        # The blob should include nodes 4, 5, 6, 7 (biconnected component)
        assert {4, 5, 6, 7} in blobs_list

    def test_filtering_trivial_false(self) -> None:
        """Test filtering with trivial=False."""
        net = DirectedPhyNetwork(
            edges=[
                (8, 5), (8, 6),
                (5, 1), (5, 2),
                (6, 3), (6, 9),
                (5, 4), (6, 4),
                (4, 7),
                (7, 10), (7, 11),
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (9, {'label': 'D'}),
                (10, {'label': 'E'}),
                (11, {'label': 'F'}),
            ],
        )
        # trivial=False requires leaves=False
        non_trivial = list(blobs(net, trivial=False, leaves=False))
        assert len(non_trivial) == 1
        assert all(len(b) > 1 for b in non_trivial)

    def test_filtering_leaves_false(self) -> None:
        """Test filtering with leaves=False."""
        net = DirectedPhyNetwork(
            edges=[
                (8, 5), (8, 6),
                (5, 1), (5, 2),
                (6, 3), (6, 9),
                (5, 4), (6, 4),
                (4, 7),
                (7, 10), (7, 11),
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (9, {'label': 'D'}),
                (10, {'label': 'E'}),
                (11, {'label': 'F'}),
            ],
        )
        without_leaves = list(blobs(net, leaves=False))
        # Should exclude blobs containing only leaves
        assert len(without_leaves) == 2
        assert all(not b.issubset(net.leaves) for b in without_leaves)

    def test_error_trivial_false_leaves_true(self) -> None:
        """Test that trivial=False and leaves=True raises error."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 1), (5, 2), (6, 3), (6, 4)],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (4, {'label': 'D'}),
            ],
        )
        with pytest.raises(ValueError, match="Cannot have trivial=False and leaves=True"):
            list(blobs(net, trivial=False, leaves=True))

    def test_all_nodes_covered(self) -> None:
        """Test that all nodes are covered by blobs."""
        net = DirectedPhyNetwork(
            edges=[
                (8, 5), (8, 6),
                (5, 1), (5, 2),
                (6, 3), (6, 9),
                (5, 4), (6, 4),
                (4, 7),
                (7, 10), (7, 11),
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (9, {'label': 'D'}),
                (10, {'label': 'E'}),
                (11, {'label': 'F'}),
            ],
        )
        blobs_list = list(blobs(net))
        all_nodes = set()
        for blob in blobs_list:
            all_nodes.update(blob)
        assert all_nodes == set(net._graph.nodes())

