"""
Tests for LSA-related utilities on DirectedPhyNetwork.
"""

import warnings

import pytest

from phylozoo.core.network import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.operations import (
    find_lsa_node,
    to_LSA_network,
)
from phylozoo.core.network.dnetwork.classifications import is_LSA_network


class TestFindLSANode:
    """Tests for find_lsa_node."""

    def test_simple_tree_root_is_lsa(self) -> None:
        """
        A simple binary tree should have its root as the LSA.
        """
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})])
        assert find_lsa_node(net) == 3

    def test_hybrid_lsa_below_root(self) -> None:
        """
        Network where all leaves are below a hybrid node; the deepest common ancestor is the LSA.
        """
        edges = [
            (7, 5),
            (7, 6),
            (5, 4, 0),
            (5, 4, 1),
            (6, 4, 0),
            (6, 4, 1),
            (4, 10),
            (10, 8),
            (10, 9),
        ]
        net = DirectedPhyNetwork(edges=edges, nodes=[(8, {"label": "A"}), (9, {"label": "B"})])
        # Common ancestors of all leaves are {7, 4, 10}; depth(10) > depth(4) > depth(7)
        assert find_lsa_node(net) == 10

    def test_empty_network_raises(self) -> None:
        """
        Empty network has no LSA node.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)  # Ignore validation warning
            with pytest.warns(UserWarning, match="empty edges list"):
                with pytest.raises(ValueError, match="empty network"):
                    find_lsa_node(DirectedPhyNetwork(edges=[]))


class TestToLSANetwork:
    """Tests for to_LSA_network."""

    def test_returns_copy_if_root_is_lsa(self) -> None:
        """
        When root is already LSA, to_LSA_network returns a copy rooted at the same node.
        """
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})])
        lsa_net = to_LSA_network(net)
        assert lsa_net.root_node == 3
        assert set(lsa_net.leaves) == {1, 2}

    def test_trims_above_lsa(self) -> None:
        """
        to_LSA_network removes nodes/edges above LSA and preserves leaves and attributes below.
        """
        edges = [
            (7, 5),
            (7, 6),
            (5, 4, 0),
            (5, 4, 1),
            (6, 4, 0),
            (6, 4, 1),
            (4, 10),
            (10, 8),
            (10, 9),
        ]
        net = DirectedPhyNetwork(edges=edges, nodes=[(8, {"label": "A"}), (9, {"label": "B"})])

        lsa_net = to_LSA_network(net)

        # LSA becomes the root (deepest common ancestor is node 10)
        assert lsa_net.root_node == 10
        # Leaves preserved
        assert set(lsa_net.leaves) == {8, 9}
        # All nodes in the trimmed network are descendants of LSA
        assert 7 not in lsa_net._graph.nodes
        assert 5 not in lsa_net._graph.nodes
        assert 6 not in lsa_net._graph.nodes
        assert 4 not in lsa_net._graph.nodes


class TestIsLSANetwork:
    """Tests for is_LSA_network classification."""

    def test_empty_network_true(self) -> None:
        """
        Empty network counts as LSA network.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)  # Ignore validation warning
            with pytest.warns(UserWarning, match="empty edges list"):
                net = DirectedPhyNetwork(edges=[])
        assert is_LSA_network(net) is True

    def test_root_is_lsa_true(self) -> None:
        """
        Tree with root as LSA returns True.
        """
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})])
        assert is_LSA_network(net) is True

    def test_lsa_below_root_false(self) -> None:
        """
        Network with LSA below root returns False.
        """
        edges = [
            (7, 5),
            (7, 6),
            (5, 4, 0),
            (5, 4, 1),
            (6, 4, 0),
            (6, 4, 1),
            (4, 10),
            (10, 8),
            (10, 9),
        ]
        net = DirectedPhyNetwork(edges=edges, nodes=[(8, {"label": "A"}), (9, {"label": "B"})])
        assert is_LSA_network(net) is False

