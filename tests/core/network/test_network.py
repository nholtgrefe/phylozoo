"""
Tests for the network module.
"""

import pytest
from phylozoo.core.network import (
    DirectedPhyNetwork,
    SemiDirectedNetwork,
    random_semi_directed_network,
)


class TestDirectedPhyNetwork:
    """Test cases for the DirectedPhyNetwork class."""

    def test_directed_network_creation(self) -> None:
        """Test creating a directed network."""
        network = DirectedPhyNetwork(edges=[])
        assert network.number_of_nodes() == 0
        assert network.number_of_edges() == 0

    def test_directed_network_with_edges(self) -> None:
        """Test creating a directed network with edges and taxa."""
        network = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"})
        assert network.number_of_nodes() == 3
        assert network.number_of_edges() == 2
        assert 1 in network
        assert 2 in network
        assert 3 in network
        assert network.has_edge(3, 1)
        assert network.has_edge(3, 2)
        assert network.taxa == {"A", "B"}

    def test_directed_network_repr(self) -> None:
        """Test string representation of a directed network."""
        network = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"})
        repr_str = repr(network)
        assert "DirectedPhyNetwork" in repr_str
        assert "A" in repr_str or "B" in repr_str

    def test_network_is_immutable(self) -> None:
        """Test that network is immutable (no mutation methods)."""
        network = DirectedPhyNetwork(edges=[(3, 1)], taxa={1: "A"})
        # Verify mutation methods don't exist
        assert not hasattr(network, "add_node")
        assert not hasattr(network, "add_edge")
        assert not hasattr(network, "remove_node")
        assert not hasattr(network, "clear")
        assert not hasattr(network, "set_label")
        assert not hasattr(network, "no_validation")

    def test_cached_property_leaves(self) -> None:
        """Test that leaves property is cached."""
        network = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"})
        leaves1 = network.leaves
        leaves2 = network.leaves
        assert leaves1 == leaves2
        assert leaves1 == {1, 2}


class TestSemiDirectedNetwork:
    """Test cases for the SemiDirectedNetwork class."""

    def test_semi_directed_network_creation(self) -> None:
        """Test creating a semi-directed network."""
        network = SemiDirectedNetwork()
        assert network.number_of_nodes() == 0
        assert network.number_of_edges() == 0
        assert len(network.undirected_edges) == 0

    def test_semi_directed_network_is_immutable(self) -> None:
        """Test that semi-directed network is immutable."""
        network = SemiDirectedNetwork()
        # Verify mutation methods don't exist
        assert not hasattr(network, "add_undirected_edge")
        assert not hasattr(network, "add_node")
        assert not hasattr(network, "add_edge")

    def test_random_semi_directed_network(self) -> None:
        """Test random semi-directed network generation."""
        network = random_semi_directed_network(5)
        assert isinstance(network, SemiDirectedNetwork)
