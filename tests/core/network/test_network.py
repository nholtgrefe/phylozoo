"""
Tests for the network module.
"""

import pytest
from phylozoo.core.network import (
    DirectedNetwork,
    SemiDirectedNetwork,
    random_semi_directed_network,
)


class TestDirectedNetwork:
    """Test cases for the DirectedNetwork class."""

    def test_directed_network_creation(self) -> None:
        """Test creating a directed network."""
        network = DirectedNetwork()
        assert len(network.nodes) == 0
        assert len(network.edges) == 0

    def test_add_node(self) -> None:
        """Test adding a node to the network."""
        network = DirectedNetwork()
        network.add_node("A")
        assert "A" in network.nodes

    def test_add_edge(self) -> None:
        """Test adding an edge to the network."""
        network = DirectedNetwork()
        network.add_edge("A", "B")
        assert "A" in network.nodes
        assert "B" in network.nodes
        assert ("A", "B") in network.edges

    def test_directed_network_repr(self) -> None:
        """Test string representation of a directed network."""
        network = DirectedNetwork()
        repr_str = repr(network)
        assert "DirectedNetwork" in repr_str


class TestSemiDirectedNetwork:
    """Test cases for the SemiDirectedNetwork class."""

    def test_semi_directed_network_creation(self) -> None:
        """Test creating a semi-directed network."""
        network = SemiDirectedNetwork()
        assert len(network.nodes) == 0
        assert len(network.edges) == 0
        assert len(network.undirected_edges) == 0

    def test_add_undirected_edge(self) -> None:
        """Test adding an undirected edge to the network."""
        network = SemiDirectedNetwork()
        network.add_undirected_edge("A", "B")
        assert "A" in network.nodes
        assert "B" in network.nodes
        assert ("A", "B") in network.undirected_edges

    def test_random_semi_directed_network(self) -> None:
        """Test generating a random semi-directed network."""
        network = random_semi_directed_network(n_leaves=4, seed=42)
        assert isinstance(network, SemiDirectedNetwork)
