"""
Tests for the network module.

This module contains tests for SemiDirectedNetwork.
Tests for DirectedPhyNetwork are in tests/core/network/dnetwork/.
"""

import warnings

import pytest
from phylozoo.core.network import (
    SemiDirectedNetwork,
    random_semi_directed_network,
)


class TestSemiDirectedNetwork:
    """Test cases for the SemiDirectedNetwork class."""

    def test_semi_directed_network_creation(self) -> None:
        """Test creating a semi-directed network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            network = SemiDirectedNetwork()
        assert network.number_of_nodes() == 0
        assert network.number_of_edges() == 0
        assert len(network.undirected_edges) == 0

    def test_semi_directed_network_is_immutable(self) -> None:
        """Test that semi-directed network is immutable."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            network = SemiDirectedNetwork()
        # Verify mutation methods don't exist
        assert not hasattr(network, "add_undirected_edge")
        assert not hasattr(network, "add_node")
        assert not hasattr(network, "add_edge")

    def test_random_semi_directed_network(self) -> None:
        """Test random semi-directed network generation."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            network = random_semi_directed_network(5)
        assert isinstance(network, SemiDirectedNetwork)


