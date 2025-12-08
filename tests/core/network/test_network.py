"""
Tests for the network module.

This module contains tests for SemiDirectedPhyNetwork.
Tests for DirectedPhyNetwork are in tests/core/network/dnetwork/.
"""

import warnings

import pytest
from phylozoo.core.network import (
    SemiDirectedPhyNetwork,
)


class TestSemiDirectedPhyNetwork:
    """Test cases for the SemiDirectedPhyNetwork class."""

    def test_semi_directed_network_creation(self) -> None:
        """Test creating a semi-directed network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            network = SemiDirectedPhyNetwork(
                directed_edges=[],
                undirected_edges=[]
            )
        assert network.number_of_nodes() == 0
        assert network.number_of_edges() == 0
        assert len(network.tree_edges) == 0

    def test_semi_directed_network_is_immutable(self) -> None:
        """Test that semi-directed network is immutable."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            network = SemiDirectedPhyNetwork(
                directed_edges=[],
                undirected_edges=[]
            )
        # Verify mutation methods don't exist
        assert not hasattr(network, "add_undirected_edge")
        assert not hasattr(network, "add_node")
        assert not hasattr(network, "add_edge")


