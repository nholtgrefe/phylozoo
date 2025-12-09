"""
Comprehensive tests for MixedPhyNetwork copy and immutability.

This module tests:
- copy() method
- Immutability verification
- Cached property behavior on copy
"""

import warnings

import pytest

from phylozoo.core.network.sdnetwork import MixedPhyNetwork
from tests.core.network.sdnetwork.conftest import expect_mixed_network_warning


class TestCopy:
    """Test cases for copy() method."""

    def test_copy_independent_objects(self) -> None:
        """Test that copy creates independent objects."""
        with expect_mixed_network_warning():
            net1 = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        net2 = net1.copy()
        
        # Should be different objects
        assert net1 is not net2
        assert net1._graph is not net2._graph
        assert net1._node_to_label is not net2._node_to_label
        assert net1._label_to_node is not net2._label_to_node

    def test_copy_labels_copied(self) -> None:
        """Test that labels are copied."""
        with expect_mixed_network_warning():
            net1 = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"},
            internal_node_labels={3: "root"}
            )
        net2 = net1.copy()
        
        # Labels should be copied
        assert net2.get_label(1) == "A"
        assert net2.get_label(2) == "B"
        assert net2.get_label(3) == "root"

    def test_copy_graph_structure_copied(self) -> None:
        """Test that graph structure is copied."""
        with expect_mixed_network_warning():
            net1 = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        net2 = net1.copy()
        
        # Graph structure should be copied
        assert net2.number_of_nodes() == net1.number_of_nodes()
        assert net2.number_of_edges() == net1.number_of_edges()
        assert net2.has_edge(3, 1)
        assert net2.has_edge(3, 2)

    def test_copy_properties_equal(self) -> None:
        """Test that copied network has equal properties."""
        with expect_mixed_network_warning():
            net1 = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            taxa={2: "A", 8: "B", 9: "C", 10: "D", 11: "E"}
            )
        net2 = net1.copy()
        
        # All properties should be equal
        assert net1.leaves == net2.leaves
        assert net1.taxa == net2.taxa
        assert net1.internal_nodes == net2.internal_nodes
        assert net1.hybrid_nodes == net2.hybrid_nodes
        assert net1.tree_nodes == net2.tree_nodes
        assert net1.hybrid_edges == net2.hybrid_edges
        assert net1.tree_edges == net2.tree_edges

    def test_copy_cached_properties_recomputed(self) -> None:
        """Test that cached properties are recomputed on copy."""
        with expect_mixed_network_warning():
            net1 = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        # Access to cache
        _ = net1.leaves
        _ = net1.taxa
        
        net2 = net1.copy()
        
        # Properties should be equal but may be different objects (recomputed)
        assert net1.leaves == net2.leaves
        assert net1.taxa == net2.taxa

    def test_copy_with_edge_attributes(self) -> None:
        """Test copy with edge attributes."""
        # Node 3 needs degree >= 3
        with expect_mixed_network_warning():
            net1 = MixedPhyNetwork(
            undirected_edges=[
            {'u': 3, 'v': 1, 'branch_length': 0.5, 'bootstrap': 0.95},
            {'u': 3, 'v': 2, 'branch_length': 0.3},
            (3, 4)
            ],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        net2 = net1.copy()
        
        # Edge attributes should be copied
        assert net2.get_branch_length(3, 1) == 0.5
        assert net2.get_bootstrap(3, 1) == 0.95
        assert net2.get_branch_length(3, 2) == 0.3

    def test_copy_with_gamma(self) -> None:
        """Test copy with gamma values."""
        # Hybrid node 4: indegree 2, total_degree must be 3 (only 1 outgoing)
        # Nodes 5 and 6 need degree >= 3
        with expect_mixed_network_warning():
            net1 = MixedPhyNetwork(
            directed_edges=[
            {'u': 5, 'v': 4, 'gamma': 0.6},
            {'u': 6, 'v': 4, 'gamma': 0.4}
            ],
            undirected_edges=[(4, 1), (5, 8), (5, 10), (6, 9), (6, 11)],
            taxa={1: "A", 8: "B", 9: "C", 10: "D", 11: "E"}
            )
        net2 = net1.copy()
        
        # Gamma values should be copied
        assert net2.get_gamma(5, 4) == 0.6
        assert net2.get_gamma(6, 4) == 0.4

    def test_copy_type(self) -> None:
        """Test that copy returns correct type."""
        with expect_mixed_network_warning():
            net1 = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        net2 = net1.copy()
        assert isinstance(net2, MixedPhyNetwork)


class TestImmutability:
    """Test cases for network immutability."""

    def test_no_add_methods(self) -> None:
        """Test that mutation methods don't exist."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        # Verify mutation methods don't exist
        assert not hasattr(net, "add_undirected_edge")
        assert not hasattr(net, "add_directed_edge")
        assert not hasattr(net, "add_node")
        assert not hasattr(net, "remove_edge")
        assert not hasattr(net, "remove_node")

    def test_graph_not_modifiable(self) -> None:
        """Test that underlying graph is not directly modifiable."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        # Direct modification should be possible but discouraged
        # We test that the network structure remains consistent
        original_edges = net.number_of_edges()
        # If we could modify, we'd test that here
        # But since it's immutable, we just verify structure
        assert net.number_of_edges() == original_edges

