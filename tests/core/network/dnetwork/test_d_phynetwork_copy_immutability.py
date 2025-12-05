"""
Comprehensive tests for DirectedPhyNetwork copy and immutability.

This module tests:
- copy() method
- Immutability verification
- Cached property behavior on copy
"""

import warnings

import pytest

from phylozoo.core.network import DirectedPhyNetwork


class TestCopy:
    """Test cases for copy() method."""

    def test_copy_independent_objects(self) -> None:
        """Test that copy creates independent objects."""
        net1 = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"})
        net2 = net1.copy()
        
        # Should be different objects
        assert net1 is not net2
        assert net1._graph is not net2._graph
        assert net1._node_to_label is not net2._node_to_label
        assert net1._label_to_node is not net2._label_to_node

    def test_copy_labels_copied(self) -> None:
        """Test that labels are copied."""
        net1 = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            taxa={1: "A", 2: "B"},
            internal_node_labels={3: "root"}
        )
        net2 = net1.copy()
        
        # Labels should be copied
        assert net2.get_label(1) == "A"
        assert net2.get_label(2) == "B"
        assert net2.get_label(3) == "root"

    def test_copy_graph_structure_copied(self) -> None:
        """Test that graph structure is copied."""
        net1 = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"})
        net2 = net1.copy()
        
        # Graph structure should be copied
        assert net2.number_of_nodes() == net1.number_of_nodes()
        assert net2.number_of_edges() == net1.number_of_edges()
        assert net2.has_edge(3, 1)
        assert net2.has_edge(3, 2)

    def test_copy_properties_equal(self) -> None:
        """Test that copied network has equal properties."""
        net1 = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            taxa={2: "A", 8: "B", 9: "C"}
        )
        net2 = net1.copy()
        
        # All properties should be equal
        assert net1.leaves == net2.leaves
        assert net1.taxa == net2.taxa
        assert net1.root_node == net2.root_node
        assert net1.internal_nodes == net2.internal_nodes
        assert net1.hybrid_nodes == net2.hybrid_nodes
        assert net1.tree_nodes == net2.tree_nodes
        assert net1.hybrid_edges == net2.hybrid_edges
        assert net1.tree_edges == net2.tree_edges

    def test_copy_cached_properties_recomputed(self) -> None:
        """Test that cached properties are recomputed on copy."""
        net1 = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"})
        # Access to cache
        _ = net1.leaves
        _ = net1.taxa
        
        net2 = net1.copy()
        
        # Properties should be equal but may be different objects (recomputed)
        assert net1.leaves == net2.leaves
        assert net1.taxa == net2.taxa

    def test_copy_with_edge_attributes(self) -> None:
        """Test copy with edge attributes."""
        net1 = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'branch_length': 0.5, 'bootstrap': 0.95},
                {'u': 3, 'v': 2, 'branch_length': 0.3}
            ],
            taxa={1: "A", 2: "B"}
        )
        net2 = net1.copy()
        
        # Edge attributes should be copied
        assert net2.get_branch_length(3, 1) == 0.5
        assert net2.get_bootstrap(3, 1) == 0.95
        assert net2.get_branch_length(3, 2) == 0.3

    def test_copy_with_hybrid_gamma(self) -> None:
        """Test copy with hybrid nodes and gamma values."""
        net1 = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),  # Root to tree nodes
                {'u': 5, 'v': 4, 'gamma': 0.6},
                {'u': 5, 'v': 8},  # Tree node 5 also has another child
                {'u': 6, 'v': 4, 'gamma': 0.4},
                {'u': 6, 'v': 9},  # Tree node 6 also has another child
                {'u': 4, 'v': 1}
            ],
            taxa={1: "A", 8: "B", 9: "C"}
        )
        net2 = net1.copy()
        
        # Gamma values should be copied
        assert net2.get_gamma(5, 4) == 0.6
        assert net2.get_gamma(6, 4) == 0.4
        assert net2.validate() is True


class TestImmutability:
    """Test cases for network immutability."""

    def test_no_mutation_methods(self) -> None:
        """Test that mutation methods don't exist."""
        net = DirectedPhyNetwork(edges=[(3, 1)], taxa={1: "A"})
        
        # Verify mutation methods don't exist
        assert not hasattr(net, "add_node")
        assert not hasattr(net, "add_edge")
        assert not hasattr(net, "remove_node")
        assert not hasattr(net, "remove_edge")
        assert not hasattr(net, "clear")
        assert not hasattr(net, "set_label")
        assert not hasattr(net, "update")
        assert not hasattr(net, "add_nodes_from")
        assert not hasattr(net, "add_edges_from")

    def test_graph_not_directly_modifiable(self) -> None:
        """Test that _graph is not directly modifiable (should raise error if tried)."""
        net = DirectedPhyNetwork(edges=[(3, 1)], taxa={1: "A"})
        
        # _graph exists but should not be modified directly
        # Attempting to modify would require accessing private attributes
        # which is discouraged, but we can verify the structure is immutable
        original_nodes = set(net._graph.nodes)
        original_edges = list(net._graph.edges())
        
        # Create copy to verify original unchanged
        net2 = net.copy()
        assert set(net._graph.nodes) == original_nodes
        assert list(net._graph.edges()) == original_edges

    def test_labels_not_directly_modifiable(self) -> None:
        """Test that label mappings are not directly modifiable."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            taxa={1: "A", 2: "B"},
            internal_node_labels={3: "root"}
        )
        
        original_labels = net._node_to_label.copy()
        original_label_to_node = net._label_to_node.copy()
        
        # Create copy
        net2 = net.copy()
        
        # Original should be unchanged
        assert net._node_to_label == original_labels
        assert net._label_to_node == original_label_to_node

    def test_immutability_after_operations(self) -> None:
        """Test that network remains immutable after operations."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"})
        
        # Perform various read operations
        _ = net.number_of_nodes()
        _ = net.number_of_edges()
        _ = net.leaves
        _ = net.taxa
        _ = net.root_node
        _ = net.has_edge(3, 1)
        _ = net.degree(3)
        _ = list(net.parents(1))
        _ = list(net.children(3))
        
        # Network should still be immutable
        assert not hasattr(net, "add_node")
        assert not hasattr(net, "add_edge")

    def test_copy_immutability_preserved(self) -> None:
        """Test that copy preserves immutability."""
        net1 = DirectedPhyNetwork(edges=[(3, 1)], taxa={1: "A"})
        net2 = net1.copy()
        
        # Copy should also be immutable
        assert not hasattr(net2, "add_node")
        assert not hasattr(net2, "add_edge")

