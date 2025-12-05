"""
Tests for the network module.
"""

import warnings

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
        # Empty network is invalid (no root), so test with a minimal valid network
        network = DirectedPhyNetwork(edges=[(1, 2)], taxa={2: "A"})
        assert network.number_of_nodes() == 2
        assert network.number_of_edges() == 1

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


class TestIncidentEdges:
    """Test cases for incident edge methods in DirectedPhyNetwork."""

    def test_incident_parent_edges_basic(self) -> None:
        """Test incident_parent_edges with basic edges."""
        # Create valid structure: root -> tree nodes (out-degree 2) -> hybrid -> leaf
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            taxa={2: "A", 8: "B", 9: "C"}
        )
        
        parent_edges = list(net.incident_parent_edges(4))  # Test on hybrid node
        assert len(parent_edges) == 2
        assert (5, 4) in parent_edges
        assert (6, 4) in parent_edges

    def test_incident_parent_edges_with_attributes(self) -> None:
        """Test incident_parent_edges with edge attributes."""
        # Valid network: root -> tree nodes (out-degree 2) -> hybrid -> leaf
        net = DirectedPhyNetwork(
            edges=[
                {'u': 7, 'v': 5},
                {'u': 7, 'v': 6},
                {'u': 5, 'v': 4, 'branch_length': 0.5, 'bootstrap': 0.95},
                {'u': 5, 'v': 8},
                {'u': 6, 'v': 4, 'branch_length': 0.3, 'bootstrap': 0.87},
                {'u': 6, 'v': 9},
                {'u': 4, 'v': 2}
            ],
            taxa={2: "A", 8: "B", 9: "C"}
        )
        
        parent_edges = list(net.incident_parent_edges(4, data=True))  # Test on hybrid node
        assert len(parent_edges) == 2
        # Check that data is included
        edge_dict = {edge[0]: edge[2] for edge in parent_edges if len(edge) == 3}
        assert 5 in edge_dict
        assert edge_dict[5]['branch_length'] == 0.5
        assert edge_dict[5]['bootstrap'] == 0.95

    def test_incident_parent_edges_with_keys_and_data(self) -> None:
        """Test incident_parent_edges with keys and data."""
        # Valid network with parallel edges to hybrid
        net = DirectedPhyNetwork(
            edges=[
                {'u': 7, 'v': 5},
                {'u': 7, 'v': 6},
                {'u': 5, 'v': 4, 'key': 0, 'branch_length': 0.5},
                {'u': 5, 'v': 4, 'key': 1, 'branch_length': 0.7},  # Parallel edge
                {'u': 5, 'v': 8},
                {'u': 6, 'v': 4},
                {'u': 6, 'v': 9},
                {'u': 4, 'v': 2}
            ],
            taxa={2: "A", 8: "B", 9: "C"}
        )
        
        parent_edges = list(net.incident_parent_edges(4, keys=True, data=True))  # Test on hybrid
        assert len(parent_edges) >= 2  # At least 2 edges from 5, plus 1 from 6
        # Check structure
        for edge in parent_edges:
            assert len(edge) == 4
            u, v, key, data = edge
            assert v == 4
            assert 'branch_length' in data or u == 6  # Edge from 6 may not have branch_length

    def test_incident_child_edges_basic(self) -> None:
        """Test incident_child_edges with basic edges."""
        # Valid tree: root -> internal -> leaves
        net = DirectedPhyNetwork(
            edges=[(1, 2), (1, 3), (1, 4)],
            taxa={2: "A", 3: "B", 4: "C"}
        )
        
        child_edges = list(net.incident_child_edges(1))
        assert len(child_edges) == 3
        assert (1, 2) in child_edges
        assert (1, 3) in child_edges
        assert (1, 4) in child_edges

    def test_incident_child_edges_with_attributes(self) -> None:
        """Test incident_child_edges with edge attributes."""
        # Valid tree: root -> leaves
        net = DirectedPhyNetwork(
            edges=[
                {'u': 1, 'v': 2, 'branch_length': 0.5, 'bootstrap': 0.95},
                {'u': 1, 'v': 3, 'branch_length': 0.3, 'bootstrap': 0.87}
            ],
            taxa={2: "A", 3: "B"}
        )
        
        child_edges = list(net.incident_child_edges(1, data=True))
        assert len(child_edges) == 2
        # Check that data is included
        edge_dict = {edge[1]: edge[2] for edge in child_edges if len(edge) == 3}
        assert 2 in edge_dict
        assert edge_dict[2]['branch_length'] == 0.5
        assert edge_dict[2]['bootstrap'] == 0.95

    def test_incident_child_edges_with_keys_and_data(self) -> None:
        """Test incident_child_edges with keys and data."""
        # Valid network: root -> tree nodes -> hybrid (with parallel edges from one parent)
        net = DirectedPhyNetwork(
            edges=[
                {'u': 7, 'v': 5},
                {'u': 7, 'v': 6},
                {'u': 5, 'v': 4, 'key': 0, 'branch_length': 0.5},
                {'u': 5, 'v': 4, 'key': 1, 'branch_length': 0.7},  # Parallel edge to hybrid
                {'u': 5, 'v': 8},
                {'u': 6, 'v': 4},
                {'u': 6, 'v': 9},
                {'u': 4, 'v': 2}
            ],
            taxa={2: "A", 8: "B", 9: "C"}
        )
        
        child_edges = list(net.incident_child_edges(5, keys=True, data=True))  # Test on tree node 5
        # Should have 2 parallel edges to hybrid node 4, plus 1 to 8
        parallel_to_4 = [e for e in child_edges if e[1] == 4]
        assert len(parallel_to_4) == 2  # Exactly 2 parallel edges to 4
        for edge in parallel_to_4:
            assert len(edge) == 4
            u, v, key, data = edge
            assert u == 5
            assert v == 4
            assert key in [0, 1]
            assert 'branch_length' in data

    def test_incident_edges_for_hybrid_node(self) -> None:
        """Test incident edges for a hybrid node with gamma values."""
        # Valid network: root -> tree nodes -> hybrid -> leaf
        net = DirectedPhyNetwork(
            edges=[
                {'u': 7, 'v': 5},
                {'u': 7, 'v': 6},
                {'u': 5, 'v': 4, 'gamma': 0.6},
                {'u': 5, 'v': 8},
                {'u': 6, 'v': 4, 'gamma': 0.4},
                {'u': 6, 'v': 9},
                {'u': 4, 'v': 2}
            ],
            taxa={2: "A", 8: "B", 9: "C"}
        )
        
        # Check parent edges (incoming to hybrid node)
        parent_edges = list(net.incident_parent_edges(4, data=True))
        assert len(parent_edges) == 2
        gamma_values = []
        for edge in parent_edges:
            if len(edge) == 3:
                _, _, data = edge
                if 'gamma' in data:
                    gamma_values.append(data['gamma'])
        assert len(gamma_values) == 2
        assert 0.6 in gamma_values
        assert 0.4 in gamma_values
        
        # Check child edges (outgoing from hybrid node)
        child_edges = list(net.incident_child_edges(4))
        assert len(child_edges) == 1
        assert (4, 2) in child_edges

    def test_incident_edges_empty(self) -> None:
        """Test incident edges for nodes with no edges."""
        net = DirectedPhyNetwork(
            edges=[(1, 2)],
            taxa={2: "A"}
        )
        
        # Node 1 has no incoming edges
        parent_edges_1 = list(net.incident_parent_edges(1))
        assert len(parent_edges_1) == 0
        
        # Node 2 has no outgoing edges
        child_edges_2 = list(net.incident_child_edges(2))
        assert len(child_edges_2) == 0

    def test_incident_edges_consistency_with_parents_children(self) -> None:
        """Test that incident edges are consistent with parents/children methods."""
        # Valid network: root -> tree nodes (out-degree 2) -> hybrid -> leaf
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 2), (5, 8), (6, 2), (6, 9), (2, 4)],
            taxa={4: "A", 8: "B", 9: "C"}
        )
        
        # Check parent edges match parents (node 2 is hybrid)
        parents = set(net.parents(2))
        parent_edges = list(net.incident_parent_edges(2))
        parent_nodes_from_edges = {edge[0] for edge in parent_edges}
        assert parents == parent_nodes_from_edges
        
        # Check child edges match children
        children = set(net.children(2))
        child_edges = list(net.incident_child_edges(2))
        child_nodes_from_edges = {edge[1] for edge in child_edges}
        assert children == child_nodes_from_edges
