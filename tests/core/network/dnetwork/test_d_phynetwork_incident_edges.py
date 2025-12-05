"""
Comprehensive tests for DirectedPhyNetwork incident edge methods.

This module tests all incident edge methods including:
- incident_parent_edges()
- incident_child_edges()
- Parallel edges handling
- Edge attributes in incident edges
- Consistency with parents()/children()
"""

import pytest

from phylozoo.core.network import DirectedPhyNetwork


class TestIncidentParentEdges:
    """Test cases for incident_parent_edges() method."""

    def test_incident_parent_edges_root(self) -> None:
        """Test incident_parent_edges for root (empty)."""
        net = DirectedPhyNetwork(edges=[(3, 1)], taxa={1: "A"})
        parent_edges = list(net.incident_parent_edges(3))
        assert len(parent_edges) == 0

    def test_incident_parent_edges_leaf(self) -> None:
        """Test incident_parent_edges for leaf."""
        net = DirectedPhyNetwork(edges=[(3, 1)], taxa={1: "A"})
        parent_edges = list(net.incident_parent_edges(1))
        assert len(parent_edges) == 1
        assert (3, 1) in parent_edges

    def test_incident_parent_edges_tree_node(self) -> None:
        """Test incident_parent_edges for tree node."""
        net = DirectedPhyNetwork(
            edges=[(4, 3), (3, 1), (3, 2)],
            taxa={1: "A", 2: "B"}
        )
        parent_edges = list(net.incident_parent_edges(3))
        assert len(parent_edges) == 1
        assert (4, 3) in parent_edges

    def test_incident_parent_edges_hybrid(self) -> None:
        """Test incident_parent_edges for hybrid node."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            taxa={2: "A", 8: "B", 9: "C"}
        )
        parent_edges = list(net.incident_parent_edges(4))
        assert len(parent_edges) == 2
        assert (5, 4) in parent_edges
        assert (6, 4) in parent_edges

    def test_incident_parent_edges_with_data(self) -> None:
        """Test incident_parent_edges with data=True."""
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
        parent_edges = list(net.incident_parent_edges(4, data=True))
        assert len(parent_edges) == 2
        # Check data is included
        for edge in parent_edges:
            assert len(edge) == 3  # (u, v, data)
            u, v, data = edge
            assert v == 4
            assert 'branch_length' in data or 'bootstrap' in data

    def test_incident_parent_edges_with_keys(self) -> None:
        """Test incident_parent_edges with keys=True."""
        net = DirectedPhyNetwork(
            edges=[
                (7, 5),
                (7, 6),
                (5, 4, 0), (5, 4, 1),  # Parallel edges
                (5, 8),
                (6, 4),
                (6, 9),
                (4, 2)
            ],
            taxa={2: "A", 8: "B", 9: "C"}
        )
        parent_edges = list(net.incident_parent_edges(4, keys=True))
        assert len(parent_edges) >= 2
        # Check keys are included
        for edge in parent_edges:
            assert len(edge) >= 3  # (u, v, key) or (u, v, key, data)

    def test_incident_parent_edges_with_keys_and_data(self) -> None:
        """Test incident_parent_edges with keys=True and data=True."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 7, 'v': 5},
                {'u': 7, 'v': 6},
                {'u': 5, 'v': 4, 'key': 0, 'branch_length': 0.5},
                {'u': 5, 'v': 4, 'key': 1, 'branch_length': 0.7},
                {'u': 5, 'v': 8},
                {'u': 6, 'v': 4},
                {'u': 6, 'v': 9},
                {'u': 4, 'v': 2}
            ],
            taxa={2: "A", 8: "B", 9: "C"}
        )
        parent_edges = list(net.incident_parent_edges(4, keys=True, data=True))
        assert len(parent_edges) >= 2
        for edge in parent_edges:
            assert len(edge) == 4  # (u, v, key, data)
            u, v, key, data = edge
            assert v == 4

    def test_incident_parent_edges_parallel(self) -> None:
        """Test incident_parent_edges with parallel edges."""
        net = DirectedPhyNetwork(
            edges=[
                (7, 5),
                (7, 6),
                (5, 4, 0), (5, 4, 1), (5, 4, 2),  # 3 parallel edges
                (5, 8),
                (6, 4),
                (6, 9),
                (4, 2)
            ],
            taxa={2: "A", 8: "B", 9: "C"}
        )
        parent_edges = list(net.incident_parent_edges(4, keys=True))
        # Should have 3 edges from 5, 1 from 6
        edges_from_5 = [e for e in parent_edges if e[0] == 5]
        assert len(edges_from_5) == 3


class TestIncidentChildEdges:
    """Test cases for incident_child_edges() method."""

    def test_incident_child_edges_leaf(self) -> None:
        """Test incident_child_edges for leaf (empty)."""
        net = DirectedPhyNetwork(edges=[(3, 1)], taxa={1: "A"})
        child_edges = list(net.incident_child_edges(1))
        assert len(child_edges) == 0

    def test_incident_child_edges_root(self) -> None:
        """Test incident_child_edges for root."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"})
        child_edges = list(net.incident_child_edges(3))
        assert len(child_edges) == 2
        assert (3, 1) in child_edges
        assert (3, 2) in child_edges

    def test_incident_child_edges_tree_node(self) -> None:
        """Test incident_child_edges for tree node."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            taxa={2: "A", 8: "B", 9: "C"}
        )
        child_edges = list(net.incident_child_edges(5))
        assert len(child_edges) == 2
        assert (5, 4) in child_edges
        assert (5, 8) in child_edges

    def test_incident_child_edges_hybrid(self) -> None:
        """Test incident_child_edges for hybrid node."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            taxa={2: "A", 8: "B", 9: "C"}
        )
        child_edges = list(net.incident_child_edges(4))
        assert len(child_edges) == 1
        assert (4, 2) in child_edges

    def test_incident_child_edges_with_data(self) -> None:
        """Test incident_child_edges with data=True."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 1, 'v': 2, 'branch_length': 0.5, 'bootstrap': 0.95},
                {'u': 1, 'v': 3, 'branch_length': 0.3, 'bootstrap': 0.87}
            ],
            taxa={2: "A", 3: "B"}
        )
        child_edges = list(net.incident_child_edges(1, data=True))
        assert len(child_edges) == 2
        for edge in child_edges:
            assert len(edge) == 3  # (u, v, data)
            u, v, data = edge
            assert u == 1

    def test_incident_child_edges_with_keys(self) -> None:
        """Test incident_child_edges with keys=True."""
        net = DirectedPhyNetwork(
            edges=[
                (4, 1),
                (1, 2, 0), (1, 2, 1),  # Parallel edges
                (1, 3),
                (2, 5)  # Make 2 a hybrid node
            ],
            taxa={3: "B", 5: "A"}
        )
        child_edges = list(net.incident_child_edges(1, keys=True))
        assert len(child_edges) >= 2
        for edge in child_edges:
            assert len(edge) >= 3  # (u, v, key) or more

    def test_incident_child_edges_parallel(self) -> None:
        """Test incident_child_edges with parallel edges."""
        net = DirectedPhyNetwork(
            edges=[
                (5, 4, 0), (5, 4, 1), (5, 4, 2),  # 3 parallel edges
                (5, 8),
                (4, 2)
            ],
            taxa={2: "A", 8: "B"}
        )
        child_edges = list(net.incident_child_edges(5, keys=True))
        # Should have 3 edges to 4, 1 to 8
        edges_to_4 = [e for e in child_edges if e[1] == 4]
        assert len(edges_to_4) == 3


class TestIncidentEdgesConsistency:
    """Test consistency between incident edges and parents/children."""

    def test_parent_edges_match_parents(self) -> None:
        """Test that incident_parent_edges matches parents()."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            taxa={2: "A", 8: "B", 9: "C"}
        )
        for node in net._graph.nodes:
            parents = set(net.parents(node))
            parent_edges = list(net.incident_parent_edges(node))
            parent_nodes_from_edges = {edge[0] for edge in parent_edges}
            assert parents == parent_nodes_from_edges

    def test_child_edges_match_children(self) -> None:
        """Test that incident_child_edges matches children()."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            taxa={2: "A", 8: "B", 9: "C"}
        )
        for node in net._graph.nodes:
            children = set(net.children(node))
            child_edges = list(net.incident_child_edges(node))
            child_nodes_from_edges = {edge[1] for edge in child_edges}
            assert children == child_nodes_from_edges

    def test_hybrid_node_incident_edges(self) -> None:
        """Test incident edges for hybrid node with gamma."""
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
        # Check parent edges (incoming to hybrid)
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
        
        # Check child edges (outgoing from hybrid)
        child_edges = list(net.incident_child_edges(4))
        assert len(child_edges) == 1
        assert (4, 2) in child_edges


class TestIncidentEdgesEdgeCases:
    """Test edge cases for incident edge methods."""

    def test_incident_edges_empty_network(self) -> None:
        """Test incident edges on empty network."""
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(edges=[])
        # Should handle gracefully (no nodes to test)

    def test_incident_edges_nonexistent_node(self) -> None:
        """Test incident edges for non-existent node."""
        net = DirectedPhyNetwork(edges=[(3, 1)], taxa={1: "A"})
        # Should return empty iterator
        parent_edges = list(net.incident_parent_edges(999))
        child_edges = list(net.incident_child_edges(999))
        assert len(parent_edges) == 0
        assert len(child_edges) == 0

    def test_incident_edges_many_parallel(self) -> None:
        """Test incident edges with many parallel edges."""
        edges = [(7, 5), (7, 6)]
        edges.extend([(5, 4, i) for i in range(10)])  # 10 parallel edges
        edges.append((5, 8))
        edges.append((6, 4))
        edges.append((6, 9))
        edges.append((4, 2))
        net = DirectedPhyNetwork(edges=edges, taxa={2: "A", 8: "B", 9: "C"})
        parent_edges = list(net.incident_parent_edges(4, keys=True))
        # Should have 10 edges from 5, 1 from 6
        edges_from_5 = [e for e in parent_edges if e[0] == 5]
        assert len(edges_from_5) == 10

    def test_incident_edges_all_combinations(self) -> None:
        """Test all combinations of keys and data parameters."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 7, 'v': 5},
                {'u': 7, 'v': 6},
                {'u': 5, 'v': 4, 'key': 0, 'branch_length': 0.5},
                {'u': 5, 'v': 4, 'key': 1, 'branch_length': 0.7},
                {'u': 5, 'v': 8},
                {'u': 6, 'v': 4},
                {'u': 6, 'v': 9},
                {'u': 4, 'v': 2}
            ],
            taxa={2: "A", 8: "B", 9: "C"}
        )
        # Test all 4 combinations
        edges_00 = list(net.incident_parent_edges(4, keys=False, data=False))
        edges_01 = list(net.incident_parent_edges(4, keys=False, data=True))
        edges_10 = list(net.incident_parent_edges(4, keys=True, data=False))
        edges_11 = list(net.incident_parent_edges(4, keys=True, data=True))
        
        # All should return edges
        assert len(edges_00) >= 2
        assert len(edges_01) >= 2
        assert len(edges_10) >= 2
        assert len(edges_11) >= 2
        
        # Check structure
        for edge in edges_00:
            assert len(edge) == 2  # (u, v)
        for edge in edges_01:
            assert len(edge) == 3  # (u, v, data)
        for edge in edges_10:
            assert len(edge) == 3  # (u, v, key)
        for edge in edges_11:
            assert len(edge) == 4  # (u, v, key, data)

