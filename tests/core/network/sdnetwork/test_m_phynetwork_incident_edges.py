"""
Comprehensive tests for MixedPhyNetwork incident edge methods.

This module tests all incident edge methods including:
- incident_parent_edges()
- incident_child_edges()
- incident_undirected_edges()
- Parallel edges handling
- Edge attributes in incident edges
"""

import pytest

from phylozoo.core.network.sdnetwork import MixedPhyNetwork
from tests.core.network.sdnetwork.conftest import expect_mixed_network_warning


class TestIncidentParentEdges:
    """Test cases for incident_parent_edges() method."""

    def test_incident_parent_edges_leaf(self) -> None:
        """Test incident_parent_edges for leaf."""
        # Node 3 needs degree >= 3
        # Node 4 is hybrid (indegree >= 2), node 1 is leaf (indegree 0)
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(3, 4), (5, 4)],  # Hybrid 4
            undirected_edges=[(4, 1), (3, 2), (3, 6), (5, 7), (5, 8)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (6, {'label': 'C'}), (7, {'label': 'D'}), (8, {'label': 'E'})]
            )
        # Node 1 is a leaf with indegree 0, so no parent edges
        parent_edges = list(net.incident_parent_edges(1))
        assert len(parent_edges) == 0

    def test_incident_parent_edges_hybrid(self) -> None:
        """Test incident_parent_edges for hybrid node."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        parent_edges = list(net.incident_parent_edges(4))
        assert len(parent_edges) == 2
        assert (5, 4) in parent_edges
        assert (6, 4) in parent_edges

    def test_incident_parent_edges_tree_node(self) -> None:
        """Test incident_parent_edges for tree node."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        # Tree nodes have indegree 0, so no parent edges
        parent_edges = list(net.incident_parent_edges(3))
        assert len(parent_edges) == 0

    def test_incident_parent_edges_with_data(self) -> None:
        """Test incident_parent_edges with data=True."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[
            {'u': 5, 'v': 4, 'branch_length': 0.5, 'bootstrap': 0.95},
            {'u': 6, 'v': 4, 'branch_length': 0.3, 'bootstrap': 0.87}
            ],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
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
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[
            (5, 4, 0), (5, 4, 1),  # Parallel edges
            (6, 4)
            ],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (10, {'label': 'D'}), (11, {'label': 'E'})]
            )
        parent_edges = list(net.incident_parent_edges(4, keys=True))
        assert len(parent_edges) >= 2
        # Check keys are included
        for edge in parent_edges:
            assert len(edge) >= 3  # (u, v, key) or (u, v, key, data)


class TestIncidentChildEdges:
    """Test cases for incident_child_edges() method."""

    def test_incident_child_edges_tree_node(self) -> None:
        """Test incident_child_edges for tree node."""
        # Node 3 is tree node, node 4 is hybrid (indegree >= 2)
        # Node 4 needs total_degree = indegree + 1 = 3 (only 1 outgoing)
        # incident_child_edges only returns directed outgoing edges, not undirected
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(3, 4), (5, 4)],
            undirected_edges=[(3, 1), (3, 2), (4, 6), (5, 7), (5, 8)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (6, {'label': 'C'}), (7, {'label': 'D'}), (8, {'label': 'E'})]
            )
        child_edges = list(net.incident_child_edges(3))
        assert len(child_edges) == 1  # Only (3, 4) directed edge
        assert (3, 4) in child_edges

    def test_incident_child_edges_leaf(self) -> None:
        """Test incident_child_edges for leaf."""
        # Node 3 needs degree >= 3
        # Node 1 is a leaf (indegree 0, degree 1)
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        child_edges = list(net.incident_child_edges(1))
        assert len(child_edges) == 0  # Leaves have no outgoing directed edges

    def test_incident_child_edges_with_data(self) -> None:
        """Test incident_child_edges with data=True."""
        # Node 3 is tree node, node 4 is hybrid (indegree >= 2)
        # Node 4 needs total_degree = indegree + 1 = 3 (only 1 outgoing)
        # incident_child_edges only returns directed outgoing edges
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[
            {'u': 3, 'v': 4, 'branch_length': 0.5},
            {'u': 5, 'v': 4, 'branch_length': 0.3}
            ],
            undirected_edges=[(3, 1), (3, 2), (4, 6), (5, 7), (5, 8)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (6, {'label': 'C'}), (7, {'label': 'D'}), (8, {'label': 'E'})]
            )
        child_edges = list(net.incident_child_edges(3, data=True))
        assert len(child_edges) == 1  # Only (3, 4) directed edge
        for edge in child_edges:
            assert len(edge) == 3  # (u, v, data)
            u, v, data = edge
            assert u == 3
            assert v == 4
            assert 'branch_length' in data


class TestIncidentUndirectedEdges:
    """Test cases for incident_undirected_edges() method."""

    def test_incident_undirected_edges_leaf(self) -> None:
        """Test incident_undirected_edges for leaf."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        undir_edges = list(net.incident_undirected_edges(1))
        assert len(undir_edges) == 1
        assert (3, 1) in undir_edges or (1, 3) in undir_edges

    def test_incident_undirected_edges_internal(self) -> None:
        """Test incident_undirected_edges for internal node."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        undir_edges = list(net.incident_undirected_edges(3))
        assert len(undir_edges) == 3

    def test_incident_undirected_edges_with_data(self) -> None:
        """Test incident_undirected_edges with data=True."""
        # Node 3 needs degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[
            {'u': 3, 'v': 1, 'branch_length': 0.5},
            {'u': 3, 'v': 2, 'branch_length': 0.3},
            (3, 4)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        undir_edges = list(net.incident_undirected_edges(3, data=True))
        assert len(undir_edges) == 3
        for edge in undir_edges:
            assert len(edge) == 3  # (u, v, data)

    def test_incident_undirected_edges_with_keys(self) -> None:
        """Test incident_undirected_edges with keys=True."""
        # Node 3 needs degree >= 3
        # Use parallel edges between internal nodes, not leaves (leaves have degree == 1)
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[
            (3, 4, 0), (3, 4, 1),  # Parallel edges between internal nodes
            (3, 1), (3, 2), (4, 5), (4, 6)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'})]
            )
        undir_edges = list(net.incident_undirected_edges(3, keys=True))
        assert len(undir_edges) >= 3  # (3, 4) with 2 keys, (3, 1), (3, 2)
        for edge in undir_edges:
            assert len(edge) >= 3  # (u, v, key) or (u, v, key, data)


class TestIncidentEdgesConsistency:
    """Test consistency between different incident edge methods."""

    def test_parent_child_edges_consistency(self) -> None:
        """Test that parent/child edges are consistent."""
        # Node 3 needs degree >= 3
        # Use undirected edges only to avoid indegree constraint issues
        # Actually, we need directed edges for this test, so use a hybrid node structure
        # Node 4 is hybrid (indegree >= 2), node 1 is leaf (indegree 0)
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(3, 4), (5, 4)],  # Hybrid 4
            undirected_edges=[(4, 1), (3, 2), (3, 6), (5, 7), (5, 8)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (6, {'label': 'C'}), (7, {'label': 'D'}), (8, {'label': 'E'})]
            )
        # Node 4 has parent edges from 3 and 5
        parent_edges_4 = list(net.incident_parent_edges(4))
        child_edges_3 = list(net.incident_child_edges(3))
        assert (3, 4) in parent_edges_4
        assert (3, 4) in child_edges_3

    def test_all_incident_edges_cover_degree(self) -> None:
        """Test that all incident edges cover total degree."""
        # Node 3 needs degree >= 3
        # Node 1 should have indegree 0 (not 1) to satisfy constraint
        # Use undirected edges only to avoid indegree issues
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        # Node 3: 0 incoming directed, 0 outgoing directed, 3 undirected
        parent_edges = list(net.incident_parent_edges(3))
        child_edges = list(net.incident_child_edges(3))
        undir_edges = list(net.incident_undirected_edges(3))
        
        # Total should match degree
        total_incident = len(parent_edges) + len(child_edges) + len(undir_edges)
        assert total_incident == net.degree(3)

