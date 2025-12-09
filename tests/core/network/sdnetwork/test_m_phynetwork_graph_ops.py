"""
Comprehensive tests for MixedPhyNetwork graph operations.

This module tests all graph operation methods including:
- number_of_nodes()
- number_of_edges()
- has_edge()
- degree(), indegree(), outdegree(), undirected_degree()
- neighbors()
- __contains__(), __iter__(), __len__()
- __repr__()
"""

import warnings

import pytest

from phylozoo.core.network.sdnetwork import MixedPhyNetwork
from tests.core.network.sdnetwork.conftest import expect_mixed_network_warning


class TestNumberOfNodes:
    """Test cases for number_of_nodes() method."""

    def test_number_of_nodes_empty(self) -> None:
        """Test number_of_nodes on empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            # Empty networks don't raise validity warning (validation is skipped)
            net = MixedPhyNetwork(directed_edges=[], undirected_edges=[])
        assert net.number_of_nodes() == 0

    def test_number_of_nodes_small(self) -> None:
        """Test number_of_nodes on small network."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        assert net.number_of_nodes() == 4

    def test_number_of_nodes_large(self) -> None:
        """Test number_of_nodes on large network."""
        edges = [(100, i) for i in range(1, 100)]
        taxa = {i: f"Taxon{i}" for i in range(1, 100)}
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(undirected_edges=edges, taxa=taxa)
        assert net.number_of_nodes() == 100


class TestNumberOfEdges:
    """Test cases for number_of_edges() method."""

    def test_number_of_edges_empty(self) -> None:
        """Test number_of_edges on empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            # Empty networks don't raise validity warning (validation is skipped)
            net = MixedPhyNetwork(directed_edges=[], undirected_edges=[])
        assert net.number_of_edges() == 0

    def test_number_of_edges_simple(self) -> None:
        """Test number_of_edges on simple network."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        assert net.number_of_edges() == 3

    def test_number_of_edges_mixed(self) -> None:
        """Test number_of_edges with both directed and undirected."""
        # Nodes 5 and 6 need degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            taxa={2: "A", 8: "B", 9: "C", 10: "D", 11: "E"}
            )
        assert net.number_of_edges() == 7

    def test_number_of_edges_parallel(self) -> None:
        """Test number_of_edges with parallel edges."""
        # Nodes 5 and 6 need degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[
            (5, 4, 0), (5, 4, 1), (5, 4, 2),  # 3 parallel edges
            (6, 4)
            ],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            taxa={2: "A", 8: "B", 9: "C", 10: "D", 11: "E"}
            )
        assert net.number_of_edges() == 9


class TestHasEdge:
    """Test cases for has_edge() method."""

    def test_has_edge_undirected_existing(self) -> None:
        """Test has_edge for existing undirected edge."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        assert net.has_edge(3, 1) is True
        assert net.has_edge(1, 3) is True  # Undirected is symmetric

    def test_has_edge_directed_existing(self) -> None:
        """Test has_edge for existing directed edge (hybrid edge)."""
        # Hybrid node 4 with directed edges from 3 and 5
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(3, 4), (5, 4)],
            undirected_edges=[(4, 1), (3, 2), (3, 6), (5, 7), (5, 8)],
            taxa={1: "A", 2: "B", 6: "C", 7: "D", 8: "E"}
            )
        assert net.has_edge(3, 4) is True
        assert net.has_edge(4, 3) is False  # Directed is not symmetric

    def test_has_edge_missing(self) -> None:
        """Test has_edge for missing edge."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        assert net.has_edge(1, 2) is False  # No direct edge between 1 and 2

    def test_has_edge_with_key(self) -> None:
        """Test has_edge with specific key."""
        # Parallel edges between internal nodes
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[
            (3, 4, 0), (3, 4, 1),  # Parallel edges
            (3, 1), (3, 2),  # Additional edges from 3
            (4, 5), (4, 6)   # Additional edges from 4
            ],
            taxa={1: "A", 2: "B", 5: "C", 6: "D"}
            )
        assert net.has_edge(3, 4, key=0) is True
        assert net.has_edge(3, 4, key=1) is True
        assert net.has_edge(3, 4, key=2) is False

    def test_has_edge_directed_parameter(self) -> None:
        """Test has_edge with directed parameter."""
        # Hybrid node 5 with directed edges from 3 and 6
        # Node 2 needs degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(3, 5), (6, 5)],
            undirected_edges=[(5, 1), (2, 3), (2, 4), (2, 11), (3, 7), (3, 8), (6, 9), (6, 10)],
            taxa={1: "A", 4: "B", 7: "C", 8: "D", 9: "E", 10: "F", 11: "G"}
            )
        assert net.has_edge(3, 5, directed=True) is True
        assert net.has_edge(3, 5, directed=False) is False
        assert net.has_edge(2, 4, directed=False) is True
        assert net.has_edge(2, 4, directed=True) is False


class TestDegree:
    """Test cases for degree() method."""

    def test_degree_leaf(self) -> None:
        """Test degree of leaf node."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        assert net.degree(1) == 1

    def test_degree_internal_node(self) -> None:
        """Test degree of internal node."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        assert net.degree(3) == 3

    def test_degree_hybrid_node(self) -> None:
        """Test degree of hybrid node."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            taxa={2: "A", 8: "B", 9: "C", 10: "D", 11: "E"}
            )
        # Hybrid 4: indegree 2, undirected degree 1, total = 3
        assert net.degree(4) == 3

    def test_degree_mixed_edges(self) -> None:
        """Test degree with both directed and undirected edges."""
        # Node 4 needs to be hybrid (indegree >= 2)
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(3, 4), (5, 4)],  # Node 4 is hybrid
            undirected_edges=[(3, 1), (3, 2), (4, 6), (5, 7), (5, 8)],
            taxa={1: "A", 2: "B", 6: "C", 7: "D", 8: "E"}
            )
        # Node 3: outdegree 1, undirected degree 2, total = 3
        assert net.degree(3) == 3


class TestIndegree:
    """Test cases for indegree() method."""

    def test_indegree_leaf(self) -> None:
        """Test indegree of leaf node."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        assert net.indegree(1) == 0  # Undirected edges don't count as indegree

    def test_indegree_hybrid(self) -> None:
        """Test indegree of hybrid node."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            taxa={2: "A", 8: "B", 9: "C", 10: "D", 11: "E"}
            )
        assert net.indegree(4) == 2

    def test_indegree_tree_node(self) -> None:
        """Test indegree of tree node."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        assert net.indegree(3) == 0


class TestOutdegree:
    """Test cases for outdegree() method."""

    def test_outdegree_leaf(self) -> None:
        """Test outdegree of leaf node."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        assert net.outdegree(1) == 0  # Leaves have no outgoing directed edges

    def test_outdegree_tree_node(self) -> None:
        """Test outdegree of tree node."""
        # Nodes 4 and 5 need to be hybrids (indegree >= 2)
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(3, 4), (6, 4), (3, 5), (7, 5)],  # 4 and 5 are hybrids
            undirected_edges=[(3, 1), (4, 8), (5, 9), (6, 10), (6, 11), (7, 12), (7, 13)],
            taxa={1: "A", 8: "B", 9: "C", 10: "D", 11: "E", 12: "F", 13: "G"}
            )
        assert net.outdegree(3) == 2


class TestUndirectedDegree:
    """Test cases for undirected_degree() method."""

    def test_undirected_degree_leaf(self) -> None:
        """Test undirected_degree of leaf node."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        assert net.undirected_degree(1) == 1

    def test_undirected_degree_internal(self) -> None:
        """Test undirected_degree of internal node."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        assert net.undirected_degree(3) == 3

    def test_undirected_degree_hybrid(self) -> None:
        """Test undirected_degree of hybrid node."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(4, 2), (5, 8), (5, 10), (6, 9), (6, 11)],
            taxa={2: "A", 8: "B", 9: "C", 10: "D", 11: "E"}
            )
        assert net.undirected_degree(4) == 1


class TestNeighbors:
    """Test cases for neighbors() method."""

    def test_neighbors_undirected(self) -> None:
        """Test neighbors with undirected edges."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        neighbors = set(net.neighbors(3))
        assert neighbors == {1, 2, 4}

    def test_neighbors_directed(self) -> None:
        """Test neighbors with directed edges (hybrid edges)."""
        # Nodes 1 and 2 need to be hybrids (indegree >= 2)
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(3, 1), (4, 1), (3, 2), (5, 2)],
            undirected_edges=[(1, 6), (2, 7), (3, 8), (3, 9), (4, 10), (4, 11), (5, 12), (5, 13)],
            taxa={6: "A", 7: "B", 8: "C", 9: "D", 10: "E", 11: "F", 12: "G", 13: "H"}
            )
        neighbors = set(net.neighbors(3))
        assert 1 in neighbors
        assert 2 in neighbors

    def test_neighbors_mixed(self) -> None:
        """Test neighbors with both directed and undirected edges."""
        # Node 1 is hybrid (indegree >= 2), total_degree must be 3 (only 1 outgoing)
        # Nodes 3 and 4 need degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            directed_edges=[(3, 1), (4, 1)],
            undirected_edges=[(1, 8), (3, 2), (3, 5), (4, 6), (4, 7)],
            taxa={2: "A", 5: "B", 6: "C", 7: "D", 8: "E"}
            )
        neighbors = set(net.neighbors(3))
        assert 1 in neighbors
        assert 2 in neighbors
        assert 5 in neighbors


class TestSpecialMethods:
    """Test cases for special methods (__contains__, __iter__, __len__)."""

    def test_contains_existing_node(self) -> None:
        """Test __contains__ for existing node."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        assert 3 in net
        assert 1 in net

    def test_contains_nonexistent_node(self) -> None:
        """Test __contains__ for non-existent node."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        assert 999 not in net

    def test_iter_nodes(self) -> None:
        """Test __iter__ returns all nodes."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        nodes = set(net)
        assert nodes == {1, 2, 3, 4}

    def test_len_nodes(self) -> None:
        """Test __len__ returns number of nodes."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        assert len(net) == 4

    def test_repr(self) -> None:
        """Test __repr__ method."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
            )
        repr_str = repr(net)
        assert isinstance(repr_str, str)
        assert "MixedPhyNetwork" in repr_str

