"""
Comprehensive tests for DirectedPhyNetwork graph operations.

This module tests all graph operation methods including:
- number_of_nodes()
- number_of_edges()
- has_edge()
- degree(), indegree(), outdegree()
- parents(), children(), neighbors()
- __contains__(), __iter__(), __len__()
- __repr__()
"""

import warnings

import pytest

from phylozoo.core.network import DirectedPhyNetwork


class TestNumberOfNodes:
    """Test cases for number_of_nodes() method."""

    def test_number_of_nodes_empty(self) -> None:
        """Test number_of_nodes on empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(edges=[])
        assert net.number_of_nodes() == 0

    def test_number_of_nodes_small(self) -> None:
        """Test number_of_nodes on small network."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        assert net.number_of_nodes() == 3

    def test_number_of_nodes_large(self) -> None:
        """Test number_of_nodes on large network."""
        edges = [(100, i) for i in range(1, 100)]
        nodes = [(i, {"label": f"Taxon{i}"}) for i in range(1, 100)]
        net = DirectedPhyNetwork(edges=edges, nodes=nodes)
        assert net.number_of_nodes() == 100


class TestNumberOfEdges:
    """Test cases for number_of_edges() method."""

    def test_number_of_edges_empty(self) -> None:
        """Test number_of_edges on empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(edges=[])
        assert net.number_of_edges() == 0

    def test_number_of_edges_simple(self) -> None:
        """Test number_of_edges on simple network."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        assert net.number_of_edges() == 2

    def test_number_of_edges_parallel(self) -> None:
        """Test number_of_edges with parallel edges."""
        # Use parallel edges to a hybrid node (valid structure)
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),  # Root to tree nodes
                (5, 4, 0), (5, 4, 1), (5, 4, 2),  # 3 parallel edges to hybrid
                (5, 8),  # Tree node 5 also has another child
                (6, 4),  # Tree node 6 also points to hybrid
                (6, 9),  # Tree node 6 also has another child
                (4, 2)   # Hybrid to leaf
            ],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        assert net.number_of_edges() == 9  # 2 (root->5, root->6) + 3 (parallel 5->4) + 1 (5->8) + 1 (6->4) + 1 (6->9) + 1 (4->2) = 9


class TestHasEdge:
    """Test cases for has_edge() method."""

    def test_has_edge_existing(self) -> None:
        """Test has_edge for existing edge."""
        net = DirectedPhyNetwork(edges=[(3, 1)], nodes=[(1, {'label': 'A'})])
        assert net.has_edge(3, 1) is True

    def test_has_edge_missing(self) -> None:
        """Test has_edge for missing edge."""
        net = DirectedPhyNetwork(edges=[(3, 1)], nodes=[(1, {'label': 'A'})])
        assert net.has_edge(3, 2) is False

    def test_has_edge_with_key(self) -> None:
        """Test has_edge with specific key."""
        # Use parallel edges to a hybrid node (valid structure)
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),  # Root to tree nodes
                (5, 4, 0), (5, 4, 1),  # Parallel edges to hybrid
                (5, 8),  # Tree node 5 also has another child
                (6, 4),  # Tree node 6 also points to hybrid
                (6, 9),  # Tree node 6 also has another child
                (4, 2)   # Hybrid to leaf
            ],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        assert net.has_edge(5, 4, key=0) is True
        assert net.has_edge(5, 4, key=1) is True
        assert net.has_edge(5, 4, key=2) is False

    def test_has_edge_parallel_without_key(self) -> None:
        """Test has_edge for parallel edges without key."""
        # Use parallel edges to a hybrid node (valid structure)
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),  # Root to tree nodes
                (5, 4, 0), (5, 4, 1),  # Parallel edges to hybrid
                (5, 8),  # Tree node 5 also has another child
                (6, 4),  # Tree node 6 also points to hybrid
                (6, 9),  # Tree node 6 also has another child
                (4, 2)   # Hybrid to leaf
            ],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        # Should return True if any edge exists (NetworkX behavior)
        assert net.has_edge(5, 4) is True


class TestDegree:
    """Test cases for degree(), indegree(), outdegree() methods."""

    def test_degree_root(self) -> None:
        """Test degree of root node."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        assert net.indegree(3) == 0
        assert net.outdegree(3) == 2
        assert net.degree(3) == 2

    def test_degree_leaf(self) -> None:
        """Test degree of leaf node."""
        net = DirectedPhyNetwork(edges=[(3, 1)], nodes=[(1, {'label': 'A'})])
        assert net.indegree(1) == 1
        assert net.outdegree(1) == 0
        assert net.degree(1) == 1

    def test_degree_tree_node(self) -> None:
        """Test degree of tree node."""
        net = DirectedPhyNetwork(
            edges=[(4, 3), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        assert net.indegree(3) == 1
        assert net.outdegree(3) == 2
        assert net.degree(3) == 3

    def test_degree_hybrid_node(self) -> None:
        """Test degree of hybrid node."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        assert net.indegree(4) == 2
        assert net.outdegree(4) == 1
        assert net.degree(4) == 3

    def test_degree_parallel_edges(self) -> None:
        """Test degree with parallel edges."""
        # Use parallel edges to a hybrid node (valid structure)
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),  # Root to tree nodes
                (5, 4, 0), (5, 4, 1), (5, 4, 2),  # 3 parallel edges to hybrid
                (5, 8),  # Tree node 5 also has another child
                (6, 4),  # Tree node 6 also points to hybrid
                (6, 9),  # Tree node 6 also has another child
                (4, 2)   # Hybrid to leaf
            ],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        # In-degree counts all parallel edges
        assert net.indegree(4) == 4  # 3 from 5, 1 from 6
        assert net.outdegree(5) == 4  # 3 to 4, 1 to 8


class TestParents:
    """Test cases for parents() method."""

    def test_parents_root(self) -> None:
        """Test parents of root (empty)."""
        net = DirectedPhyNetwork(edges=[(3, 1)], nodes=[(1, {'label': 'A'})])
        assert list(net.parents(3)) == []

    def test_parents_leaf(self) -> None:
        """Test parents of leaf."""
        net = DirectedPhyNetwork(edges=[(3, 1)], nodes=[(1, {'label': 'A'})])
        assert list(net.parents(1)) == [3]

    def test_parents_hybrid(self) -> None:
        """Test parents of hybrid node."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        parents = set(net.parents(4))
        assert parents == {5, 6}

    def test_parents_parallel_edges(self) -> None:
        """Test parents with parallel edges."""
        # Use parallel edges to a hybrid node (valid structure)
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),  # Root to tree nodes
                (5, 4, 0), (5, 4, 1),  # Parallel edges to hybrid
                (5, 8),  # Tree node 5 also has another child
                (6, 4),  # Tree node 6 also points to hybrid
                (6, 9),  # Tree node 6 also has another child
                (4, 2)   # Hybrid to leaf
            ],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        # Parents should only appear once (unique nodes)
        parents = set(net.parents(4))
        assert parents == {5, 6}


class TestChildren:
    """Test cases for children() method."""

    def test_children_leaf(self) -> None:
        """Test children of leaf (empty)."""
        net = DirectedPhyNetwork(edges=[(3, 1)], nodes=[(1, {'label': 'A'})])
        assert list(net.children(1)) == []

    def test_children_root(self) -> None:
        """Test children of root."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        children = set(net.children(3))
        assert children == {1, 2}

    def test_children_tree_node(self) -> None:
        """Test children of tree node."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        children = set(net.children(5))
        assert children == {4, 8}

    def test_children_hybrid(self) -> None:
        """Test children of hybrid node."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        children = list(net.children(4))
        assert children == [2]


class TestNeighbors:
    """Test cases for neighbors() method."""

    def test_neighbors_root(self) -> None:
        """Test neighbors of root."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        neighbors = set(net.neighbors(3))
        assert neighbors == {1, 2}  # Only children (no parents)

    def test_neighbors_leaf(self) -> None:
        """Test neighbors of leaf."""
        net = DirectedPhyNetwork(edges=[(3, 1)], nodes=[(1, {'label': 'A'})])
        neighbors = set(net.neighbors(1))
        assert neighbors == {3}  # Only parent (no children)

    def test_neighbors_hybrid(self) -> None:
        """Test neighbors of hybrid node."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        neighbors = set(net.neighbors(4))
        assert neighbors == {2, 5, 6}  # Parents + children

    def test_neighbors_consistency(self) -> None:
        """Test that neighbors = parents + children."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        for node in net._graph.nodes:
            neighbors = set(net.neighbors(node))
            parents_children = set(net.parents(node)) | set(net.children(node))
            assert neighbors == parents_children


class TestContains:
    """Test cases for __contains__() method."""

    def test_contains_existing_node(self) -> None:
        """Test __contains__ for existing node."""
        net = DirectedPhyNetwork(edges=[(3, 1)], nodes=[(1, {'label': 'A'})])
        assert 3 in net
        assert 1 in net

    def test_contains_missing_node(self) -> None:
        """Test __contains__ for missing node."""
        net = DirectedPhyNetwork(edges=[(3, 1)], nodes=[(1, {'label': 'A'})])
        assert 999 not in net

    def test_contains_all_nodes(self) -> None:
        """Test __contains__ for all nodes."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        for node in net._graph.nodes:
            assert node in net


class TestIter:
    """Test cases for __iter__() method."""

    def test_iter_all_nodes(self) -> None:
        """Test __iter__ includes all nodes."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        nodes = set(net)
        assert nodes == {1, 2, 3}

    def test_iter_empty(self) -> None:
        """Test __iter__ on empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(edges=[])
        nodes = list(net)
        assert len(nodes) == 0

    def test_iter_order(self) -> None:
        """Test that __iter__ order matches _graph.nodes."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            nodes=[(2, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
        )
        iter_nodes = list(net)
        graph_nodes = list(net._graph.nodes)
        assert set(iter_nodes) == set(graph_nodes)


class TestLen:
    """Test cases for __len__() method."""

    def test_len_consistency(self) -> None:
        """Test that __len__() matches number_of_nodes()."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        assert len(net) == net.number_of_nodes()

    def test_len_empty(self) -> None:
        """Test __len__ on empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(edges=[])
        assert len(net) == 0

    def test_len_large(self) -> None:
        """Test __len__ on large network."""
        edges = [(100, i) for i in range(1, 100)]
        nodes = [(i, {"label": f"Taxon{i}"}) for i in range(1, 100)]
        net = DirectedPhyNetwork(edges=edges, nodes=nodes)
        assert len(net) == 100


class TestRepr:
    """Test cases for __repr__() method."""

    def test_repr_format(self) -> None:
        """Test __repr__ format."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        repr_str = repr(net)
        assert "DirectedPhyNetwork" in repr_str
        assert "nodes=" in repr_str
        assert "edges=" in repr_str
        assert "taxa=" in repr_str

    def test_repr_truncation(self) -> None:
        """Test __repr__ truncates taxa list at 10."""
        edges = [(100, i) for i in range(1, 21)]
        nodes = [(i, {"label": f"Taxon{i}"}) for i in range(1, 21)]
        net = DirectedPhyNetwork(edges=edges, nodes=nodes)
        repr_str = repr(net)
        # Should show first 10 taxa and "..."
        assert "..." in repr_str
        # Should show taxa count (format may vary, check for number 20)
        assert "20" in repr_str or "taxa=20" in repr_str or "taxa:20" in repr_str

    def test_repr_empty(self) -> None:
        """Test __repr__ on empty network."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            net = DirectedPhyNetwork(edges=[])
        repr_str = repr(net)
        assert "nodes=0" in repr_str
        assert "edges=0" in repr_str
        assert "taxa=0" in repr_str

    def test_repr_small_taxa_list(self) -> None:
        """Test __repr__ with small taxa list (no truncation)."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        repr_str = repr(net)
        # Should show both taxa
        assert "A" in repr_str or "B" in repr_str

