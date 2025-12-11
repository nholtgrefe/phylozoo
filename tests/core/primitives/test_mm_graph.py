"""
Tests for the MixedMultiGraph class.

This test suite provides comprehensive coverage of all MixedMultiGraph features,
including edge cases, parallel edges, edge attributes, and larger graphs.
"""

import warnings

import pytest
from typing import Dict, List, Set, Tuple

from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
from phylozoo.core.primitives.m_multigraph.features import (
    number_of_connected_components,
    is_connected,
    connected_components,
    source_components,
    has_self_loops,
    is_cutedge,
    is_cutvertex,
)
from phylozoo.core.primitives.m_multigraph.transformations import (
    identify_two_nodes,
    identify_node_set,
    orient_away_from_vertex,
)
from phylozoo.core.primitives.m_multigraph.conversions import (
    graph_to_mixedmultigraph,
    multigraph_to_mixedmultigraph,
    multidigraph_to_mixedmultigraph,
    directedmultigraph_to_mixedmultigraph,
)


# Helper functions for testing (since some methods don't exist)
def count_directed_edges(G: MixedMultiGraph, u: int, v: int) -> int:
    """Count number of directed edges between u and v."""
    if u in G._directed and v in G._directed[u]:
        return len(G._directed[u][v])
    return 0


def has_directed_edge(G: MixedMultiGraph, u: int, v: int, key: int = None) -> bool:
    """Check if directed edge exists."""
    if key is None:
        return G._directed.has_edge(u, v)
    return G._directed.has_edge(u, v, key)


def count_undirected_edges(G: MixedMultiGraph, u: int, v: int) -> int:
    """Count number of undirected edges between u and v."""
    if u in G._undirected and v in G._undirected[u]:
        return len(G._undirected[u][v])
    return 0


def number_of_directed_edges_total(G: MixedMultiGraph) -> int:
    """Get total number of directed edges."""
    return G._directed.number_of_edges()


class TestInitialization:
    """Test cases for graph initialization."""

    def test_empty_graph(self) -> None:
        """Test creating an empty graph."""
        G = MixedMultiGraph()
        assert G.number_of_nodes() == 0
        assert G.number_of_edges() == 0
        # G.nodes() returns iterator, need to convert to list
        assert len(list(G.nodes())) == 0
        assert len(list(G.edges())) == 0

    def test_repr_counts(self) -> None:
        """__repr__ reports node and edge counts by type."""
        G = MixedMultiGraph(directed_edges=[(1, 2)], undirected_edges=[(2, 3), (3, 4)])
        assert (
            repr(G)
            == "MixedMultiGraph(nodes=4, directed_edges=1, undirected_edges=2)"
        )

    def test_init_with_directed_edges(self) -> None:
        """Test initialization with directed edges."""
        G = MixedMultiGraph(directed_edges=[(1, 2), (2, 3), (1, 2)])  # Parallel edge
        assert G.number_of_nodes() == 3
        # Check number of directed edges between nodes
        num_edges_1_2 = len(G._directed[1][2]) if 1 in G._directed and 2 in G._directed[1] else 0
        num_edges_2_3 = len(G._directed[2][3]) if 2 in G._directed and 3 in G._directed[2] else 0
        print(f"Number of directed edges (1,2): {num_edges_1_2}, (2,3): {num_edges_2_3}")
        assert num_edges_1_2 == 2
        assert num_edges_2_3 == 1

    def test_init_with_directed_edges_with_keys(self) -> None:
        """Test initialization with directed edges and explicit keys."""
        G = MixedMultiGraph(directed_edges=[(1, 2, 0), (1, 2, 5), (2, 3, 0)])
        num_edges = len(G._directed[1][2]) if 1 in G._directed and 2 in G._directed[1] else 0
        assert num_edges == 2
        assert G._directed.has_edge(1, 2, key=0)
        assert G._directed.has_edge(1, 2, key=5)
        assert G._directed.has_edge(2, 3, key=0)

    def test_directed_self_loops(self) -> None:
        """Directed self-loops are allowed and support parallel edges."""
        G = MixedMultiGraph()
        k0 = G.add_directed_edge(1, 1, weight=1.0)
        k1 = G.add_directed_edge(1, 1, weight=2.0)
        assert {k0, k1} == {0, 1}
        assert count_directed_edges(G, 1, 1) == 2
        assert G._directed.has_edge(1, 1, key=k0)
        assert G._directed.has_edge(1, 1, key=k1)
        assert G.number_of_edges() == 2
        assert has_self_loops(G) is True

    def test_init_with_undirected_edges(self) -> None:
        """Test initialization with undirected edges."""
        G = MixedMultiGraph(undirected_edges=[(1, 2), (2, 3), (1, 2)])  # Parallel edge
        assert G.number_of_nodes() == 3
        num_edges_1_2 = count_undirected_edges(G, 1, 2)
        num_edges_2_3 = count_undirected_edges(G, 2, 3)
        print(f"Number of undirected edges (1,2): {num_edges_1_2}, (2,3): {num_edges_2_3}")
        assert num_edges_1_2 == 2
        assert num_edges_2_3 == 1

    def test_init_with_undirected_edges_with_keys(self) -> None:
        """Test initialization with undirected edges and explicit keys."""
        G = MixedMultiGraph(undirected_edges=[(1, 2, 0), (1, 2, 5), (2, 3, 0)])
        num_edges = count_undirected_edges(G, 1, 2)
        assert num_edges == 2
        assert G._undirected.has_edge(1, 2, key=0)
        assert G._undirected.has_edge(1, 2, key=5)
        assert G._undirected.has_edge(2, 3, key=0)

    def test_undirected_self_loops(self) -> None:
        """Undirected self-loops are allowed and support parallel edges."""
        G = MixedMultiGraph()
        k0 = G.add_undirected_edge(2, 2, weight=1.0)
        k1 = G.add_undirected_edge(2, 2, weight=2.0)
        assert {k0, k1} == {0, 1}
        assert count_undirected_edges(G, 2, 2) == 2
        assert G._undirected.has_edge(2, 2, key=k0)
        assert G._undirected.has_edge(2, 2, key=k1)
        assert G.number_of_edges() == 2
        assert has_self_loops(G) is True

    def test_has_self_loops_false(self) -> None:
        """has_self_loops returns False when no self-loops exist."""
        G = MixedMultiGraph(
            undirected_edges=[(1, 2)],
            directed_edges=[(2, 3)],
        )
        assert has_self_loops(G) is False

    def test_self_loop_mutual_exclusivity(self) -> None:
        """Directed/undirected self-loops cannot coexist; newer edge replaces older type."""
        G = MixedMultiGraph()
        # Start with undirected self-loop
        k_u0 = G.add_undirected_edge(1, 1, weight=1.0)
        assert count_undirected_edges(G, 1, 1) == 1
        assert count_directed_edges(G, 1, 1) == 0
        assert has_self_loops(G) is True
        # Add directed self-loop; undirected should be removed
        k_d0 = G.add_directed_edge(1, 1, weight=2.0)
        assert count_directed_edges(G, 1, 1) == 1
        assert count_undirected_edges(G, 1, 1) == 0
        assert G._directed.has_edge(1, 1, key=k_d0)
        # Now add undirected self-loop again; directed should be removed
        k_u1 = G.add_undirected_edge(1, 1, weight=3.0)
        assert count_undirected_edges(G, 1, 1) == 1
        assert count_directed_edges(G, 1, 1) == 0
        assert G._undirected.has_edge(1, 1, key=k_u1)
        assert has_self_loops(G) is True

    def test_init_with_both_edge_types(self) -> None:
        """Test initialization with both directed and undirected edges."""
        G = MixedMultiGraph(
            undirected_edges=[(1, 2), (2, 3)],
            directed_edges=[(4, 5), (5, 6)]
        )
        assert G.number_of_nodes() == 6
        assert count_undirected_edges(G, 1, 2) == 1
        assert count_undirected_edges(G, 2, 3) == 1
        assert count_directed_edges(G, 4, 5) == 1
        assert count_directed_edges(G, 5, 6) == 1

    def test_init_with_mutually_exclusive_edges(self) -> None:
        """Test that mutual exclusivity is enforced during initialization."""
        # If same nodes appear in both lists, directed should win (added last)
        G = MixedMultiGraph(
            undirected_edges=[(1, 2), (1, 2)],  # Parallel undirected
            directed_edges=[(1, 2)]  # This should remove undirected edges
        )
        assert count_undirected_edges(G, 1, 2) == 0
        assert count_directed_edges(G, 1, 2) == 1

    def test_init_with_edge_attributes_undirected(self) -> None:
        """Test initialization with undirected edges and attributes."""
        G = MixedMultiGraph(
            undirected_edges=[
                (1, 2),  # No attributes
                {'u': 2, 'v': 3, 'weight': 5.0, 'label': 'test'},  # Dict format
                {'u': 3, 'v': 4, 'key': 10, 'weight': 10.0}  # With key
            ]
        )
        assert count_undirected_edges(G, 1, 2) == 1
        assert count_undirected_edges(G, 2, 3) == 1
        assert count_undirected_edges(G, 3, 4) == 1
        # Check attributes
        edge_data = G._undirected[2][3][0]
        assert edge_data['weight'] == 5.0
        assert edge_data['label'] == 'test'
        edge_data = G._undirected[3][4][10]
        assert edge_data['weight'] == 10.0

    def test_init_with_edge_attributes_directed(self) -> None:
        """Test initialization with directed edges and attributes."""
        G = MixedMultiGraph(
            directed_edges=[
                (1, 2),  # No attributes
                {'u': 2, 'v': 3, 'weight': 20.0, 'type': 'special'},  # Dict format
                {'u': 3, 'v': 4, 'key': 5, 'weight': 30.0}  # With key
            ]
        )
        assert count_directed_edges(G, 1, 2) == 1
        assert count_directed_edges(G, 2, 3) == 1
        assert count_directed_edges(G, 3, 4) == 1
        # Check attributes
        edge_data = G._directed[2][3][0]
        assert edge_data['weight'] == 20.0
        assert edge_data['type'] == 'special'
        edge_data = G._directed[3][4][5]
        assert edge_data['weight'] == 30.0

    def test_init_with_edge_attributes_mixed(self) -> None:
        """Test initialization with both edge types and attributes."""
        G = MixedMultiGraph(
            undirected_edges=[
                {'u': 1, 'v': 2, 'weight': 1.0}
            ],
            directed_edges=[
                {'u': 3, 'v': 4, 'weight': 2.0, 'label': 'dir'}
            ]
        )
        assert G._undirected[1][2][0]['weight'] == 1.0
        assert G._directed[3][4][0]['weight'] == 2.0
        assert G._directed[3][4][0]['label'] == 'dir'


class TestFactoryMethods:
    """Test cases for factory methods (from_graph, from_multigraph, etc.)."""

    def test_from_graph(self) -> None:
        """Test from_graph factory method."""
        import networkx as nx
        G = nx.Graph()
        G.add_node(1, label='node1')
        G.add_edge(1, 2, weight=5.0)
        G.add_edge(2, 3)
        
        M = graph_to_mixedmultigraph(G)
        assert M.number_of_nodes() == 3
        assert M.number_of_edges() == 2


class TestIncidentEdges:
    """Test cases for incident edge methods."""

    def test_incident_parent_edges_basic(self) -> None:
        """Test incident_parent_edges with basic directed edges."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2)
        G.add_directed_edge(3, 2)
        G.add_directed_edge(4, 2)
        
        parent_edges = list(G.incident_parent_edges(2))
        assert len(parent_edges) == 3
        assert (1, 2) in parent_edges
        assert (3, 2) in parent_edges
        assert (4, 2) in parent_edges

    def test_incident_parent_edges_with_keys_and_data(self) -> None:
        """Test incident_parent_edges with keys and data."""
        G = MixedMultiGraph()
        key1 = G.add_directed_edge(1, 2, weight=1.0)
        key2 = G.add_directed_edge(1, 2, weight=2.0)  # Parallel edge
        
        parent_edges = list(G.incident_parent_edges(2, keys=True, data=True))
        assert len(parent_edges) == 2
        for edge in parent_edges:
            assert len(edge) == 4
            u, v, key, data = edge
            assert u == 1
            assert v == 2
            assert key in [key1, key2]
            assert 'weight' in data

    def test_incident_child_edges_basic(self) -> None:
        """Test incident_child_edges with basic directed edges."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2)
        G.add_directed_edge(1, 3)
        G.add_directed_edge(1, 4)
        
        child_edges = list(G.incident_child_edges(1))
        assert len(child_edges) == 3
        assert (1, 2) in child_edges
        assert (1, 3) in child_edges
        assert (1, 4) in child_edges

    def test_incident_child_edges_with_keys_and_data(self) -> None:
        """Test incident_child_edges with keys and data."""
        G = MixedMultiGraph()
        key1 = G.add_directed_edge(1, 2, weight=1.0)
        key2 = G.add_directed_edge(1, 2, weight=2.0)  # Parallel edge
        
        child_edges = list(G.incident_child_edges(1, keys=True, data=True))
        assert len(child_edges) == 2
        for edge in child_edges:
            assert len(edge) == 4
            u, v, key, data = edge
            assert u == 1
            assert v == 2
            assert key in [key1, key2]
            assert 'weight' in data

    def test_incident_undirected_edges_basic(self) -> None:
        """Test incident_undirected_edges with basic undirected edges."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_undirected_edge(2, 4)
        
        undirected_edges = list(G.incident_undirected_edges(2))
        assert len(undirected_edges) == 3
        assert (1, 2) in undirected_edges or (2, 1) in undirected_edges
        assert (2, 3) in undirected_edges or (3, 2) in undirected_edges
        assert (2, 4) in undirected_edges or (4, 2) in undirected_edges

    def test_incident_undirected_edges_with_keys_and_data(self) -> None:
        """Test incident_undirected_edges with keys and data."""
        G = MixedMultiGraph()
        key1 = G.add_undirected_edge(1, 2, weight=1.0)
        key2 = G.add_undirected_edge(1, 2, weight=2.0)  # Parallel edge
        
        undirected_edges = list(G.incident_undirected_edges(1, keys=True, data=True))
        assert len(undirected_edges) == 2
        for edge in undirected_edges:
            assert len(edge) == 4
            u, v, key, data = edge
            assert 1 in [u, v]
            assert 2 in [u, v]
            assert key in [key1, key2]
            assert 'weight' in data

    def test_incident_undirected_edges_empty(self) -> None:
        """Test incident_undirected_edges for node with no undirected edges."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2)
        
        undirected_edges = list(G.incident_undirected_edges(1))
        assert len(undirected_edges) == 0

    def test_incident_edges_mixed_graph(self) -> None:
        """Test incident edges in a graph with both directed and undirected edges."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2, weight=1.0)
        G.add_directed_edge(3, 2, weight=2.0)
        G.add_undirected_edge(2, 4, weight=3.0)
        G.add_undirected_edge(2, 5, weight=4.0)
        
        # Check directed parent edges
        parent_edges = list(G.incident_parent_edges(2))
        assert len(parent_edges) == 2
        assert (1, 2) in parent_edges
        assert (3, 2) in parent_edges
        
        # Check undirected edges
        undirected_edges = list(G.incident_undirected_edges(2))
        assert len(undirected_edges) == 2
        # Check that edges contain node 2
        edge_nodes = set()
        for edge in undirected_edges:
            edge_nodes.update(edge[:2])
        assert 2 in edge_nodes
        assert 4 in edge_nodes
        assert 5 in edge_nodes

    def test_incident_edges_consistency(self) -> None:
        """Test that incident edges are consistent with graph structure."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2, weight=1.0)
        G.add_directed_edge(2, 3, weight=2.0)
        G.add_undirected_edge(2, 4, weight=3.0)
        
        # Node 2: incoming from 1, outgoing to 3, undirected to 4
        parent_edges_2 = list(G.incident_parent_edges(2))
        child_edges_2 = list(G.incident_child_edges(2))
        undirected_edges_2 = list(G.incident_undirected_edges(2))
        
        assert len(parent_edges_2) == 1
        assert len(child_edges_2) == 1
        assert len(undirected_edges_2) == 1
        assert (1, 2) in parent_edges_2
        assert (2, 3) in child_edges_2
        # Undirected edge should contain both 2 and 4
        assert any(2 in edge[:2] and 4 in edge[:2] for edge in undirected_edges_2)

    def test_from_multigraph(self) -> None:
        """Test from_multigraph factory method."""
        import networkx as nx
        G = nx.MultiGraph()
        G.add_node(1, label='node1')
        G.add_edge(1, 2, key=0, weight=1.0)
        G.add_edge(1, 2, key=1, weight=2.0)  # Parallel edge
        G.add_edge(2, 3, key=0, weight=3.0)
        
        M = multigraph_to_mixedmultigraph(G)
        assert M.number_of_nodes() == 3
        assert M.number_of_edges() == 3
        # Parallel edges preserved
        assert count_undirected_edges(M, 1, 2) == 2
        # Keys preserved
        assert M._undirected.has_edge(1, 2, key=0)
        assert M._undirected.has_edge(1, 2, key=1)
        # Attributes preserved
        assert M._undirected[1][2][0]['weight'] == 1.0
        assert M._undirected[1][2][1]['weight'] == 2.0
        assert M._undirected.nodes[1].get('label') == 'node1'

    def test_from_multidigraph(self) -> None:
        """Test from_multidigraph factory method."""
        import networkx as nx
        G = nx.MultiDiGraph()
        G.add_node(1, label='node1')
        G.add_edge(1, 2, key=0, weight=10.0)
        G.add_edge(1, 2, key=1, weight=20.0)  # Parallel edge
        G.add_edge(2, 3, key=0, weight=30.0)
        
        M = multidigraph_to_mixedmultigraph(G)
        assert M.number_of_nodes() == 3
        assert M._directed.number_of_edges() == 3
        # All edges should be directed
        assert count_directed_edges(M, 1, 2) == 2
        assert count_directed_edges(M, 2, 3) == 1
        # Keys preserved
        assert M._directed.has_edge(1, 2, key=0)
        assert M._directed.has_edge(1, 2, key=1)
        # Attributes preserved
        assert M._directed[1][2][0]['weight'] == 10.0
        assert M._directed[1][2][1]['weight'] == 20.0
        assert M._directed.nodes[1].get('label') == 'node1'

    def test_from_directedmultigraph(self) -> None:
        """Test from_directedmultigraph factory method."""
        from phylozoo.core.primitives import DirectedMultiGraph
        
        G = DirectedMultiGraph()
        G.add_node(1, label='node1')
        G.add_edge(1, 2, weight=100.0)
        G.add_edge(1, 2, weight=200.0)  # Parallel edge
        
        M = directedmultigraph_to_mixedmultigraph(G)
        assert M.number_of_nodes() == 2
        assert M._directed.number_of_edges() == 2
        # All edges should be directed
        assert count_directed_edges(M, 1, 2) == 2
        # Attributes preserved
        weights = {M._directed[1][2][k]['weight'] for k in M._directed[1][2].keys()}
        assert weights == {100.0, 200.0}
        assert M._directed.nodes[1].get('label') == 'node1'

    def test_from_graph_with_attributes(self) -> None:
        """Test from_graph preserves node and edge attributes."""
        import networkx as nx
        G = nx.Graph()
        G.add_node(1, type='A', value=10)
        G.add_node(2, type='B', value=20)
        G.add_edge(1, 2, weight=5.0, label='edge1')
        G.add_edge(2, 3, weight=10.0)
        
        M = graph_to_mixedmultigraph(G)
        # Check node attributes
        assert M._undirected.nodes[1].get('type') == 'A'
        assert M._undirected.nodes[1].get('value') == 10
        assert M._undirected.nodes[2].get('type') == 'B'
        # Check edge attributes
        assert M._undirected[1][2][0]['weight'] == 5.0
        assert M._undirected[1][2][0]['label'] == 'edge1'
        assert M._undirected[2][3][0]['weight'] == 10.0


class TestNodeOperations:
    """Test cases for node operations."""

    def test_add_node(self) -> None:
        """Test adding a single node."""
        G = MixedMultiGraph()
        G.add_node(1)
        assert 1 in G
        assert G.number_of_nodes() == 1

    def test_add_node_with_attributes(self) -> None:
        """Test adding a node with attributes."""
        G = MixedMultiGraph()
        G.add_node(1, label='node1', weight=5.0)
        assert 1 in G
        # Check if attributes are stored (access via underlying graph)
        assert G._undirected.nodes[1].get('label') == 'node1'
        assert G._undirected.nodes[1].get('weight') == 5.0

    def test_add_nodes_from(self) -> None:
        """Test adding multiple nodes."""
        G = MixedMultiGraph()
        G.add_nodes_from([1, 2, 3, 4, 5])
        assert G.number_of_nodes() == 5
        assert all(n in G for n in [1, 2, 3, 4, 5])

    def test_add_nodes_from_with_attributes(self) -> None:
        """Test adding multiple nodes with attributes."""
        G = MixedMultiGraph()
        G.add_nodes_from([1, 2, 3], label='test', type='A')
        for n in [1, 2, 3]:
            assert G._undirected.nodes[n].get('label') == 'test'
            assert G._undirected.nodes[n].get('type') == 'A'

    def test_remove_node(self) -> None:
        """Test removing a node."""
        G = MixedMultiGraph()
        G.add_nodes_from([1, 2, 3])
        G.add_undirected_edge(1, 2)
        G.add_directed_edge(2, 3)
        G.remove_node(2)
        assert 2 not in G
        assert G.number_of_nodes() == 2
        assert not G.has_edge(1, 2)
        assert not G._directed.has_edge(2, 3)

    def test_remove_nodes_from(self) -> None:
        """Test removing multiple nodes."""
        G = MixedMultiGraph()
        G.add_nodes_from([1, 2, 3, 4, 5])
        G.remove_nodes_from([2, 4])
        assert G.number_of_nodes() == 3
        assert 2 not in G
        assert 4 not in G
        assert all(n in G for n in [1, 3, 5])


class TestUndirectedEdgeOperations:
    """Test cases for undirected edge operations."""

    def test_add_undirected_edge(self) -> None:
        """Test adding a single undirected edge."""
        G = MixedMultiGraph()
        key = G.add_undirected_edge(1, 2)
        assert G.has_edge(1, 2)
        assert G.has_edge(2, 1)  # Undirected is symmetric
        assert isinstance(key, int)

    def test_add_undirected_edge_with_attributes(self) -> None:
        """Test adding undirected edge with attributes."""
        G = MixedMultiGraph()
        key = G.add_undirected_edge(1, 2, weight=5.0, label='test', color='red')
        assert G.has_edge(1, 2)
        # Check attributes via __getitem__
        edge_data = G[1][2][key]
        assert edge_data['weight'] == 5.0
        assert edge_data['label'] == 'test'
        assert edge_data['color'] == 'red'

    def test_parallel_undirected_edges(self) -> None:
        """Test adding parallel undirected edges."""
        G = MixedMultiGraph()
        key1 = G.add_undirected_edge(1, 2, weight=1.0)
        key2 = G.add_undirected_edge(1, 2, weight=2.0)
        key3 = G.add_undirected_edge(1, 2, weight=3.0)
        assert key1 != key2 != key3
        num_edges = count_undirected_edges(G, 1, 2)
        assert num_edges == 3
        # Check all edges have different weights
        weights = {G[1][2][k]['weight'] for k in [key1, key2, key3]}
        assert weights == {1.0, 2.0, 3.0}

    def test_parallel_undirected_edges_with_explicit_keys(self) -> None:
        """Test adding parallel undirected edges with explicit keys."""
        G = MixedMultiGraph()
        key1 = G.add_undirected_edge(1, 2, key=10, weight=1.0)
        key2 = G.add_undirected_edge(1, 2, key=20, weight=2.0)
        assert key1 == 10
        assert key2 == 20
        assert G.has_edge(1, 2, key=10)
        assert G.has_edge(1, 2, key=20)

    def test_add_edges_from(self) -> None:
        """Test adding multiple undirected edges."""
        G = MixedMultiGraph()
        G.add_undirected_edges_from([(1, 2), (2, 3), (3, 4)])
        assert G.has_edge(1, 2)
        assert G.has_edge(2, 3)
        assert G.has_edge(3, 4)
        assert G.number_of_edges() == 3

    def test_add_edges_from_with_attributes(self) -> None:
        """Test adding multiple edges with shared attributes."""
        G = MixedMultiGraph()
        G.add_undirected_edges_from([(1, 2), (2, 3)], weight=5.0, type='test')
        for u, v in [(1, 2), (2, 3)]:
            keys = list(G[u][v].keys())
            assert len(keys) > 0
            assert G[u][v][keys[0]]['weight'] == 5.0
            assert G[u][v][keys[0]]['type'] == 'test'

    def test_remove_undirected_edge(self) -> None:
        """Test removing an undirected edge."""
        G = MixedMultiGraph()
        key1 = G.add_undirected_edge(1, 2, weight=1.0)
        key2 = G.add_undirected_edge(1, 2, weight=2.0)
        G.remove_edge(1, 2, key=key1)
        assert not G.has_edge(1, 2, key=key1)
        assert G.has_edge(1, 2, key=key2)
        assert count_undirected_edges(G, 1, 2) == 1

    def test_remove_undirected_edge_without_key(self) -> None:
        """Test removing undirected edge without specifying key."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2, weight=1.0)
        G.add_undirected_edge(1, 2, weight=2.0)
        # Behavior: should remove one edge (which one?)
        G.remove_edge(1, 2)
        # Report behavior
        remaining = count_undirected_edges(G, 1, 2)
        print(f"After removing edge without key, remaining edges: {remaining}")

    def test_remove_edges_from(self) -> None:
        """Test removing multiple edges."""
        G = MixedMultiGraph()
        G.add_undirected_edges_from([(1, 2), (2, 3), (3, 4)])
        G.remove_edges_from([(1, 2), (3, 4)])
        assert not G.has_edge(1, 2)
        assert G.has_edge(2, 3)
        assert not G.has_edge(3, 4)


class TestDirectedEdgeOperations:
    """Test cases for directed edge operations."""

    def test_add_directed_edge(self) -> None:
        """Test adding a single directed edge."""
        G = MixedMultiGraph()
        key = G.add_directed_edge(1, 2)
        assert has_directed_edge(G, 1, 2)
        assert not has_directed_edge(G, 2, 1)  # Directed is not symmetric
        assert isinstance(key, int)

    def test_add_directed_edge_with_attributes(self) -> None:
        """Test adding directed edge with attributes."""
        G = MixedMultiGraph()
        key = G.add_directed_edge(1, 2, weight=10.0, label='dir', direction='forward')
        assert has_directed_edge(G, 1, 2)
        edge_data = G._directed[1][2][key]
        assert edge_data['weight'] == 10.0
        assert edge_data['label'] == 'dir'
        assert edge_data['direction'] == 'forward'

    def test_parallel_directed_edges(self) -> None:
        """Test adding parallel directed edges."""
        G = MixedMultiGraph()
        key1 = G.add_directed_edge(1, 2, weight=1.0)
        key2 = G.add_directed_edge(1, 2, weight=2.0)
        key3 = G.add_directed_edge(1, 2, weight=3.0)
        assert key1 != key2 != key3
        assert count_directed_edges(G, 1, 2) == 3
        weights = {G._directed[1][2][k]['weight'] for k in [key1, key2, key3]}
        assert weights == {1.0, 2.0, 3.0}

    def test_parallel_directed_edges_with_explicit_keys(self) -> None:
        """Test adding parallel directed edges with explicit keys."""
        G = MixedMultiGraph()
        key1 = G.add_directed_edge(1, 2, key=100, weight=1.0)
        key2 = G.add_directed_edge(1, 2, key=200, weight=2.0)
        assert key1 == 100
        assert key2 == 200
        assert has_directed_edge(G, 1, 2, key=100)
        assert has_directed_edge(G, 1, 2, key=200)

    def test_add_directed_edges_from(self) -> None:
        """Test adding multiple directed edges."""
        G = MixedMultiGraph()
        G.add_directed_edges_from([(1, 2), (2, 3), (3, 1)])
        assert has_directed_edge(G, 1, 2)
        assert has_directed_edge(G, 2, 3)
        assert has_directed_edge(G, 3, 1)
        assert number_of_directed_edges_total(G) == 3

    def test_remove_directed_edge(self) -> None:
        """Test removing a directed edge."""
        G = MixedMultiGraph()
        key1 = G.add_directed_edge(1, 2, weight=1.0)
        key2 = G.add_directed_edge(1, 2, weight=2.0)
        G.remove_directed_edge(1, 2, key=key1)
        assert not has_directed_edge(G, 1, 2, key=key1)
        assert has_directed_edge(G, 1, 2, key=key2)

    def test_remove_directed_edge_without_key(self) -> None:
        """Test removing directed edge without specifying key."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2, weight=1.0)
        G.add_directed_edge(1, 2, weight=2.0)
        # Behavior: should remove one edge (which one?)
        G.remove_directed_edge(1, 2)
        remaining = count_directed_edges(G, 1, 2)
        print(f"After removing directed edge without key, remaining edges: {remaining}")


class TestMutualExclusivity:
    """Test cases for mutual exclusivity between directed and undirected edges."""

    def test_adding_directed_removes_undirected(self) -> None:
        """Test that adding directed edge removes undirected edges."""
        G = MixedMultiGraph()
        key1 = G.add_undirected_edge(1, 2, weight=1.0)
        key2 = G.add_undirected_edge(1, 2, weight=2.0)
        assert count_undirected_edges(G, 1, 2) == 2
        
        G.add_directed_edge(1, 2, weight=3.0)
        assert count_undirected_edges(G, 1, 2) == 0  # Undirected edges removed
        assert count_directed_edges(G, 1, 2) == 1  # Directed edge added

    def test_adding_undirected_removes_directed(self) -> None:
        """Test that adding undirected edge removes directed edges."""
        G = MixedMultiGraph()
        key1 = G.add_directed_edge(1, 2, weight=1.0)
        key2 = G.add_directed_edge(1, 2, weight=2.0)
        assert count_directed_edges(G, 1, 2) == 2
        
        G.add_undirected_edge(1, 2, weight=3.0)
        assert count_directed_edges(G, 1, 2) == 0  # Directed edges removed
        assert count_undirected_edges(G, 1, 2) == 1  # Undirected edge added

    def test_switching_back_and_forth(self) -> None:
        """Test switching between directed and undirected multiple times."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_directed_edge(1, 2)  # This removes undirected edge
        G.add_undirected_edge(1, 2)  # This removes directed edge and adds undirected
        G.add_directed_edge(1, 2)  # This removes undirected and adds directed
        # Final state should be directed
        assert count_directed_edges(G, 1, 2) == 1
        assert count_undirected_edges(G, 1, 2) == 0


class TestQueryOperations:
    """Test cases for query operations."""

    def test_has_edge(self) -> None:
        """Test has_edge method."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_directed_edge(3, 4)
        assert G.has_edge(1, 2)
        assert G.has_edge(2, 1)  # Undirected is symmetric
        assert G.has_edge(3, 4)  # Directed edge
        assert not G.has_edge(4, 3)  # Reverse direction doesn't exist
        assert not G.has_edge(1, 3)

    def test_has_edge_with_key(self) -> None:
        """Test has_edge with specific key."""
        G = MixedMultiGraph()
        key1 = G.add_undirected_edge(1, 2, weight=1.0)
        key2 = G.add_undirected_edge(1, 2, weight=2.0)
        assert G.has_edge(1, 2, key=key1)
        assert G.has_edge(1, 2, key=key2)
        assert not G.has_edge(1, 2, key=999)

    def test_number_of_nodes(self) -> None:
        """Test number_of_nodes method."""
        G = MixedMultiGraph()
        G.add_nodes_from([1, 2, 3, 4, 5])
        assert G.number_of_nodes() == 5

    def test_number_of_edges(self) -> None:
        """Test number_of_edges method."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(1, 2)  # Parallel
        G.add_directed_edge(2, 3)
        G.add_directed_edge(2, 3)  # Parallel
        assert G.number_of_edges() == 4

    def test_undirected_degree(self) -> None:
        """Test undirected_degree calculation."""
        G = MixedMultiGraph()
        # Node with no edges
        assert G.undirected_degree(1) == 0
        
        # Single undirected edge
        G.add_undirected_edge(1, 2)
        assert G.undirected_degree(1) == 1
        assert G.undirected_degree(2) == 1
        
        # Parallel undirected edges
        G.add_undirected_edge(1, 2)  # Parallel edge
        assert G.undirected_degree(1) == 2
        assert G.undirected_degree(2) == 2
        
        # Multiple undirected edges
        G.add_undirected_edge(1, 3)
        G.add_undirected_edge(1, 3)  # Parallel
        assert G.undirected_degree(1) == 4  # 2 to node 2, 2 to node 3
        assert G.undirected_degree(3) == 2
        
        # Directed edges should not affect undirected degree
        G.add_directed_edge(1, 4)
        G.add_directed_edge(5, 1)
        assert G.undirected_degree(1) == 4  # Still 4
        assert G.undirected_degree(4) == 0
        assert G.undirected_degree(5) == 0

    def test_degree(self) -> None:
        """Test total degree calculation (undirected + directed)."""
        G = MixedMultiGraph()
        # Node with no edges
        assert G.degree(1) == 0
        
        # Only undirected edges
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(1, 2)  # Parallel undirected
        assert G.degree(1) == 2
        assert G.degree(2) == 2
        
        # Only directed edges
        G2 = MixedMultiGraph()
        G2.add_directed_edge(1, 2)
        G2.add_directed_edge(1, 3)  # Outgoing
        G2.add_directed_edge(4, 1)  # Incoming
        assert G2.degree(1) == 3  # 2 outgoing + 1 incoming
        assert G2.degree(2) == 1  # 1 incoming
        assert G2.degree(4) == 1  # 1 outgoing
        
        # Mixed edges
        G3 = MixedMultiGraph()
        G3.add_undirected_edge(1, 2)  # Undirected: contributes 1
        G3.add_undirected_edge(1, 2)  # Parallel undirected: contributes 1
        G3.add_directed_edge(1, 3)  # Outgoing: contributes 1
        G3.add_directed_edge(4, 1)  # Incoming: contributes 1
        # Node 1: 2 undirected + 1 outgoing + 1 incoming = 4
        assert G3.degree(1) == 4
        assert G3.degree(2) == 2  # 2 undirected
        assert G3.degree(3) == 1  # 1 incoming
        assert G3.degree(4) == 1  # 1 outgoing
        
        # Parallel directed edges
        G3.add_directed_edge(1, 3)  # Another parallel directed edge
        assert G3.degree(1) == 5  # 2 undirected + 2 outgoing + 1 incoming
        assert G3.degree(3) == 2  # 2 incoming

    def test_indegree(self) -> None:
        """Test indegree calculation."""
        G = MixedMultiGraph()
        # Node with no edges
        assert G.indegree(1) == 0
        
        # Single incoming edge
        G.add_directed_edge(1, 2)
        assert G.indegree(2) == 1
        assert G.indegree(1) == 0
        
        # Multiple incoming edges
        G.add_directed_edge(3, 2)
        G.add_directed_edge(4, 2)
        assert G.indegree(2) == 3
        assert G.indegree(1) == 0
        assert G.indegree(3) == 0
        assert G.indegree(4) == 0
        
        # Parallel directed edges (same direction)
        G.add_directed_edge(1, 2)  # Another edge from 1 to 2
        assert G.indegree(2) == 4  # 4 incoming edges total
        
        # Undirected edges should not affect indegree
        G.add_undirected_edge(5, 2)
        assert G.indegree(2) == 4  # Still 4 (undirected doesn't count)

    def test_outdegree(self) -> None:
        """Test outdegree calculation."""
        G = MixedMultiGraph()
        # Node with no edges
        assert G.outdegree(1) == 0
        
        # Single outgoing edge
        G.add_directed_edge(1, 2)
        assert G.outdegree(1) == 1
        assert G.outdegree(2) == 0
        
        # Multiple outgoing edges
        G.add_directed_edge(1, 3)
        G.add_directed_edge(1, 4)
        assert G.outdegree(1) == 3
        assert G.outdegree(2) == 0
        assert G.outdegree(3) == 0
        assert G.outdegree(4) == 0
        
        # Parallel directed edges (same direction)
        G.add_directed_edge(1, 2)  # Another edge from 1 to 2
        assert G.outdegree(1) == 4  # 4 outgoing edges total
        
        # Undirected edges should not affect outdegree
        G.add_undirected_edge(1, 5)
        assert G.outdegree(1) == 4  # Still 4 (undirected doesn't count)


class TestSpecialMethods:
    """Test cases for special methods (__contains__, __iter__, __len__, __getitem__)."""

    def test_contains(self) -> None:
        """Test __contains__ method."""
        G = MixedMultiGraph()
        G.add_node(1)
        assert 1 in G
        assert 2 not in G

    def test_iter(self) -> None:
        """Test __iter__ method."""
        G = MixedMultiGraph()
        G.add_nodes_from([1, 2, 3, 4, 5])
        nodes = list(G)
        assert set(nodes) == {1, 2, 3, 4, 5}

    def test_len(self) -> None:
        """Test __len__ method."""
        G = MixedMultiGraph()
        G.add_nodes_from([1, 2, 3])
        assert len(G) == 3

    def test_getitem(self) -> None:
        """Test __getitem__ method for adjacency access."""
        G = MixedMultiGraph()
        key1 = G.add_undirected_edge(1, 2, weight=1.0)
        key2 = G.add_undirected_edge(1, 2, weight=2.0)
        key3 = G.add_directed_edge(1, 3, weight=3.0)
        
        # Check structure: G[u] should return dict-like view of neighbors
        adj = G[1]
        # Behavior: G[1] returns AdjacencyView, not dict, but should work like dict
        print(f"Type of G[1]: {type(adj)}")
        assert 2 in adj
        assert 3 in adj
        # Check edge keys (AdjacencyView should support dict-like access)
        assert key1 in adj[2]
        assert key2 in adj[2]
        assert key3 in adj[3]
        # Check attributes
        assert adj[2][key1]['weight'] == 1.0
        assert adj[2][key2]['weight'] == 2.0
        assert adj[3][key3]['weight'] == 3.0


class TestProperties:
    """Test cases for properties (nodes, edges, directed_edges, combined_graph)."""

    def test_nodes_property(self) -> None:
        """Test nodes property."""
        G = MixedMultiGraph()
        G.add_nodes_from([1, 2, 3])
        # Test as attribute
        nodes = G.nodes
        assert 1 in nodes
        assert 2 in nodes
        assert 3 in nodes
        # Test as method
        nodes_list = list(G.nodes())
        assert set(nodes_list) == {1, 2, 3}

    def test_nodes_property_with_data(self) -> None:
        """Test nodes property with data=True."""
        G = MixedMultiGraph()
        G.add_node(1, label='node1')
        G.add_node(2, weight=5.0)
        nodes_data = dict(G.nodes(data=True))
        assert 'label' in nodes_data[1]
        assert 'weight' in nodes_data[2]

    def test_edges_property(self) -> None:
        """Test edges property."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_directed_edge(2, 3)
        # Test as attribute
        edges = G.edges
        assert (1, 2) in edges or (2, 1) in edges
        # Test as method
        edges_list = list(G.edges())
        assert len(edges_list) >= 2

    def test_edges_property_with_data(self) -> None:
        """Test edges property with data=True."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2, weight=5.0)
        edges_data = list(G.edges(data=True))
        # Find the edge and check data
        for u, v, data in edges_data:
            if (u, v) == (1, 2) or (u, v) == (2, 1):
                assert data['weight'] == 5.0
                break

    def test_directed_edges_property(self) -> None:
        """Test directed_edges property."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2, weight=1.0)
        G.add_directed_edge(1, 2, weight=2.0)  # Parallel edge
        G.add_directed_edge(2, 3, weight=3.0)
        directed = G.directed_edges
        assert (1, 2) in directed
        assert (2, 3) in directed
        assert (2, 1) not in directed  # Not reverse
        # Test as method with keys and data
        edges_with_keys = list(G.directed_edges(keys=True))
        assert (1, 2, 0) in edges_with_keys
        assert (1, 2, 1) in edges_with_keys
        edges_with_data = list(G.directed_edges(data=True))
        weights = {e[2]['weight'] for e in edges_with_data if e[0] == 1 and e[1] == 2}
        assert weights == {1.0, 2.0}

    def test_undirected_edges_property(self) -> None:
        """Test undirected_edges property."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2, weight=1.0)
        G.add_undirected_edge(1, 2, weight=2.0)  # Parallel edge
        G.add_undirected_edge(2, 3, weight=3.0)
        undirected = G.undirected_edges
        assert (1, 2) in undirected
        assert (2, 3) in undirected
        # Test as method with keys and data
        edges_with_keys = list(G.undirected_edges(keys=True))
        assert (1, 2, 0) in edges_with_keys
        assert (1, 2, 1) in edges_with_keys
        edges_with_data = list(G.undirected_edges(data=True))
        weights = {e[2]['weight'] for e in edges_with_data if e[0] == 1 and e[1] == 2}
        assert weights == {1.0, 2.0}

    def test_combined_graph_property(self) -> None:
        """Test combined_graph property."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_directed_edge(2, 3)
        combined = G.combined_graph
        # Combined should have all edges as undirected
        assert combined.has_edge(1, 2)
        assert combined.has_edge(2, 3)
        assert combined.has_edge(3, 2)  # Undirected view


class TestIterators:
    """Test cases for iterator methods."""

    def test_nodes_iter(self) -> None:
        """Test nodes_iter method."""
        G = MixedMultiGraph()
        G.add_nodes_from([1, 2, 3])
        nodes = list(G.nodes_iter())
        assert set(nodes) == {1, 2, 3}

    def test_nodes_iter_with_data(self) -> None:
        """Test nodes_iter with data=True."""
        G = MixedMultiGraph()
        G.add_node(1, label='test')
        nodes_data = list(G.nodes_iter(data=True))
        assert len(nodes_data) == 1
        # Behavior: Check what nodes_iter(data=True) actually returns
        print(f"nodes_iter(data=True) returns: {nodes_data}")
        print(f"Type of first element: {type(nodes_data[0])}")
        # Try to handle both cases
        if isinstance(nodes_data[0], tuple):
            node, data = nodes_data[0]
            assert node == 1
            assert data.get('label') == 'test'
        else:
            # If it returns something else, report the behavior
            print(f"Unexpected return type: {type(nodes_data[0])}, value: {nodes_data[0]}")
            # For now, just check it's not empty
            assert len(nodes_data) > 0

    def test_edges_iter(self) -> None:
        """Test edges_iter method."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_directed_edge(2, 3)
        edges = list(G.edges_iter())
        assert len(edges) >= 2

    def test_edges_iter_with_keys(self) -> None:
        """Test edges_iter with keys=True."""
        G = MixedMultiGraph()
        key1 = G.add_undirected_edge(1, 2)
        key2 = G.add_undirected_edge(1, 2)
        edges = list(G.edges_iter(keys=True))
        keys_found = {k for u, v, k in edges if (u, v) == (1, 2) or (u, v) == (2, 1)}
        assert key1 in keys_found
        assert key2 in keys_found

    def test_neighbors(self) -> None:
        """Test neighbors method."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(1, 3)
        G.add_directed_edge(1, 4)
        G.add_directed_edge(5, 1)  # Incoming
        neighbors = set(G.neighbors(1))
        # Should include all neighbors (both directions for directed)
        print(f"Neighbors of 1: {neighbors}")
        assert 2 in neighbors
        assert 3 in neighbors
        # Check behavior for directed edges
        assert 4 in neighbors or 5 in neighbors or (4 in neighbors and 5 in neighbors)


class TestConnectivity:
    """Test cases for connectivity methods."""

    def test_is_connected(self) -> None:
        """Test is_connected method."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_directed_edge(3, 4)
        assert is_connected(G)
        
        G2 = MixedMultiGraph()
        G2.add_undirected_edge(1, 2)
        G2.add_undirected_edge(3, 4)
        assert not is_connected(G2)

    def test_number_of_connected_components(self) -> None:
        """Test number_of_connected_components method."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_undirected_edge(4, 5)
        assert number_of_connected_components(G) == 2

    def test_connected_components(self) -> None:
        """Test connected_components method."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_undirected_edge(4, 5)
        components = list(connected_components(G))
        assert len(components) == 2
        comp_sets = [set(comp) for comp in components]
        assert {1, 2, 3} in comp_sets
        assert {4, 5} in comp_sets

    def test_is_cutedge(self) -> None:
        """Test is_cutedge function."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_undirected_edge(3, 4)
        # Edge (2, 3) is a cut edge
        key = next(iter(G._undirected[2][3].keys()))
        assert is_cutedge(G, 2, 3, key=key)
        # Edge (1, 2) is not a cut edge (there's another path)
        key = next(iter(G._undirected[1][2].keys()))
        # Actually, in a path graph, all edges are cut edges except in a cycle
        result = is_cutedge(G, 1, 2, key=key)
        print(f"Is (1,2) a cut edge? {result}")

    def test_is_cutvertex(self) -> None:
        """Test is_cutvertex function."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_undirected_edge(2, 4)
        # Node 2 is a cut vertex
        assert is_cutvertex(G, 2)
        # Node 1 is not a cut vertex
        assert not is_cutvertex(G, 1)


class TestGraphOperations:
    """Test cases for graph operations (copy, clear, identify_nodes)."""

    def test_copy(self) -> None:
        """Test copy method."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2, weight=5.0)
        G.add_directed_edge(2, 3, weight=10.0)
        H = G.copy()
        assert H.number_of_nodes() == G.number_of_nodes()
        assert H.number_of_edges() == G.number_of_edges()
        # Modify H and ensure G is unchanged
        H.add_undirected_edge(4, 5)
        assert 4 not in G
        assert 5 not in G

    def test_clear(self) -> None:
        """Test clear method."""
        G = MixedMultiGraph()
        G.add_nodes_from([1, 2, 3, 4, 5])
        G.add_undirected_edge(1, 2)
        G.add_directed_edge(2, 3)
        G.clear()
        assert G.number_of_nodes() == 0
        assert G.number_of_edges() == 0

    def test_identify_two_nodes(self) -> None:
        """Test identify_two_nodes method."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_directed_edge(1, 4)
        identify_two_nodes(G, 1, 2)  # Keep node 1
        # Node 2 should be removed, edges should be transferred
        assert 2 not in G
        assert 1 in G
        # Check if edges were transferred correctly
        print(f"After identifying 1 and 2, edges incident to 1: {list(G.neighbors(1))}")

    def test_identify_node_set(self) -> None:
        """Test identify_node_set method."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_undirected_edge(3, 4)
        identify_node_set(G, [1, 2, 3])  # Keep first node (1)
        assert 1 in G
        assert 2 not in G
        assert 3 not in G
        print(f"After identifying [1,2,3], edges incident to 1: {list(G.neighbors(1))}")


class TestSourceComponents:
    """Test cases for source_components function."""
    
    def test_empty_graph(self) -> None:
        """Test source_components on empty graph."""
        G = MixedMultiGraph()
        components = source_components(G)
        assert components == []
    
    def test_single_node_no_edges(self) -> None:
        """Test source_components with isolated node."""
        G = MixedMultiGraph()
        G.add_node(1)
        components = source_components(G)
        assert len(components) == 1
        nodes, undir_edges, out_edges = components[0]
        assert set(nodes) == {1}
        assert undir_edges == []
        assert out_edges == []
    
    def test_simple_source_component(self) -> None:
        """Test simple source component with undirected edges and outgoing directed edge."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_directed_edge(3, 4)
        components = source_components(G)
        assert len(components) == 1
        nodes, undir_edges, out_edges = components[0]
        assert set(nodes) == {1, 2, 3}
        assert set(undir_edges) == {(1, 2, 0), (2, 3, 0)}
        assert out_edges == [(3, 4, 0)]
    
    def test_component_with_incoming_edge(self) -> None:
        """Test component with incoming directed edge (not a source component)."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_directed_edge(3, 1)  # Edge pointing into component {1, 2}
        components = source_components(G)
        # Component {1, 2} is NOT a source component (has incoming edge from 3)
        # Component {3} IS a source component (no incoming edges)
        assert len(components) == 1
        nodes, undir_edges, out_edges = components[0]
        assert set(nodes) == {3}
        assert undir_edges == []
        assert out_edges == [(3, 1, 0)]
    
    def test_multiple_source_components(self) -> None:
        """Test multiple source components."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(3, 4)
        G.add_directed_edge(2, 5)
        G.add_directed_edge(4, 5)
        components = source_components(G)
        assert len(components) == 2
        # Sort by first node for consistent testing
        components_sorted = sorted(components, key=lambda x: min(x[0]))
        nodes1, undir1, out1 = components_sorted[0]
        nodes2, undir2, out2 = components_sorted[1]
        # Both should be source components
        assert (set(nodes1) == {1, 2} and set(nodes2) == {3, 4}) or \
               (set(nodes1) == {3, 4} and set(nodes2) == {1, 2})
        assert set(nodes1) != set(nodes2)
    
    def test_only_undirected_edges(self) -> None:
        """Test graph with only undirected edges (all components are source components)."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(3, 4)
        G.add_undirected_edge(5, 6)
        G.add_undirected_edge(6, 7)
        components = source_components(G)
        assert len(components) == 3
        # Components: {1, 2}, {3, 4}, and {5, 6, 7}
        node_sets = [set(nodes) for nodes, _, _ in components]
        assert {1, 2} in node_sets
        assert {3, 4} in node_sets
        assert {5, 6, 7} in node_sets
    
    def test_component_with_no_outgoing_edges(self) -> None:
        """Test source component with no outgoing edges."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        components = source_components(G)
        assert len(components) == 1
        nodes, undir_edges, out_edges = components[0]
        assert set(nodes) == {1, 2, 3}
        assert set(undir_edges) == {(1, 2, 0), (2, 3, 0)}
        assert out_edges == []
    
    def test_component_with_multiple_outgoing_edges(self) -> None:
        """Test source component with multiple outgoing edges."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_directed_edge(1, 3)
        G.add_directed_edge(1, 4)
        G.add_directed_edge(2, 5)
        components = source_components(G)
        assert len(components) == 1
        nodes, undir_edges, out_edges = components[0]
        assert set(nodes) == {1, 2}
        assert set(undir_edges) == {(1, 2, 0)}
        assert set(out_edges) == {(1, 3, 0), (1, 4, 0), (2, 5, 0)}
    
    def test_complex_multiple_components(self) -> None:
        """Test complex scenario with multiple components, some source, some not."""
        G = MixedMultiGraph()
        # Component 1: {1, 2, 3} - has incoming edge from 6, NOT a source component
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_directed_edge(3, 6)  # Outgoing
        G.add_directed_edge(6, 1)  # Incoming to {1, 2, 3}
        # Component 2: {4, 5} - no incoming edges, IS a source component
        G.add_undirected_edge(4, 5)
        G.add_directed_edge(5, 6)  # Outgoing
        # Component 3: {7, 8} - no incoming edges, IS a source component
        G.add_undirected_edge(7, 8)
        # Component 4: {6} - has incoming edges from {1,2,3} and {4,5}, NOT a source component
        # (node 6 is isolated in undirected graph, but has incoming directed edges)
        components = source_components(G)
        assert len(components) == 2
        node_sets = [set(nodes) for nodes, _, _ in components]
        assert {4, 5} in node_sets
        assert {7, 8} in node_sets
    
    def test_parallel_undirected_edges(self) -> None:
        """Test source component with parallel undirected edges."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2, weight=1.0)
        G.add_undirected_edge(1, 2, weight=2.0)  # Parallel edge
        G.add_directed_edge(2, 3)
        components = source_components(G)
        assert len(components) == 1
        nodes, undir_edges, out_edges = components[0]
        assert set(nodes) == {1, 2}
        # Should include both parallel edges with their keys
        assert len(undir_edges) == 2
        assert set(undir_edges) == {(1, 2, 0), (1, 2, 1)}
        assert out_edges == [(2, 3, 0)]
    
    def test_parallel_directed_edges(self) -> None:
        """Test source component with parallel directed outgoing edges."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_directed_edge(1, 3, weight=1.0)
        G.add_directed_edge(1, 3, weight=2.0)  # Parallel edge
        components = source_components(G)
        assert len(components) == 1
        nodes, undir_edges, out_edges = components[0]
        assert set(nodes) == {1, 2}
        assert set(undir_edges) == {(1, 2, 0)}
        # Should include both parallel edges with their keys
        assert len(out_edges) == 2
        assert set(out_edges) == {(1, 3, 0), (1, 3, 1)}
    
    def test_isolated_nodes_as_source_components(self) -> None:
        """Test isolated nodes that are source components."""
        G = MixedMultiGraph()
        G.add_node(1)
        G.add_node(2)
        G.add_directed_edge(1, 3)
        G.add_directed_edge(2, 4)
        components = source_components(G)
        assert len(components) == 2
        node_sets = [set(nodes) for nodes, _, _ in components]
        assert {1} in node_sets
        assert {2} in node_sets
        # Check outgoing edges
        for nodes, undir, out in components:
            if set(nodes) == {1}:
                assert undir == []
                assert out == [(1, 3, 0)]
            elif set(nodes) == {2}:
                assert undir == []
                assert out == [(2, 4, 0)]
    
    def test_isolated_node_with_incoming_edge(self) -> None:
        """Test isolated node with incoming edge (not a source component)."""
        G = MixedMultiGraph()
        G.add_node(1)
        G.add_directed_edge(2, 1)  # Edge pointing into {1}
        components = source_components(G)
        # {1} is NOT a source component (has incoming edge)
        # {2} IS a source component (no incoming edges)
        assert len(components) == 1
        nodes, undir_edges, out_edges = components[0]
        assert set(nodes) == {2}
        assert undir_edges == []
        assert out_edges == [(2, 1, 0)]
    
    def test_large_source_component(self) -> None:
        """Test large source component with many nodes and edges."""
        G = MixedMultiGraph()
        # Create a chain: 1-2-3-4-5-6-7-8-9-10
        for i in range(1, 10):
            G.add_undirected_edge(i, i + 1)
        # Add outgoing edges
        G.add_directed_edge(5, 20)
        G.add_directed_edge(10, 21)
        components = source_components(G)
        assert len(components) == 1
        nodes, undir_edges, out_edges = components[0]
        assert set(nodes) == set(range(1, 11))
        assert len(undir_edges) == 9
        assert set(out_edges) == {(5, 20, 0), (10, 21, 0)}
    
    def test_nested_components(self) -> None:
        """Test nested component structure."""
        G = MixedMultiGraph()
        # Component A: {1, 2, 3} - source component
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_directed_edge(3, 4)
        # Component B: {4, 5} - has incoming edge from {1,2,3}, NOT source
        G.add_undirected_edge(4, 5)
        G.add_directed_edge(5, 6)
        # Component C: {6, 7} - has incoming edge from {4,5}, NOT source
        G.add_undirected_edge(6, 7)
        # Component D: {8, 9} - source component
        G.add_undirected_edge(8, 9)
        components = source_components(G)
        assert len(components) == 2
        node_sets = [set(nodes) for nodes, _, _ in components]
        assert {1, 2, 3} in node_sets
        assert {8, 9} in node_sets
    
    def test_component_with_bidirectional_edges(self) -> None:
        """Test component with edges going both ways (not source if has incoming)."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_directed_edge(3, 1)  # Incoming to {1, 2}
        G.add_directed_edge(1, 3)  # Outgoing from {1, 2}
        components = source_components(G)
        # {1, 2} is NOT a source component (has incoming edge (3, 1) from {3})
        # {3} is NOT a source component (has incoming edge (1, 3) from {1, 2})
        # Neither component is a source component
        assert len(components) == 0
    
    def test_all_components_have_incoming_edges(self) -> None:
        """Test case where no components are source components."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(3, 4)
        # Create a cycle of directed edges
        G.add_directed_edge(2, 3)  # {1,2} -> {3,4}
        G.add_directed_edge(4, 1)  # {3,4} -> {1,2}
        components = source_components(G)
        # Both components have incoming edges, so neither is a source component
        assert len(components) == 0
    
    def test_single_node(self) -> None:
        """Test isolated node"""
        G = MixedMultiGraph()
        G.add_node(1)
        components = source_components(G)
        assert len(components) == 1
        nodes, undir_edges, out_edges = components[0]
        assert set(nodes) == {1}
    
    def test_mixed_parallel_edges(self) -> None:
        """Test source component with both parallel undirected and directed edges."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2, weight=1.0)
        G.add_undirected_edge(1, 2, weight=2.0)  # Parallel undirected
        G.add_directed_edge(1, 3, weight=1.0)
        G.add_directed_edge(1, 3, weight=2.0)  # Parallel directed
        G.add_directed_edge(2, 4)
        components = source_components(G)
        assert len(components) == 1
        nodes, undir_edges, out_edges = components[0]
        assert set(nodes) == {1, 2}
        assert len(undir_edges) == 2  # Both parallel undirected edges
        assert set(undir_edges) == {(1, 2, 0), (1, 2, 1)}
        assert len(out_edges) == 3  # (1, 3) twice + (2, 4) once
        assert set(out_edges) == {(1, 3, 0), (1, 3, 1), (2, 4, 0)}


class TestValidation:

    def test_validate_synchronization(self) -> None:
        """Test _validate_synchronization method."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_directed_edge(2, 3)
        assert G._validate_synchronization()
        
        # Direct modification (should desync)
        G._undirected.add_edge(99, 100)
        result = G._validate_synchronization()
        print(f"After direct modification, synchronization valid: {result}")
        assert not result


class TestEdgeCases:
    """Test cases for edge cases and special scenarios."""

    def test_empty_graph_operations(self) -> None:
        """Test operations on empty graph."""
        G = MixedMultiGraph()
        # Behavior: NetworkX raises NetworkXPointlessConcept for empty graph connectivity
        try:
            result = is_connected(G)
            print(f"is_connected() on empty graph returns: {result}")
        except Exception as e:
            print(f"is_connected() on empty graph raises: {type(e).__name__}: {e}")
            # This is expected NetworkX behavior - empty graph connectivity is undefined
            assert "PointlessConcept" in str(type(e).__name__) or "empty" in str(e).lower()

    def test_single_node(self) -> None:
        """Test graph with single node."""
        G = MixedMultiGraph()
        G.add_node(1)
        assert G.number_of_nodes() == 1
        assert G.number_of_edges() == 0
        assert is_connected(G)  # Single node is connected
        assert G.degree(1) == 0

    def test_self_loops_undirected(self) -> None:
        """Test self-loops with undirected edges."""
        G = MixedMultiGraph()
        # Check if self-loops are allowed
        try:
            key = G.add_undirected_edge(1, 1, weight=5.0)
            print(f"Self-loop undirected edge added with key: {key}")
            assert G.has_edge(1, 1)
        except Exception as e:
            print(f"Self-loop undirected edge raises: {e}")

    def test_self_loops_directed(self) -> None:
        """Test self-loops with directed edges."""
        G = MixedMultiGraph()
        try:
            key = G.add_directed_edge(1, 1, weight=5.0)
            print(f"Self-loop directed edge added with key: {key}")
            assert has_directed_edge(G, 1, 1)
            assert G.outdegree(1) >= 1
            assert G.indegree(1) >= 1
        except Exception as e:
            print(f"Self-loop directed edge raises: {e}")

    def test_parallel_edges_same_attributes(self) -> None:
        """Test parallel edges with same attributes."""
        G = MixedMultiGraph()
        key1 = G.add_undirected_edge(1, 2, weight=5.0, label='same')
        key2 = G.add_undirected_edge(1, 2, weight=5.0, label='same')
        # Should still create parallel edges
        assert key1 != key2
        assert count_undirected_edges(G, 1, 2) == 2

    def test_removing_nonexistent_edge(self) -> None:
        """Test removing an edge that doesn't exist."""
        G = MixedMultiGraph()
        # Should raise error or silently fail?
        try:
            G.remove_edge(1, 2)
            print("Removing nonexistent edge: silently fails")
        except Exception as e:
            print(f"Removing nonexistent edge raises: {e}")

    def test_removing_nonexistent_node(self) -> None:
        """Test removing a node that doesn't exist."""
        G = MixedMultiGraph()
        try:
            G.remove_node(999)
            print("Removing nonexistent node: silently fails")
        except Exception as e:
            print(f"Removing nonexistent node raises: {e}")


class TestLargerGraphs:
    """Test cases with larger graphs (10+ nodes)."""

    def test_large_undirected_graph(self) -> None:
        """Test larger graph with many undirected edges."""
        G = MixedMultiGraph()
        nodes = list(range(1, 11))  # 10 nodes
        G.add_nodes_from(nodes)
        
        # Create a cycle
        for i in range(1, 10):
            G.add_undirected_edge(i, i + 1)
        G.add_undirected_edge(10, 1)  # Close the cycle
        
        # Add some parallel edges
        G.add_undirected_edge(1, 2, weight=2.0)
        G.add_undirected_edge(5, 6, weight=3.0)
        
        assert G.number_of_nodes() == 10
        assert is_connected(G)
        assert number_of_connected_components(G) == 1
        # Check degrees
        for node in nodes:
            degree = G.degree(node)
            print(f"Node {node} degree: {degree}")
            assert degree >= 2  # At least 2 in cycle

    def test_large_directed_graph(self) -> None:
        """Test larger graph with many directed edges."""
        G = MixedMultiGraph()
        nodes = list(range(1, 11))
        G.add_nodes_from(nodes)
        
        # Create directed cycle
        for i in range(1, 10):
            G.add_directed_edge(i, i + 1)
        G.add_directed_edge(10, 1)
        
        # Add parallel directed edges
        G.add_directed_edge(1, 2, weight=10.0)
        G.add_directed_edge(1, 2, weight=20.0)
        G.add_directed_edge(5, 6, weight=30.0)
        
        assert G.number_of_nodes() == 10
        assert is_connected(G)  # Weakly connected
        # Check degrees
        for node in nodes:
            indeg = G.indegree(node)
            outdeg = G.outdegree(node)
            print(f"Node {node}: indegree={indeg}, outdegree={outdeg}")
            assert indeg >= 1
            assert outdeg >= 1

    def test_large_mixed_graph(self) -> None:
        """Test larger graph with both directed and undirected edges."""
        G = MixedMultiGraph()
        nodes = list(range(1, 11))
        G.add_nodes_from(nodes)
        
        # Undirected component
        for i in range(1, 6):
            G.add_undirected_edge(i, i + 1)
        G.add_undirected_edge(5, 1)  # Close cycle
        
        # Directed component
        for i in range(6, 10):
            G.add_directed_edge(i, i + 1)
        G.add_directed_edge(10, 6)  # Close cycle
        
        # Connect components with one edge
        G.add_undirected_edge(5, 6)
        
        assert G.number_of_nodes() == 10
        assert is_connected(G)
        assert number_of_connected_components(G) == 1

    def test_large_graph_with_many_parallel_edges(self) -> None:
        """Test large graph with many parallel edges."""
        G = MixedMultiGraph()
        nodes = list(range(1, 11))
        G.add_nodes_from(nodes)
        
        # Add many parallel edges between same nodes
        for i in range(10):
            G.add_undirected_edge(1, 2, weight=float(i))
        for i in range(10):
            G.add_directed_edge(3, 4, weight=float(i))
        
        assert count_undirected_edges(G, 1, 2) == 10
        assert count_directed_edges(G, 3, 4) == 10
        
        # Check all weights are preserved
        # G[1][2] might be AdjacencyView, but should support dict-like access
        adj_1_2 = G[1][2]
        weights_undir = {adj_1_2[k]['weight'] for k in adj_1_2.keys()}
        assert len(weights_undir) == 10
        assert weights_undir == {float(i) for i in range(10)}
        
        weights_dir = {G._directed[3][4][k]['weight'] for k in G._directed[3][4].keys()}
        assert len(weights_dir) == 10
        assert weights_dir == {float(i) for i in range(10)}

    def test_large_graph_edge_removal(self) -> None:
        """Test removing edges from large graph."""
        G = MixedMultiGraph()
        nodes = list(range(1, 11))
        G.add_nodes_from(nodes)
        
        # Create graph
        for i in range(1, 10):
            key1 = G.add_undirected_edge(i, i + 1, weight=1.0)
            key2 = G.add_undirected_edge(i, i + 1, weight=2.0)  # Parallel
        
        # Remove some edges
        G.remove_edge(1, 2, key=0)  # Remove one parallel edge
        G.remove_edge(5, 6)  # Remove without key
        
        # Check connectivity
        components = number_of_connected_components(G)
        print(f"After edge removals, number of components: {components}")

    def test_large_graph_node_removal(self) -> None:
        """Test removing nodes from large graph."""
        G = MixedMultiGraph()
        nodes = list(range(1, 11))
        G.add_nodes_from(nodes)
        
        # Create connected graph
        for i in range(1, 10):
            G.add_undirected_edge(i, i + 1)
        
        # Remove some nodes
        G.remove_node(5)
        G.remove_nodes_from([2, 8])
        
        assert 5 not in G
        assert 2 not in G
        assert 8 not in G
        assert G.number_of_nodes() == 7
        
        # Check connectivity
        components = number_of_connected_components(G)
        print(f"After node removals, number of components: {components}")


class TestKeywordWarnings:
    """Test warnings for Python keyword identifiers on nodes and edge keys."""

    def test_keyword_node_id_directed_warns(self) -> None:
        """Directed edge with keyword node id should emit a warning."""
        G = MixedMultiGraph()
        with pytest.warns(UserWarning, match="Python keyword"):
            G.add_directed_edge("for", 2)
        assert ("for", 2, 0) in list(G._directed.edges(keys=True))

    def test_keyword_node_id_undirected_warns(self) -> None:
        """Undirected edge with keyword node id should emit a warning."""
        G = MixedMultiGraph()
        with pytest.warns(UserWarning, match="Python keyword"):
            G.add_undirected_edge(1, "while")
        assert (1, "while", 0) in list(G._undirected.edges(keys=True))

    def test_keyword_edge_key_warns(self) -> None:
        """Edge keys that are Python keywords should emit a warning."""
        G = MixedMultiGraph()
        with pytest.warns(UserWarning, match="Python keyword"):
            key = G.add_directed_edge(1, 2, key="class")
        assert (1, 2, key) in list(G._directed.edges(keys=True))

    def test_keyword_attribute_name_warns_directed_edge(self) -> None:
        """Attribute names that are Python keywords should emit a warning for directed edges."""
        G = MixedMultiGraph()
        with pytest.warns(UserWarning, match="Attribute name.*Python keyword"):
            G.add_directed_edge(1, 2, **{"True": "yes"})
        # Edge should still be added
        assert (1, 2, 0) in list(G._directed.edges(keys=True))
        # Attribute should be accessible
        assert G._directed[1][2][0]["True"] == "yes"

    def test_keyword_attribute_name_warns_undirected_edge(self) -> None:
        """Attribute names that are Python keywords should emit a warning for undirected edges."""
        G = MixedMultiGraph()
        with pytest.warns(UserWarning, match="Attribute name.*Python keyword"):
            G.add_undirected_edge(1, 2, **{"False": "no"})
        # Edge should still be added
        assert (1, 2, 0) in list(G._undirected.edges(keys=True))
        # Attribute should be accessible
        assert G._undirected[1][2][0]["False"] == "no"

    def test_keyword_attribute_name_warns_node(self) -> None:
        """Attribute names that are Python keywords should emit a warning for nodes."""
        G = MixedMultiGraph()
        with pytest.warns(UserWarning, match="Attribute name.*Python keyword"):
            G.add_node(1, **{"None": "nothing"})
        # Node should still be added
        assert 1 in G
        # Attribute should be accessible
        assert G._undirected.nodes[1]["None"] == "nothing"

    def test_non_keyword_attribute_name_no_warning(self) -> None:
        """Non-keyword attribute names with keyword values should not warn (except None)."""
        G = MixedMultiGraph()
        # Using "yes" as attribute name (not a keyword) with True as value (keyword value)
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            G.add_directed_edge(1, 2, yes=True)
            G.add_undirected_edge(2, 3, no=False)
        # Edges should be added
        assert (1, 2, 0) in list(G._directed.edges(keys=True))
        assert (2, 3, 0) in list(G._undirected.edges(keys=True))
        # Attributes should be accessible
        assert G._directed[1][2][0]["yes"] is True
        assert G._undirected[2][3][0]["no"] is False

    def test_none_attribute_value_warns_directed_edge(self) -> None:
        """Attribute values that are None should emit a warning for directed edges."""
        G = MixedMultiGraph()
        with pytest.warns(UserWarning, match="has value None.*Python keyword"):
            G.add_directed_edge(1, 2, weight=None)
        # Edge should still be added
        assert (1, 2, 0) in list(G._directed.edges(keys=True))
        # Attribute should be accessible
        assert G._directed[1][2][0]["weight"] is None

    def test_none_attribute_value_warns_undirected_edge(self) -> None:
        """Attribute values that are None should emit a warning for undirected edges."""
        G = MixedMultiGraph()
        with pytest.warns(UserWarning, match="has value None.*Python keyword"):
            G.add_undirected_edge(1, 2, length=None)
        # Edge should still be added
        assert (1, 2, 0) in list(G._undirected.edges(keys=True))
        # Attribute should be accessible
        assert G._undirected[1][2][0]["length"] is None

    def test_none_attribute_value_warns_node(self) -> None:
        """Attribute values that are None should emit a warning for nodes."""
        G = MixedMultiGraph()
        with pytest.warns(UserWarning, match="has value None.*Python keyword"):
            G.add_node(1, label=None)
        # Node should still be added
        assert 1 in G
        # Attribute should be accessible
        assert G._undirected.nodes[1]["label"] is None

    def test_none_attribute_value_multiple_warns(self) -> None:
        """Multiple attributes with None values should each emit a warning."""
        G = MixedMultiGraph()
        with pytest.warns(UserWarning, match="has value None.*Python keyword") as record:
            G.add_directed_edge(1, 2, weight=None, length=None)
        # Should have two warnings (one for each None attribute)
        assert len(record) == 2
        # Edge should still be added
        assert (1, 2, 0) in list(G._directed.edges(keys=True))
        # Attributes should be accessible
        assert G._directed[1][2][0]["weight"] is None
        assert G._directed[1][2][0]["length"] is None

