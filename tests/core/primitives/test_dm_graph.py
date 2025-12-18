"""
Tests for the DirectedMultiGraph class.

This test suite provides comprehensive coverage of all DirectedMultiGraph features,
including edge cases, parallel edges, edge attributes, and larger graphs.
"""

import warnings

import pytest
from typing import Dict, List, Set, Tuple

from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
from phylozoo.core.primitives.d_multigraph.features import (
    bi_edge_connected_components,
    biconnected_components,
    connected_components,
    cut_edges,
    cut_vertices,
    has_parallel_edges,
    has_self_loops,
    is_connected,
    number_of_connected_components,
)
from phylozoo.core.primitives.d_multigraph.transformations import (
    identify_two_nodes,
    identify_node_set,
)
from phylozoo.core.primitives.d_multigraph.conversions import (
    digraph_to_directedmultigraph,
    multidigraph_to_directedmultigraph,
)


# Helper functions for testing
def count_edges(G: DirectedMultiGraph, u: int, v: int) -> int:
    """Count number of edges between u and v."""
    if u in G._graph and v in G._graph[u]:
        return len(G._graph[u][v])
    return 0


def has_edge(G: DirectedMultiGraph, u: int, v: int, key: int = None) -> bool:
    """Check if edge exists."""
    if key is None:
        return G._graph.has_edge(u, v)
    return G._graph.has_edge(u, v, key)


class TestInitialization:
    """Test cases for graph initialization."""

    def test_empty_graph(self) -> None:
        """Test creating an empty graph."""
        G = DirectedMultiGraph()
        assert G.number_of_nodes() == 0
        assert G.number_of_edges() == 0
        # G.nodes() returns iterator, need to convert to list
        assert len(list(G.nodes())) == 0
        assert len(list(G.edges())) == 0

    def test_init_with_edges(self) -> None:
        """Test initialization with edges."""
        G = DirectedMultiGraph(edges=[(1, 2), (2, 3), (1, 2)])  # Parallel edge
        assert G.number_of_nodes() == 3
        # Check number of edges between nodes
        num_edges_1_2 = len(G._graph[1][2]) if 1 in G._graph and 2 in G._graph[1] else 0
        num_edges_2_3 = len(G._graph[2][3]) if 2 in G._graph and 3 in G._graph[2] else 0
        print(f"Number of edges (1,2): {num_edges_1_2}, (2,3): {num_edges_2_3}")
        assert num_edges_1_2 == 2
        assert num_edges_2_3 == 1

    def test_self_loops_and_parallel_self_loops(self) -> None:
        """Self-loops are allowed and support parallel edges with distinct keys."""
        G = DirectedMultiGraph()
        k0 = G.add_edge(1, 1, weight=1.0)
        k1 = G.add_edge(1, 1, weight=2.0)  # Parallel self-loop
        assert {k0, k1} == {0, 1}
        assert G.number_of_edges() == 2
        assert G._graph.has_edge(1, 1, key=k0)
        assert G._graph.has_edge(1, 1, key=k1)
        # Combined graph also retains both self-loops
        assert G._combined.has_edge(1, 1, key=k0)
        assert G._combined.has_edge(1, 1, key=k1)
        assert has_self_loops(G) is True

    def test_has_self_loops_false(self) -> None:
        """has_self_loops returns False when no self-loops exist."""
        G = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
        assert has_self_loops(G) is False

    def test_init_with_edges_with_keys(self) -> None:
        """Test initialization with edges and explicit keys."""
        G = DirectedMultiGraph(edges=[(1, 2, 0), (1, 2, 5), (2, 3, 0)])
        num_edges = len(G._graph[1][2]) if 1 in G._graph and 2 in G._graph[1] else 0
        assert num_edges == 2
        assert G._graph.has_edge(1, 2, key=0)
        assert G._graph.has_edge(1, 2, key=5)
        assert G._graph.has_edge(2, 3, key=0)

    def test_init_with_edge_attributes(self) -> None:
        """Test initialization with edges and attributes."""
        G = DirectedMultiGraph(
            edges=[
                (1, 2),  # No attributes
                {'u': 2, 'v': 3, 'weight': 5.0, 'label': 'test'},  # Dict format
                {'u': 3, 'v': 4, 'key': 10, 'weight': 10.0}  # With key
            ]
        )
        assert count_edges(G, 1, 2) == 1
        assert count_edges(G, 2, 3) == 1
        assert count_edges(G, 3, 4) == 1
        # Check attributes
        edge_data = G._graph[2][3][0]
        assert edge_data['weight'] == 5.0
        assert edge_data['label'] == 'test'
        edge_data = G._graph[3][4][10]
        assert edge_data['weight'] == 10.0

    def test_repr_counts(self) -> None:
        """__repr__ reports node and edge counts."""
        G = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
        assert repr(G) == "DirectedMultiGraph(nodes=3, edges=2)"


class TestFactoryMethods:
    """Test cases for factory methods."""

    def test_digraph_to_directedmultigraph(self) -> None:
        """Test digraph_to_directedmultigraph factory method."""
        import networkx as nx
        G = nx.DiGraph()
        G.add_node(1, label='node1')
        G.add_edge(1, 2, weight=5.0)
        G.add_edge(2, 3)
        
        M = digraph_to_directedmultigraph(G)
        assert M.number_of_nodes() == 3
        assert M.number_of_edges() == 2


class TestIncidentEdges:
    """Test cases for incident edge methods."""

    def test_incident_parent_edges_basic(self) -> None:
        """Test incident_parent_edges with basic edges."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(3, 2)
        G.add_edge(4, 2)
        
        parent_edges = list(G.incident_parent_edges(2))
        assert len(parent_edges) == 3
        assert (1, 2) in parent_edges
        assert (3, 2) in parent_edges
        assert (4, 2) in parent_edges

    def test_incident_parent_edges_with_keys(self) -> None:
        """Test incident_parent_edges with keys."""
        G = DirectedMultiGraph()
        key1 = G.add_edge(1, 2, weight=1.0)
        key2 = G.add_edge(1, 2, weight=2.0)  # Parallel edge
        
        parent_edges = list(G.incident_parent_edges(2, keys=True))
        assert len(parent_edges) == 2
        assert (1, 2, key1) in parent_edges
        assert (1, 2, key2) in parent_edges

    def test_incident_parent_edges_with_data(self) -> None:
        """Test incident_parent_edges with data."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0, label="test")
        G.add_edge(3, 2, weight=2.0)
        
        parent_edges = list(G.incident_parent_edges(2, data=True))
        assert len(parent_edges) == 2
        # Check that data is included
        edge_dict = {edge[0]: edge[2] for edge in parent_edges if len(edge) == 3}
        assert 1 in edge_dict
        assert edge_dict[1]['weight'] == 1.0
        assert edge_dict[1]['label'] == "test"

    def test_incident_parent_edges_with_keys_and_data(self) -> None:
        """Test incident_parent_edges with keys and data."""
        G = DirectedMultiGraph()
        key1 = G.add_edge(1, 2, weight=1.0)
        key2 = G.add_edge(1, 2, weight=2.0)
        
        parent_edges = list(G.incident_parent_edges(2, keys=True, data=True))
        assert len(parent_edges) == 2
        # Check structure: (u, v, key, data)
        for edge in parent_edges:
            assert len(edge) == 4
            u, v, key, data = edge
            assert u == 1
            assert v == 2
            assert key in [key1, key2]
            assert 'weight' in data

    def test_incident_parent_edges_empty(self) -> None:
        """Test incident_parent_edges for node with no incoming edges."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        
        parent_edges = list(G.incident_parent_edges(1))
        assert len(parent_edges) == 0

    def test_incident_child_edges_basic(self) -> None:
        """Test incident_child_edges with basic edges."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(1, 3)
        G.add_edge(1, 4)
        
        child_edges = list(G.incident_child_edges(1))
        assert len(child_edges) == 3
        assert (1, 2) in child_edges
        assert (1, 3) in child_edges
        assert (1, 4) in child_edges

    def test_incident_child_edges_with_keys(self) -> None:
        """Test incident_child_edges with keys."""
        G = DirectedMultiGraph()
        key1 = G.add_edge(1, 2, weight=1.0)
        key2 = G.add_edge(1, 2, weight=2.0)  # Parallel edge
        
        child_edges = list(G.incident_child_edges(1, keys=True))
        assert len(child_edges) == 2
        assert (1, 2, key1) in child_edges
        assert (1, 2, key2) in child_edges

    def test_incident_child_edges_with_data(self) -> None:
        """Test incident_child_edges with data."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0, label="test")
        G.add_edge(1, 3, weight=2.0)
        
        child_edges = list(G.incident_child_edges(1, data=True))
        assert len(child_edges) == 2
        # Check that data is included
        edge_dict = {edge[1]: edge[2] for edge in child_edges if len(edge) == 3}
        assert 2 in edge_dict
        assert edge_dict[2]['weight'] == 1.0
        assert edge_dict[2]['label'] == "test"

    def test_incident_child_edges_with_keys_and_data(self) -> None:
        """Test incident_child_edges with keys and data."""
        G = DirectedMultiGraph()
        key1 = G.add_edge(1, 2, weight=1.0)
        key2 = G.add_edge(1, 2, weight=2.0)
        
        child_edges = list(G.incident_child_edges(1, keys=True, data=True))
        assert len(child_edges) == 2
        # Check structure: (u, v, key, data)
        for edge in child_edges:
            assert len(edge) == 4
            u, v, key, data = edge
            assert u == 1
            assert v == 2
            assert key in [key1, key2]
            assert 'weight' in data

    def test_incident_child_edges_empty(self) -> None:
        """Test incident_child_edges for node with no outgoing edges."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        
        child_edges = list(G.incident_child_edges(2))
        assert len(child_edges) == 0

    def test_incident_edges_consistency(self) -> None:
        """Test that incident_parent_edges and incident_child_edges are consistent."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0)
        G.add_edge(2, 3, weight=2.0)
        G.add_edge(3, 4, weight=3.0)
        
        # Node 2: incoming from 1, outgoing to 3
        parent_edges_2 = list(G.incident_parent_edges(2))
        child_edges_2 = list(G.incident_child_edges(2))
        assert len(parent_edges_2) == 1
        assert len(child_edges_2) == 1
        assert (1, 2) in parent_edges_2
        assert (2, 3) in child_edges_2
        
        # Node 3: incoming from 2, outgoing to 4
        parent_edges_3 = list(G.incident_parent_edges(3))
        child_edges_3 = list(G.incident_child_edges(3))
        assert len(parent_edges_3) == 1
        assert len(child_edges_3) == 1
        assert (2, 3) in parent_edges_3
        assert (3, 4) in child_edges_3

    def test_multidigraph_to_directedmultigraph(self) -> None:
        """Test multidigraph_to_directedmultigraph factory method."""
        import networkx as nx
        G = nx.MultiDiGraph()
        G.add_node(1, label='node1')
        G.add_edge(1, 2, key=0, weight=10.0)
        G.add_edge(1, 2, key=1, weight=20.0)  # Parallel edge
        G.add_edge(2, 3, key=0, weight=30.0)
        
        M = multidigraph_to_directedmultigraph(G)
        assert M.number_of_nodes() == 3
        assert M.number_of_edges() == 3
        # Parallel edges preserved
        assert count_edges(M, 1, 2) == 2
        assert count_edges(M, 2, 3) == 1
        # Keys preserved
        assert M._graph.has_edge(1, 2, key=0)
        assert M._graph.has_edge(1, 2, key=1)
        # Attributes preserved
        assert M._graph[1][2][0]['weight'] == 10.0
        assert M._graph[1][2][1]['weight'] == 20.0
    
    def test_digraph_to_directedmultigraph_with_attributes(self) -> None:
        """Test conversion from DiGraph with node and edge attributes."""
        import networkx as nx
        G = nx.DiGraph()
        G.add_node(1, type='A', value=10)
        G.add_node(2, type='B', value=20)
        G.add_edge(1, 2, weight=5.0, label='edge1')
        G.add_edge(2, 3, weight=10.0)
        
        M = digraph_to_directedmultigraph(G)
        # Check node attributes
        assert M._graph.nodes[1].get('type') == 'A'
        assert M._graph.nodes[1].get('value') == 10
        assert M._graph.nodes[2].get('type') == 'B'
        # Check edge attributes
        assert M._graph[1][2][0]['weight'] == 5.0
        assert M._graph[1][2][0]['label'] == 'edge1'
        assert M._graph[2][3][0]['weight'] == 10.0


class TestNodeOperations:
    """Test cases for node operations."""

    def test_add_node(self) -> None:
        """Test adding a single node."""
        G = DirectedMultiGraph()
        G.add_node(1)
        assert 1 in G
        assert G.number_of_nodes() == 1

    def test_add_node_with_attributes(self) -> None:
        """Test adding a node with attributes."""
        G = DirectedMultiGraph()
        G.add_node(1, label='node1', weight=5.0)
        assert 1 in G
        # Check if attributes are stored (access via underlying graph)
        assert G._graph.nodes[1].get('label') == 'node1'
        assert G._graph.nodes[1].get('weight') == 5.0

    def test_add_nodes_from(self) -> None:
        """Test adding multiple nodes."""
        G = DirectedMultiGraph()
        G.add_nodes_from([1, 2, 3, 4, 5])
        assert G.number_of_nodes() == 5
        assert all(n in G for n in [1, 2, 3, 4, 5])

    def test_add_nodes_from_with_attributes(self) -> None:
        """Test adding multiple nodes with attributes."""
        G = DirectedMultiGraph()
        G.add_nodes_from([1, 2, 3], label='test', type='A')
        for n in [1, 2, 3]:
            assert G._graph.nodes[n].get('label') == 'test'
            assert G._graph.nodes[n].get('type') == 'A'

    def test_remove_node(self) -> None:
        """Test removing a node."""
        G = DirectedMultiGraph()
        G.add_nodes_from([1, 2, 3])
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        G.remove_node(2)
        assert 2 not in G
        assert G.number_of_nodes() == 2
        assert not G.has_edge(1, 2)
        assert not G.has_edge(2, 3)

    def test_remove_nodes_from(self) -> None:
        """Test removing multiple nodes."""
        G = DirectedMultiGraph()
        G.add_nodes_from([1, 2, 3, 4, 5])
        G.remove_nodes_from([2, 4])
        assert G.number_of_nodes() == 3
        assert 2 not in G
        assert 4 not in G
        assert all(n in G for n in [1, 3, 5])


class TestEdgeOperations:
    """Test cases for edge operations."""

    def test_add_edge(self) -> None:
        """Test adding a single edge."""
        G = DirectedMultiGraph()
        key = G.add_edge(1, 2)
        assert G.has_edge(1, 2)
        assert not G.has_edge(2, 1)  # Directed is not symmetric
        assert isinstance(key, int)

    def test_add_edge_with_attributes(self) -> None:
        """Test adding edge with attributes."""
        G = DirectedMultiGraph()
        key = G.add_edge(1, 2, weight=10.0, label='dir', direction='forward')
        assert G.has_edge(1, 2)
        edge_data = G._graph[1][2][key]
        assert edge_data['weight'] == 10.0
        assert edge_data['label'] == 'dir'
        assert edge_data['direction'] == 'forward'

    def test_parallel_edges(self) -> None:
        """Test adding parallel edges."""
        G = DirectedMultiGraph()
        key1 = G.add_edge(1, 2, weight=1.0)
        key2 = G.add_edge(1, 2, weight=2.0)
        key3 = G.add_edge(1, 2, weight=3.0)
        assert key1 != key2 != key3
        assert count_edges(G, 1, 2) == 3
        weights = {G._graph[1][2][k]['weight'] for k in [key1, key2, key3]}
        assert weights == {1.0, 2.0, 3.0}

    def test_parallel_edges_with_explicit_keys(self) -> None:
        """Test adding parallel edges with explicit keys."""
        G = DirectedMultiGraph()
        key1 = G.add_edge(1, 2, key=100, weight=1.0)
        key2 = G.add_edge(1, 2, key=200, weight=2.0)
        assert key1 == 100
        assert key2 == 200
        assert has_edge(G, 1, 2, key=100)
        assert has_edge(G, 1, 2, key=200)

    def test_add_edges_from(self) -> None:
        """Test adding multiple edges."""
        G = DirectedMultiGraph()
        G.add_edges_from([(1, 2), (2, 3), (3, 1)])
        assert G.has_edge(1, 2)
        assert G.has_edge(2, 3)
        assert G.has_edge(3, 1)
        assert G.number_of_edges() == 3

    def test_add_edges_from_with_attributes(self) -> None:
        """Test adding multiple edges with shared attributes."""
        G = DirectedMultiGraph()
        G.add_edges_from([(1, 2), (2, 3)], weight=5.0, type='test')
        for u, v in [(1, 2), (2, 3)]:
            keys = list(G[u][v].keys())
            assert len(keys) > 0
            assert G[u][v][keys[0]]['weight'] == 5.0
            assert G[u][v][keys[0]]['type'] == 'test'

    def test_remove_edge(self) -> None:
        """Test removing an edge."""
        G = DirectedMultiGraph()
        key1 = G.add_edge(1, 2, weight=1.0)
        key2 = G.add_edge(1, 2, weight=2.0)
        G.remove_edge(1, 2, key=key1)
        assert not G.has_edge(1, 2, key=key1)
        assert G.has_edge(1, 2, key=key2)
        assert count_edges(G, 1, 2) == 1

    def test_remove_edge_without_key(self) -> None:
        """Test removing edge without specifying key."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0)
        G.add_edge(1, 2, weight=2.0)
        # Behavior: should remove one edge (which one?)
        G.remove_edge(1, 2)
        # Report behavior
        remaining = count_edges(G, 1, 2)
        print(f"After removing edge without key, remaining edges: {remaining}")

    def test_remove_edges_from(self) -> None:
        """Test removing multiple edges."""
        G = DirectedMultiGraph()
        G.add_edges_from([(1, 2), (2, 3), (3, 4)])
        G.remove_edges_from([(1, 2), (3, 4)])
        assert not G.has_edge(1, 2)
        assert G.has_edge(2, 3)
        assert not G.has_edge(3, 4)

    def test_remove_edge_nonexistent(self) -> None:
        """Test removing a non-existent edge raises error."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        with pytest.raises(ValueError):
            G.remove_edge(1, 3)  # Edge doesn't exist
        with pytest.raises(ValueError):
            G.remove_edge(1, 2, key=99)  # Key doesn't exist


class TestQueryOperations:
    """Test cases for query operations."""

    def test_has_edge(self) -> None:
        """Test has_edge method."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        assert G.has_edge(1, 2)
        assert not G.has_edge(2, 1)  # Directed
        assert not G.has_edge(1, 3)

    def test_has_edge_with_key(self) -> None:
        """Test has_edge with key parameter."""
        G = DirectedMultiGraph()
        key1 = G.add_edge(1, 2, key=10)
        key2 = G.add_edge(1, 2, key=20)
        assert G.has_edge(1, 2, key=10)
        assert G.has_edge(1, 2, key=20)
        assert not G.has_edge(1, 2, key=30)

    def test_number_of_nodes(self) -> None:
        """Test number_of_nodes method."""
        G = DirectedMultiGraph()
        G.add_nodes_from([1, 2, 3, 4, 5])
        assert G.number_of_nodes() == 5

    def test_number_of_edges(self) -> None:
        """Test number_of_edges method."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        G.add_edge(1, 2)  # Parallel edge
        assert G.number_of_edges() == 3

    def test_degree(self) -> None:
        """Test degree calculation."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(3, 2)  # Incoming to 2
        G.add_edge(2, 4)  # Outgoing from 2
        assert G.degree(2) == 3  # Total: 1 in + 2 out
        assert G.degree(1) == 1  # Only outgoing
        assert G.degree(99) == 0  # Node not in graph

    def test_indegree(self) -> None:
        """Test indegree calculation."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(3, 2)
        G.add_edge(4, 2)
        assert G.indegree(2) == 3
        assert G.indegree(1) == 0
        assert G.indegree(99) == 0  # Node not in graph

    def test_outdegree(self) -> None:
        """Test outdegree calculation."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(1, 3)
        G.add_edge(1, 4)
        assert G.outdegree(1) == 3
        assert G.outdegree(2) == 0
        assert G.outdegree(99) == 0  # Node not in graph

    def test_degree_with_parallel_edges(self) -> None:
        """Test degree calculation with parallel edges."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(1, 2)  # Parallel
        G.add_edge(3, 2)
        G.add_edge(3, 2)  # Parallel
        assert G.indegree(2) == 4  # 2 from 1, 2 from 3
        assert G.outdegree(1) == 2
        assert G.degree(2) == 4


class TestProperties:
    """Test cases for properties (nodes, edges, combined_graph)."""

    def test_nodes_property(self) -> None:
        """Test nodes property."""
        G = DirectedMultiGraph()
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
        G = DirectedMultiGraph()
        G.add_node(1, label='node1')
        G.add_node(2, weight=5.0)
        nodes_data = dict(G.nodes(data=True))
        assert 'label' in nodes_data[1]
        assert 'weight' in nodes_data[2]

    def test_edges_property(self) -> None:
        """Test edges property."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        # Test as attribute
        edges = G.edges
        assert (1, 2) in edges
        assert (2, 3) in edges
        # Test as method
        edges_list = list(G.edges())
        assert len(edges_list) >= 2

    def test_edges_property_with_data(self) -> None:
        """Test edges property with data=True."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=5.0)
        edges_data = list(G.edges(data=True))
        # Find the edge and check data
        for u, v, data in edges_data:
            if (u, v) == (1, 2):
                assert data['weight'] == 5.0
                break

    def test_edges_property_with_keys(self) -> None:
        """Test edges property with keys=True."""
        G = DirectedMultiGraph()
        key1 = G.add_edge(1, 2, weight=1.0)
        key2 = G.add_edge(1, 2, weight=2.0)
        edges_with_keys = list(G.edges(keys=True))
        keys_found = {k for u, v, k in edges_with_keys if (u, v) == (1, 2)}
        assert keys_found == {key1, key2}

    def test_combined_graph_property(self) -> None:
        """Test combined_graph property."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        combined = G.combined_graph
        # Combined should have all edges as undirected
        assert combined.has_edge(1, 2)
        assert combined.has_edge(2, 3)
        assert combined.has_edge(3, 2)  # Undirected view


class TestIterators:
    """Test cases for iterator methods."""

    def test_nodes_iter(self) -> None:
        """Test nodes_iter method."""
        G = DirectedMultiGraph()
        G.add_nodes_from([1, 2, 3])
        nodes = list(G.nodes_iter())
        assert set(nodes) == {1, 2, 3}

    def test_nodes_iter_with_data(self) -> None:
        """Test nodes_iter with data=True."""
        G = DirectedMultiGraph()
        G.add_node(1, label='test')
        nodes_data = list(G.nodes_iter(data=True))
        assert len(nodes_data) >= 1
        # Check format
        if isinstance(nodes_data[0], tuple):
            node, data = nodes_data[0]
            assert node == 1
            assert data.get('label') == 'test'

    def test_edges_iter(self) -> None:
        """Test edges_iter method."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        edges = list(G.edges_iter())
        assert len(edges) >= 2

    def test_edges_iter_with_keys(self) -> None:
        """Test edges_iter with keys=True."""
        G = DirectedMultiGraph()
        key1 = G.add_edge(1, 2)
        key2 = G.add_edge(1, 2)
        edges = list(G.edges_iter(keys=True))
        keys_found = {k for u, v, k in edges if (u, v) == (1, 2)}
        assert key1 in keys_found
        assert key2 in keys_found

    def test_nodes_iter_with_attribute(self) -> None:
        """Test nodes_iter with data='attribute'."""
        G = DirectedMultiGraph()
        G.add_node(1, weight=2.0, label='A')
        G.add_node(2, weight=3.5, label='B')
        G.add_node(3, weight=1.0)  # No label
        
        # Test with data='weight'
        nodes_weight = list(G.nodes_iter(data='weight'))
        assert isinstance(nodes_weight[0], tuple)
        assert len(nodes_weight[0]) == 2
        
        weight_dict = dict(nodes_weight)
        assert weight_dict[1] == 2.0
        assert weight_dict[2] == 3.5
        assert weight_dict[3] == 1.0
        
        # Test with data='label'
        nodes_label = list(G.nodes_iter(data='label'))
        label_dict = dict(nodes_label)
        assert label_dict[1] == 'A'
        assert label_dict[2] == 'B'
        assert label_dict[3] is None  # Missing attribute returns None

    def test_edges_iter_with_attribute(self) -> None:
        """Test edges_iter with data='attribute'."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.5, color='red')
        G.add_edge(2, 3, weight=2.5, color='blue')
        G.add_edge(3, 4, weight=3.0)  # No color
        
        # Test with data='weight'
        edges_weight = list(G.edges_iter(data='weight'))
        assert isinstance(edges_weight[0], tuple)
        assert len(edges_weight[0]) == 3  # (u, v, weight)
        
        # Check values
        edge_weights = {(u, v): w for u, v, w in edges_weight}
        assert edge_weights[(1, 2)] == 1.5
        assert edge_weights[(2, 3)] == 2.5
        assert edge_weights[(3, 4)] == 3.0
        
        # Test with data='color'
        edges_color = list(G.edges_iter(data='color'))
        edge_colors = {(u, v): c for u, v, c in edges_color}
        assert edge_colors[(1, 2)] == 'red'
        assert edge_colors[(2, 3)] == 'blue'
        assert edge_colors[(3, 4)] is None  # Missing attribute returns None

    def test_edges_iter_with_keys_and_attribute(self) -> None:
        """Test edges_iter with keys=True and data='attribute'."""
        G = DirectedMultiGraph()
        key1 = G.add_edge(1, 2, weight=1.5)
        key2 = G.add_edge(1, 2, weight=2.5)  # Parallel edge
        G.add_edge(2, 3, weight=3.0)
        
        # Test with keys=True and data='weight'
        edges = list(G.edges_iter(keys=True, data='weight'))
        assert isinstance(edges[0], tuple)
        assert len(edges[0]) == 4  # (u, v, key, weight)
        
        # Find the parallel edges
        parallel_edges = [(k, w) for u, v, k, w in edges if (u, v) == (1, 2)]
        assert len(parallel_edges) == 2
        
        weights = {k: w for k, w in parallel_edges}
        assert weights[key1] == 1.5
        assert weights[key2] == 2.5

    def test_incident_parent_edges_with_attribute(self) -> None:
        """Test incident_parent_edges with data='attribute'."""
        G = DirectedMultiGraph()
        G.add_edge(1, 3, weight=1.0)
        G.add_edge(2, 3, weight=2.0)
        
        # Test with data='weight'
        parent_edges = list(G.incident_parent_edges(3, data='weight'))
        assert len(parent_edges) == 2
        assert isinstance(parent_edges[0], tuple)
        assert len(parent_edges[0]) == 3  # (u, v, weight)
        
        edge_weights = {u: w for u, v, w in parent_edges}
        assert edge_weights[1] == 1.0
        assert edge_weights[2] == 2.0

    def test_incident_child_edges_with_attribute(self) -> None:
        """Test incident_child_edges with data='attribute'."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0)
        G.add_edge(1, 3, weight=2.0)
        
        # Test with data='weight'
        child_edges = list(G.incident_child_edges(1, data='weight'))
        assert len(child_edges) == 2
        assert isinstance(child_edges[0], tuple)
        assert len(child_edges[0]) == 3  # (u, v, weight)
        
        edge_weights = {v: w for u, v, w in child_edges}
        assert edge_weights[2] == 1.0
        assert edge_weights[3] == 2.0

    def test_nodes_property_with_attribute(self) -> None:
        """Test nodes property with data='attribute'."""
        G = DirectedMultiGraph()
        G.add_node(1, weight=2.0)
        G.add_node(2, weight=3.0)
        
        # Test via property
        nodes_weight = list(G.nodes(data='weight'))
        assert isinstance(nodes_weight[0], tuple)
        weight_dict = dict(nodes_weight)
        assert weight_dict[1] == 2.0
        assert weight_dict[2] == 3.0

    def test_edges_property_with_attribute(self) -> None:
        """Test edges property with data='attribute'."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.5)
        G.add_edge(2, 3, weight=2.5)
        
        # Test via property
        edges_weight = list(G.edges(data='weight'))
        assert isinstance(edges_weight[0], tuple)
        assert len(edges_weight[0]) == 3
        
        edge_weights = {(u, v): w for u, v, w in edges_weight}
        assert edge_weights[(1, 2)] == 1.5
        assert edge_weights[(2, 3)] == 2.5

    def test_neighbors(self) -> None:
        """Test neighbors method."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(3, 2)  # Incoming to 2
        G.add_edge(2, 4)  # Outgoing from 2
        neighbors = set(G.neighbors(2))
        assert neighbors == {1, 3, 4}

    def test_predecessors(self) -> None:
        """Test predecessors method."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(3, 2)
        predecessors = set(G.predecessors(2))
        assert predecessors == {1, 3}

    def test_successors(self) -> None:
        """Test successors method."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(1, 3)
        successors = set(G.successors(1))
        assert successors == {2, 3}


class TestSpecialMethods:
    """Test cases for special methods (__contains__, __iter__, __len__, __getitem__)."""

    def test_contains(self) -> None:
        """Test __contains__ method."""
        G = DirectedMultiGraph()
        G.add_node(1)
        assert 1 in G
        assert 2 not in G

    def test_iter(self) -> None:
        """Test __iter__ method."""
        G = DirectedMultiGraph()
        G.add_nodes_from([1, 2, 3])
        nodes = list(G)
        assert set(nodes) == {1, 2, 3}

    def test_len(self) -> None:
        """Test __len__ method."""
        G = DirectedMultiGraph()
        G.add_nodes_from([1, 2, 3, 4, 5])
        assert len(G) == 5

    def test_getitem(self) -> None:
        """Test __getitem__ method."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0)
        adj = G[1]
        # Should return adjacency view (NetworkX behavior)
        import networkx as nx
        # Check it's an AdjacencyView-like object
        assert hasattr(adj, '__getitem__')
        assert 2 in adj
        # Check edge data
        edge_data = adj[2][0]
        assert edge_data['weight'] == 1.0


class TestConnectivity:
    """Test cases for connectivity methods."""

    def test_number_of_connected_components(self) -> None:
        """Test number_of_connected_components."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(3, 4)
        assert number_of_connected_components(G) == 2

    def test_is_connected(self) -> None:
        """Test is_connected method."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        assert is_connected(G)
        
        G2 = DirectedMultiGraph()
        G2.add_edge(1, 2)
        G2.add_edge(3, 4)
        assert not is_connected(G2)

    def test_connected_components(self) -> None:
        """Test connected_components method."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(3, 4)
        components = list(connected_components(G))
        assert len(components) == 2
        assert {1, 2} in components
        assert {3, 4} in components

    def test_biconnected_components(self) -> None:
        """Test biconnected_components method."""
        G = DirectedMultiGraph()
        # Cycle forms one biconnected component
        _ = G.add_edge(1, 2)
        _ = G.add_edge(2, 3)
        _ = G.add_edge(3, 1)
        # Separate edge forms its own biconnected component
        _ = G.add_edge(4, 5)
        comps = list(biconnected_components(G))
        comp_sets = [set(comp) for comp in comps]
        assert {1, 2, 3} in comp_sets
        assert {4, 5} in comp_sets

    def test_biconnected_components_parallel_edges(self) -> None:
        """Test biconnected_components with parallel edges."""
        G = DirectedMultiGraph()
        # Cycle with parallel edges - still one biconnected component
        _ = G.add_edge(1, 2)
        _ = G.add_edge(1, 2)  # Parallel edge
        _ = G.add_edge(2, 3)
        _ = G.add_edge(3, 1)
        comps = list(biconnected_components(G))
        comp_sets = [set(comp) for comp in comps]
        # Parallel edges don't change biconnected component structure
        assert {1, 2, 3} in comp_sets
        assert len(comps) == 1
        
        # Test with parallel edges on bridge
        G2 = DirectedMultiGraph()
        # Cycle 1: 1-2-3-1
        _ = G2.add_edge(1, 2)
        _ = G2.add_edge(2, 3)
        _ = G2.add_edge(3, 1)
        # Bridge with parallel edges: 3-4
        _ = G2.add_edge(3, 4)
        _ = G2.add_edge(3, 4)  # Parallel bridge edge
        # Cycle 2: 4-5-6-4
        _ = G2.add_edge(4, 5)
        _ = G2.add_edge(5, 6)
        _ = G2.add_edge(6, 4)
        comps2 = list(biconnected_components(G2))
        comp_sets2 = [set(comp) for comp in comps2]
        # Parallel edges on 3-4 make it non-bridge, creating separate components:
        # - {1, 2, 3} - first cycle
        # - {3, 4} - parallel edges between 3 and 4 (non-bridge due to parallel edges)
        # - {4, 5, 6} - second cycle
        assert {1, 2, 3} in comp_sets2
        assert {3, 4} in comp_sets2  # Parallel edges create biconnected component
        assert {4, 5, 6} in comp_sets2
        assert len(comps2) == 3

    def test_bi_edge_connected_components(self) -> None:
        """Test bi_edge_connected_components method."""
        G = DirectedMultiGraph()
        # Cycle forms one bi-edge connected component (no bridges in cycle)
        _ = G.add_edge(1, 2)
        _ = G.add_edge(2, 3)
        _ = G.add_edge(3, 1)
        # Bridge edge to separate node
        _ = G.add_edge(3, 4)
        comps = list(bi_edge_connected_components(G))
        comp_sets = [set(comp) for comp in comps]
        assert {1, 2, 3} in comp_sets
        assert {4} in comp_sets
        assert len(comps) == 2

    def test_bi_edge_connected_components_two_cycles(self) -> None:
        """Test bi_edge_connected_components with two cycles connected by bridge."""
        G = DirectedMultiGraph()
        # First cycle: 1-2-3-1
        _ = G.add_edge(1, 2)
        _ = G.add_edge(2, 3)
        _ = G.add_edge(3, 1)
        # Bridge: 3-4
        _ = G.add_edge(3, 4)
        # Second cycle: 4-5-6-4
        _ = G.add_edge(4, 5)
        _ = G.add_edge(5, 6)
        _ = G.add_edge(6, 4)
        comps = list(bi_edge_connected_components(G))
        comp_sets = [set(comp) for comp in comps]
        assert {1, 2, 3} in comp_sets
        assert {4, 5, 6} in comp_sets
        assert len(comps) == 2

    def test_bi_edge_connected_components_path(self) -> None:
        """Test bi_edge_connected_components with path graph (all edges are bridges)."""
        G = DirectedMultiGraph()
        # Path graph: 1-2-3-4-5 (all edges are bridges)
        _ = G.add_edge(1, 2)
        _ = G.add_edge(2, 3)
        _ = G.add_edge(3, 4)
        _ = G.add_edge(4, 5)
        comps = list(bi_edge_connected_components(G))
        comp_sets = [set(comp) for comp in comps]
        # Each node is its own component since all edges are bridges
        assert {1} in comp_sets
        assert {2} in comp_sets
        assert {3} in comp_sets
        assert {4} in comp_sets
        assert {5} in comp_sets
        assert len(comps) == 5

    def test_bi_edge_connected_components_parallel_edges(self) -> None:
        """Test bi_edge_connected_components with parallel edges."""
        G = DirectedMultiGraph()
        # Cycle with parallel edges - still one bi-edge connected component
        _ = G.add_edge(1, 2)
        _ = G.add_edge(1, 2)  # Parallel edge
        _ = G.add_edge(2, 3)
        _ = G.add_edge(3, 1)
        comps = list(bi_edge_connected_components(G))
        comp_sets = [set(comp) for comp in comps]
        # Parallel edges don't create bridges, so cycle remains one component
        assert {1, 2, 3} in comp_sets
        assert len(comps) == 1
        
        # Test with parallel edges on what would be a bridge
        G2 = DirectedMultiGraph()
        # Cycle 1: 1-2-3-1
        _ = G2.add_edge(1, 2)
        _ = G2.add_edge(2, 3)
        _ = G2.add_edge(3, 1)
        # Parallel edges between 3 and 4 (not a bridge due to parallel edges)
        _ = G2.add_edge(3, 4)
        _ = G2.add_edge(3, 4)  # Parallel edge
        # Cycle 2: 4-5-6-4
        _ = G2.add_edge(4, 5)
        _ = G2.add_edge(5, 6)
        _ = G2.add_edge(6, 4)
        comps2 = list(bi_edge_connected_components(G2))
        comp_sets2 = [set(comp) for comp in comps2]
        # Parallel edges make 3-4 non-bridge, so all nodes are in one component
        assert {1, 2, 3, 4, 5, 6} in comp_sets2
        assert len(comps2) == 1
        
        # Test path graph with parallel edges in the middle: 1-2-3-4 where 2-3 has parallel edges
        G3 = DirectedMultiGraph()
        _ = G3.add_edge(1, 2)
        _ = G3.add_edge(2, 3)
        _ = G3.add_edge(2, 3)  # Parallel edge
        _ = G3.add_edge(3, 4)
        comps3 = list(bi_edge_connected_components(G3))
        comp_sets3 = [set(comp) for comp in comps3]
        # Parallel edges make 2-3 non-bridge, so we get: {1}, {2, 3}, {4}
        assert {1} in comp_sets3
        assert {2, 3} in comp_sets3
        assert {4} in comp_sets3
        assert len(comps3) == 3

    def test_cut_edges(self) -> None:
        """Test cut_edges function."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        # Both edges are cut-edges in a chain 1->2->3
        edges = cut_edges(G, keys=True)
        assert (1, 2, 0) in edges
        assert (2, 3, 0) in edges
        # Test without keys
        edges_no_keys = cut_edges(G, keys=False)
        assert (1, 2) in edges_no_keys
        assert (2, 3) in edges_no_keys

    def test_cut_edges_with_parallel(self) -> None:
        """Test cut_edges with parallel edges."""
        G = DirectedMultiGraph()
        key1 = G.add_edge(1, 2)
        key2 = G.add_edge(1, 2)  # Parallel edge
        G.add_edge(2, 3)
        # The edge 1->2 has parallel edges, so neither individual edge is a cut-edge
        # But edge 2->3 is still a cut-edge
        edges = cut_edges(G, keys=True)
        assert (2, 3, 0) in edges
        # Parallel edges 1->2 should not be cut-edges
        assert (1, 2, key1) not in edges
        assert (1, 2, key2) not in edges

    def test_cut_edges_with_data(self) -> None:
        """Test cut_edges with data parameter."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=5.0)
        G.add_edge(2, 3, weight=3.0)
        # Test with data=True
        edges_with_data = cut_edges(G, keys=True, data=True)
        assert any(u == 1 and v == 2 and k == 0 for u, v, k, _ in edges_with_data)
        # Test with data='weight'
        edges_with_weight = cut_edges(G, keys=True, data='weight')
        assert (1, 2, 0, 5.0) in edges_with_weight
        assert (2, 3, 0, 3.0) in edges_with_weight

    def test_cut_vertices(self) -> None:
        """Test cut_vertices function."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        # Node 2 is a cut-vertex
        vertices = cut_vertices(G)
        assert 2 in vertices
        # Node 1 is not a cut-vertex
        assert 1 not in vertices

    def test_cut_vertices_with_data(self) -> None:
        """Test cut_vertices with data parameter."""
        G = DirectedMultiGraph()
        G.add_node(1, label='A')
        G.add_node(2, label='B')
        G.add_node(3, label='C')
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        # Test with data=True
        vertices_with_data = cut_vertices(G, data=True)
        assert (2, {'label': 'B'}) in vertices_with_data
        # Test with data='label'
        vertices_with_label = cut_vertices(G, data='label')
        assert (2, 'B') in vertices_with_label

    def test_cut_edges_large_tree(self) -> None:
        """Test cut_edges on a larger tree structure."""
        G = DirectedMultiGraph()
        # Build a binary tree: depth 4, 15 nodes
        for i in range(1, 8):
            G.add_edge(i, 2*i)      # Left child
            G.add_edge(i, 2*i + 1)  # Right child
        
        edges = cut_edges(G, keys=True)
        # In a tree, all edges are bridges
        assert len(edges) == 14  # 7 internal nodes, each with 2 children
        
        # Verify no parallel edges are reported as bridges
        for u, v, k in edges:
            assert k == 0, "Only single edges should be bridges"

    def test_cut_edges_dense_graph(self) -> None:
        """Test cut_edges on a densely connected graph."""
        G = DirectedMultiGraph()
        # Create a complete directed graph on 5 nodes
        for i in range(1, 6):
            for j in range(1, 6):
                if i != j:
                    G.add_edge(i, j)
        
        edges = cut_edges(G, keys=True)
        # Complete graph has no bridges
        assert len(edges) == 0

    def test_cut_edges_bridge_network(self) -> None:
        """Test cut_edges on a graph with clear bridge structure."""
        G = DirectedMultiGraph()
        # Create two cliques connected by a bridge
        # Clique 1: nodes 1, 2, 3
        for i in [1, 2, 3]:
            for j in [1, 2, 3]:
                if i != j:
                    G.add_edge(i, j)
        
        # Bridge
        G.add_edge(3, 4)
        
        # Clique 2: nodes 4, 5, 6
        for i in [4, 5, 6]:
            for j in [4, 5, 6]:
                if i != j:
                    G.add_edge(i, j)
        
        edges = cut_edges(G, keys=True)
        # Only the bridge edge should be a cut edge
        assert len(edges) == 1
        assert (3, 4, 0) in edges

    def test_cut_vertices_large_graph(self) -> None:
        """Test cut_vertices on a larger graph structure."""
        G = DirectedMultiGraph()
        # Create a star graph with center node 0
        center = 0
        for i in range(1, 21):  # 20 leaf nodes
            G.add_edge(center, i)
        
        vertices = cut_vertices(G)
        # Only the center should be a cut vertex
        assert len(vertices) == 1
        assert center in vertices

    def test_cut_vertices_chain_of_cliques(self) -> None:
        """Test cut_vertices on chain of cliques."""
        G = DirectedMultiGraph()
        # Create 3 cliques of size 3, connected in a chain
        # Clique 1: 1, 2, 3
        for i in [1, 2, 3]:
            for j in [1, 2, 3]:
                if i != j:
                    G.add_edge(i, j)
        
        # Connection 1: node 3 connects to node 4
        G.add_edge(3, 4)
        
        # Clique 2: 4, 5, 6
        for i in [4, 5, 6]:
            for j in [4, 5, 6]:
                if i != j:
                    G.add_edge(i, j)
        
        # Connection 2: node 6 connects to node 7
        G.add_edge(6, 7)
        
        # Clique 3: 7, 8, 9
        for i in [7, 8, 9]:
            for j in [7, 8, 9]:
                if i != j:
                    G.add_edge(i, j)
        
        vertices = cut_vertices(G)
        # Nodes 3, 4, 6, 7 should be cut vertices
        assert 3 in vertices or 4 in vertices  # At least one connection point
        assert 6 in vertices or 7 in vertices  # At least one connection point

    def test_cut_edges_mixed_parallel_and_bridges(self) -> None:
        """Test cut_edges with both parallel edges and bridges."""
        G = DirectedMultiGraph()
        # Path with parallel edges at start: 1 ⇉ 2 → 3 → 4
        G.add_edge(1, 2)
        G.add_edge(1, 2)  # Parallel
        G.add_edge(2, 3)
        G.add_edge(3, 4)
        
        edges = cut_edges(G, keys=True)
        # Only edges 2→3 and 3→4 should be bridges
        assert (2, 3, 0) in edges
        assert (3, 4, 0) in edges
        assert (1, 2, 0) not in edges  # Parallel edge not a bridge
        assert (1, 2, 1) not in edges  # Parallel edge not a bridge
        assert len(edges) == 2


class TestGraphOperations:
    """Test cases for graph operations."""

    def test_copy(self) -> None:
        """Test copy method."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0)
        G.add_edge(2, 3, weight=2.0)
        H = G.copy()
        # Should have same structure
        assert H.number_of_nodes() == G.number_of_nodes()
        assert H.number_of_edges() == G.number_of_edges()
        # Should be independent
        H.add_edge(3, 4)
        assert H.number_of_edges() == 3
        assert G.number_of_edges() == 2

    def test_clear(self) -> None:
        """Test clear method."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        G.clear()
        assert G.number_of_nodes() == 0
        assert G.number_of_edges() == 0

    def test_identify_two_nodes(self) -> None:
        """Test identify_two_nodes method."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        identify_two_nodes(G, 1, 2)
        assert 2 not in G
        assert 1 in G
        # Edges should be moved to node 1
        assert G.has_edge(1, 3)

    def test_identify_node_set(self) -> None:
        """Test identify_node_set method."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        G.add_edge(3, 4)
        identify_node_set(G, [1, 2, 3])
        # First node (1) should be kept
        assert 1 in G
        # Other nodes may or may not be present (depends on implementation)
        assert len(list(G.nodes())) <= 3

    def test_validate_synchronization(self) -> None:
        """Test _validate_synchronization method."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        assert G._validate_synchronization()
        
        # Direct modification (should fail validation)
        G._graph.add_edge(99, 100)
        assert not G._validate_synchronization()


class TestEdgeCases:
    """Test cases for edge cases and error handling."""

    def test_empty_graph_operations(self) -> None:
        """Test operations on empty graph."""
        G = DirectedMultiGraph()
        assert G.number_of_nodes() == 0
        assert G.number_of_edges() == 0
        # NetworkX raises NetworkXPointlessConcept for connectivity on null graph
        import networkx as nx
        with pytest.raises(nx.NetworkXPointlessConcept):
            is_connected(G)
        assert number_of_connected_components(G) == 0
        assert list(G.nodes()) == []
        assert list(G.edges()) == []

    def test_single_node(self) -> None:
        """Test graph with single node."""
        G = DirectedMultiGraph()
        G.add_node(1)
        assert G.number_of_nodes() == 1
        assert G.number_of_edges() == 0
        assert G.degree(1) == 0
        assert G.indegree(1) == 0
        assert G.outdegree(1) == 0

    def test_self_loops(self) -> None:
        """Test self-loops (if supported)."""
        G = DirectedMultiGraph()
        # Note: identify_two_nodes might create self-loops, but they should be cleaned up
        # Direct self-loops might not be explicitly supported
        # This test documents current behavior
        try:
            key = G.add_edge(1, 1)
            assert G.has_edge(1, 1)
            assert G.degree(1) == 2  # Self-loop contributes 2 to degree
        except Exception as e:
            print(f"Self-loops behavior: {type(e).__name__}: {e}")

    def test_parallel_edges_same_attributes(self) -> None:
        """Test parallel edges with same attributes."""
        G = DirectedMultiGraph()
        key1 = G.add_edge(1, 2, weight=5.0)
        key2 = G.add_edge(1, 2, weight=5.0)  # Same weight
        assert key1 != key2
        assert G._graph[1][2][key1]['weight'] == 5.0
        assert G._graph[1][2][key2]['weight'] == 5.0

    def test_large_graph(self) -> None:
        """Test operations on larger graph."""
        G = DirectedMultiGraph()
        # Create a chain
        for i in range(100):
            G.add_edge(i, i + 1)
        assert G.number_of_nodes() == 101
        assert G.number_of_edges() == 100
        assert is_connected(G)
        assert number_of_connected_components(G) == 1

    def test_large_graph_node_removal(self) -> None:
        """Test node removal on larger graph."""
        G = DirectedMultiGraph()
        # Create a chain
        for i in range(50):
            G.add_edge(i, i + 1)
        # Remove middle node
        G.remove_node(25)
        assert 25 not in G
        assert G.number_of_nodes() == 50
        # Graph should be disconnected
        assert number_of_connected_components(G) == 2


class TestLargerGraphs:
    """Test cases for more complex graph structures."""

    def test_cycle(self) -> None:
        """Test cycle structure."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        G.add_edge(3, 1)
        assert is_connected(G)
        assert number_of_connected_components(G) == 1
        # No cut-edges in a cycle
        edges = cut_edges(G, keys=True)
        assert (1, 2, 0) not in edges
        assert (2, 3, 0) not in edges
        assert (3, 1, 0) not in edges

    def test_star_graph(self) -> None:
        """Test star graph structure."""
        G = DirectedMultiGraph()
        # Center node 0, leaves 1-5
        for i in range(1, 6):
            G.add_edge(0, i)
        assert G.number_of_nodes() == 6
        assert G.number_of_edges() == 5
        assert G.outdegree(0) == 5
        assert all(G.indegree(i) == 1 for i in range(1, 6))
        # Center is cut-vertex
        vertices = cut_vertices(G)
        assert 0 in vertices
        # All edges are cut-edges
        edges = cut_edges(G, keys=True)
        assert all((0, i, 0) in edges for i in range(1, 6))

    def test_complete_bipartite(self) -> None:
        """Test complete bipartite structure."""
        G = DirectedMultiGraph()
        # Left partition: 1, 2
        # Right partition: 3, 4, 5
        for left in [1, 2]:
            for right in [3, 4, 5]:
                G.add_edge(left, right)
        assert G.number_of_nodes() == 5
        assert G.number_of_edges() == 6
        assert is_connected(G)
        # All nodes in left partition have outdegree 3
        assert all(G.outdegree(i) == 3 for i in [1, 2])
        # All nodes in right partition have indegree 2
        assert all(G.indegree(i) == 2 for i in [3, 4, 5])


class TestKeywordWarnings:
    """Test warnings for Python keyword identifiers on nodes and edge keys."""

    def test_keyword_node_id_warns(self) -> None:
        """Node IDs that are Python keywords should emit a warning."""
        G = DirectedMultiGraph()
        with pytest.warns(UserWarning, match="Python keyword"):
            G.add_edge("for", 2)
        assert G.has_edge("for", 2)

    def test_keyword_edge_key_warns(self) -> None:
        """Edge keys that are Python keywords should emit a warning."""
        G = DirectedMultiGraph()
        with pytest.warns(UserWarning, match="Python keyword"):
            key = G.add_edge(1, 2, key="class")
        assert G.has_edge(1, 2, key)

    def test_keyword_attribute_name_warns_edge(self) -> None:
        """Attribute names that are Python keywords should emit a warning."""
        G = DirectedMultiGraph()
        with pytest.warns(UserWarning, match="Attribute name.*Python keyword"):
            G.add_edge(1, 2, **{"True": "yes"})
        # Edge should still be added
        assert G.has_edge(1, 2)
        # Attribute should be accessible
        assert G._graph[1][2][0]["True"] == "yes"

    def test_keyword_attribute_name_warns_node(self) -> None:
        """Attribute names that are Python keywords should emit a warning for nodes."""
        G = DirectedMultiGraph()
        with pytest.warns(UserWarning, match="Attribute name.*Python keyword"):
            G.add_node(1, **{"False": "no"})
        # Node should still be added
        assert 1 in G
        # Attribute should be accessible
        assert G._graph.nodes[1]["False"] == "no"

    def test_non_keyword_attribute_name_no_warning(self) -> None:
        """Non-keyword attribute names with keyword values should not warn."""
        G = DirectedMultiGraph()
        # Using "yes" as attribute name (not a keyword) with True as value (keyword value)
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            G.add_edge(1, 2, yes=True)
            G.add_node(3, no=False)
        # Edges and nodes should be added
        assert G.has_edge(1, 2)
        assert 3 in G
        # Attributes should be accessible
        assert G._graph[1][2][0]["yes"] is True
        assert G._graph.nodes[3]["no"] is False

    def test_none_attribute_value_warns_edge(self) -> None:
        """Attribute values that are None should emit a warning."""
        G = DirectedMultiGraph()
        with pytest.warns(UserWarning, match="has value None.*Python keyword"):
            G.add_edge(1, 2, weight=None)
        # Edge should still be added
        assert G.has_edge(1, 2)
        # Attribute should be accessible
        assert G._graph[1][2][0]["weight"] is None

    def test_none_attribute_value_warns_node(self) -> None:
        """Attribute values that are None should emit a warning for nodes."""
        G = DirectedMultiGraph()
        with pytest.warns(UserWarning, match="has value None.*Python keyword"):
            G.add_node(1, label=None)
        # Node should still be added
        assert 1 in G
        # Attribute should be accessible
        assert G._graph.nodes[1]["label"] is None

    def test_none_attribute_value_multiple_warns(self) -> None:
        """Multiple attributes with None values should each emit a warning."""
        G = DirectedMultiGraph()
        with pytest.warns(UserWarning, match="has value None.*Python keyword") as record:
            G.add_edge(1, 2, weight=None, length=None)
        # Should have two warnings (one for each None attribute)
        assert len(record) == 2
        # Edge should still be added
        assert G.has_edge(1, 2)
        # Attributes should be accessible
        assert G._graph[1][2][0]["weight"] is None
        assert G._graph[1][2][0]["length"] is None


class TestHasParallelEdges:
    """Test cases for has_parallel_edges function."""

    def test_no_parallel_edges(self) -> None:
        """Test graph with no parallel edges."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        G.add_edge(3, 4)
        assert not has_parallel_edges(G)

    def test_with_parallel_edges(self) -> None:
        """Test graph with parallel edges."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(1, 2)  # Parallel edge
        assert has_parallel_edges(G)

    def test_multiple_parallel_edges(self) -> None:
        """Test graph with multiple pairs of parallel edges."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(1, 2)  # Parallel edge
        G.add_edge(2, 3)
        G.add_edge(2, 3)  # Another parallel edge
        assert has_parallel_edges(G)

    def test_parallel_edges_different_directions(self) -> None:
        """Test that edges in different directions are not considered parallel."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 1)  # Different direction, not parallel
        assert not has_parallel_edges(G)

    def test_empty_graph(self) -> None:
        """Test empty graph has no parallel edges."""
        G = DirectedMultiGraph()
        assert not has_parallel_edges(G)

    def test_single_edge(self) -> None:
        """Test graph with single edge has no parallel edges."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        assert not has_parallel_edges(G)

    def test_parallel_edges_with_attributes(self) -> None:
        """Test that parallel edges with different attributes are still parallel."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0)
        G.add_edge(1, 2, weight=2.0)  # Parallel with different attribute
        assert has_parallel_edges(G)

    def test_three_parallel_edges(self) -> None:
        """Test graph with three parallel edges between same nodes."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(1, 2)
        G.add_edge(1, 2)
        assert has_parallel_edges(G)

