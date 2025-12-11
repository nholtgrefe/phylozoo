"""
Tests for suppress_degree2_node function in dm_operations.
"""

import pytest

from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
from phylozoo.core.primitives.d_multigraph.transformations import suppress_degree2_node


class TestSuppressDegree2NodeValid:
    """Tests for valid degree-2 node suppression."""

    def test_indegree1_outdegree1(self) -> None:
        """Test suppression of node with indegree=1 and outdegree=1."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        
        suppress_degree2_node(G, 2)
        
        assert 2 not in G.nodes()
        assert G.has_edge(1, 3)
        assert G.number_of_edges() == 1
        assert G.number_of_nodes() == 2


class TestSuppressDegree2NodeInvalidCombinations:
    """Tests for invalid edge configurations that raise errors."""

    def test_indegree2_outdegree0_raises_error(self) -> None:
        """Test that indegree=2, outdegree=0 raises ValueError."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(3, 2)
        
        with pytest.raises(ValueError, match="indegree 2 and outdegree 0"):
            suppress_degree2_node(G, 2)

    def test_indegree0_outdegree2_raises_error(self) -> None:
        """Test that indegree=0, outdegree=2 raises ValueError."""
        G = DirectedMultiGraph()
        G.add_edge(2, 1)
        G.add_edge(2, 3)
        
        with pytest.raises(ValueError, match="indegree 0 and outdegree 2"):
            suppress_degree2_node(G, 2)


class TestSuppressDegree2NodeErrorCases:
    """Tests for error cases (node not found, wrong degree, etc.)."""

    def test_node_not_in_graph_raises_error(self) -> None:
        """Test that suppressing a non-existent node raises ValueError."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        
        with pytest.raises(ValueError, match="not found in graph"):
            suppress_degree2_node(G, 99)

    def test_node_not_degree_2_raises_error(self) -> None:
        """Test that suppressing a non-degree-2 node raises ValueError."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        G.add_edge(2, 4)  # Node 2 now has degree 3
        
        with pytest.raises(ValueError, match="has degree 3, expected degree 2"):
            suppress_degree2_node(G, 2)

    def test_node_degree_1_raises_error(self) -> None:
        """Test that suppressing a degree-1 node raises ValueError."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        
        with pytest.raises(ValueError, match="has degree 1, expected degree 2"):
            suppress_degree2_node(G, 2)

    def test_node_indegree1_outdegree0_raises_error(self) -> None:
        """Test that indegree=1, outdegree=0 (but degree=1) raises error."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        
        with pytest.raises(ValueError, match="has degree 1, expected degree 2"):
            suppress_degree2_node(G, 2)

    def test_node_indegree0_outdegree1_raises_error(self) -> None:
        """Test that indegree=0, outdegree=1 (but degree=1) raises error."""
        G = DirectedMultiGraph()
        G.add_edge(2, 1)
        
        with pytest.raises(ValueError, match="has degree 1, expected degree 2"):
            suppress_degree2_node(G, 2)


class TestSuppressDegree2NodeAttributes:
    """Tests for attribute handling during suppression."""

    def test_attributes_preserved(self) -> None:
        """Test that attributes are preserved when suppressing edges."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0, label="edge1")
        G.add_edge(2, 3, weight=2.0, label="edge2")
        
        suppress_degree2_node(G, 2)
        
        # Check that edge exists and has attributes (outgoing edge overrides incoming)
        assert G.has_edge(1, 3)
        edge_data = dict(G._graph[1][3][0])
        assert edge_data.get('weight') == 2.0  # outgoing overrides
        assert edge_data.get('label') == "edge2"  # outgoing overrides

    def test_merged_attrs_provided(self) -> None:
        """Test that provided merged_attrs are used directly."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0)
        G.add_edge(2, 3, weight=2.0)
        
        merged_attrs = {'weight': 5.0, 'label': 'merged'}
        suppress_degree2_node(G, 2, merged_attrs=merged_attrs)
        
        # Check that merged_attrs are used
        assert G.has_edge(1, 3)
        edge_data = dict(G._graph[1][3][0])
        assert edge_data.get('weight') == 5.0
        assert edge_data.get('label') == 'merged'

    def test_attributes_incoming_first(self) -> None:
        """Test that incoming edge attributes come first, then outgoing overrides."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0, label="incoming", color="red")
        G.add_edge(2, 3, weight=2.0, label="outgoing")
        
        suppress_degree2_node(G, 2)
        
        # Incoming attributes first, then outgoing overrides
        assert G.has_edge(1, 3)
        edge_data = dict(G._graph[1][3][0])
        assert edge_data.get('weight') == 2.0  # outgoing overrides
        assert edge_data.get('label') == "outgoing"  # outgoing overrides
        assert edge_data.get('color') == "red"  # from incoming (not overridden)


class TestSuppressDegree2NodeParallelEdges:
    """Tests for parallel edge handling."""

    def test_parallel_edges_created(self) -> None:
        """Test that suppression can create parallel edges."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        G.add_edge(1, 3)  # Already exists
        
        suppress_degree2_node(G, 2)
        
        # Should have parallel edges between 1 and 3
        # Count edges between 1 and 3
        edge_count = len(list(G._graph.edges(1, 3, keys=True)))
        assert edge_count >= 1  # At least one edge exists
        assert G.has_edge(1, 3)  # Edge exists between 1 and 3

    def test_same_key_preserved(self) -> None:
        """Test that same keys are preserved when possible."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, key=5)
        G.add_edge(2, 3, key=5)
        
        suppress_degree2_node(G, 2)
        
        # Key should be preserved if same
        assert G.has_edge(1, 3, key=5)


class TestSuppressDegree2NodeComplex:
    """Tests for complex scenarios."""

    def test_suppress_in_larger_graph(self) -> None:
        """Test suppression in a larger graph structure."""
        G = DirectedMultiGraph()
        # Create a path: 1->2->3->4->5
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        G.add_edge(3, 4)
        G.add_edge(4, 5)
        
        # Suppress node 2
        suppress_degree2_node(G, 2)
        assert 2 not in G.nodes()
        assert G.has_edge(1, 3)
        
        # Suppress node 4 (now degree-2 after previous suppression)
        suppress_degree2_node(G, 4)
        assert 4 not in G.nodes()
        assert G.has_edge(3, 5)
        
        # Final graph should be 1->3->5
        assert G.number_of_nodes() == 3
        assert G.number_of_edges() == 2
        assert G.has_edge(1, 3)
        assert G.has_edge(3, 5)

    def test_suppress_with_branching(self) -> None:
        """Test suppression when node is part of a branching structure."""
        G = DirectedMultiGraph()
        # Structure: 1->2->3, 1->2->4
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        G.add_edge(2, 4)  # Node 2 has outdegree 2, so not suppressible
        
        # Node 2 is not degree-2 (it's degree 3)
        with pytest.raises(ValueError, match="has degree 3, expected degree 2"):
            suppress_degree2_node(G, 2)

