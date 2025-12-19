"""
Tests for identify_vertices function in d_multigraph transformations.
"""

import pytest

from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
from phylozoo.core.primitives.d_multigraph.transformations import identify_vertices


class TestIdentifyVerticesBasic:
    """Basic tests for identify_vertices."""

    def test_empty_vertices_list_raises_error(self) -> None:
        """Test that empty vertices list raises ValueError."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        
        with pytest.raises(ValueError, match="Vertices list cannot be empty"):
            identify_vertices(G, [])

    def test_single_vertex_no_change(self) -> None:
        """Test that identifying a single vertex does nothing."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        
        initial_nodes = list(G.nodes())
        identify_vertices(G, [1])
        
        assert list(G.nodes()) == initial_nodes
        assert G.has_edge(1, 2)

    def test_vertex_not_in_graph_raises_error(self) -> None:
        """Test that identifying a non-existent vertex raises ValueError."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        
        with pytest.raises(ValueError, match="not found in graph"):
            identify_vertices(G, [1, 99])

    def test_two_vertices_basic(self) -> None:
        """Test identifying two vertices."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        
        identify_vertices(G, [1, 2])
        
        assert 2 not in G.nodes()
        assert 1 in G.nodes()
        assert G.has_edge(1, 3)

    def test_multiple_vertices(self) -> None:
        """Test identifying multiple vertices."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        G.add_edge(3, 4)
        G.add_edge(4, 5)
        
        identify_vertices(G, [1, 2, 3])
        
        assert 1 in G.nodes()
        assert 2 not in G.nodes()
        assert 3 not in G.nodes()
        assert 4 in G.nodes()
        assert 5 in G.nodes()
        assert G.has_edge(1, 4)
        assert G.has_edge(4, 5)


class TestIdentifyVerticesNoSelfLoops:
    """Tests that identify_vertices does not create self-loops."""

    def test_no_self_loops_created(self) -> None:
        """Test that self-loops are not created."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 1)  # Bidirectional edge
        
        identify_vertices(G, [1, 2])
        
        # Should not have self-loop
        assert not G.has_edge(1, 1)
        assert 1 in G.nodes()
        assert 2 not in G.nodes()

    def test_self_loops_from_merged_vertices_removed(self) -> None:
        """Test that edges that would become self-loops are removed."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 1)
        G.add_edge(2, 3)
        
        identify_vertices(G, [1, 2])
        
        assert not G.has_edge(1, 1)
        assert G.has_edge(1, 3)


class TestIdentifyVerticesBidirectionalError:
    """Tests for bidirectional edge error cases."""

    def test_bidirectional_error_simple(self) -> None:
        """Test that creating bidirectional edges raises error."""
        G = DirectedMultiGraph()
        G.add_edge(1, 3)  # 1 -> 3
        G.add_edge(2, 1)  # 2 -> 1
        
        # After merging 1 and 2: would have 1 -> 3 and 1 -> 1 (self-loop, OK)
        # But wait, 2 -> 1 becomes nothing (self-loop)
        # Actually, this should be OK since 2->1 becomes self-loop which is removed
        
        # Let's try a case that would actually create bidirectional
        G2 = DirectedMultiGraph()
        G2.add_edge(1, 3)  # 1 -> 3
        G2.add_edge(2, 3)  # 2 -> 3
        G2.add_edge(3, 1)  # 3 -> 1
        
        # After merging 1 and 2: would have 1 -> 3 and 3 -> 1 (bidirectional!)
        with pytest.raises(ValueError, match="edges in both directions"):
            identify_vertices(G2, [1, 2])

    def test_bidirectional_error_complex(self) -> None:
        """Test bidirectional error with multiple vertices."""
        G = DirectedMultiGraph()
        G.add_edge(1, 4)  # 1 -> 4
        G.add_edge(2, 5)  # 2 -> 5
        G.add_edge(3, 4)  # 3 -> 4
        G.add_edge(5, 1)  # 5 -> 1
        
        # After merging [1, 2, 3]: would have 1 -> 4, 1 -> 5, and 5 -> 1 (bidirectional!)
        with pytest.raises(ValueError, match="edges in both directions"):
            identify_vertices(G, [1, 2, 3])


class TestIdentifyVerticesAttributes:
    """Tests for node attribute handling."""

    def test_preserve_first_vertex_attributes(self) -> None:
        """Test that first vertex's attributes are preserved by default."""
        G = DirectedMultiGraph()
        G.add_node(1, label='first', weight=1.0)
        G.add_node(2, label='second', weight=2.0)
        G.add_edge(1, 2)
        
        identify_vertices(G, [1, 2])
        
        assert 1 in G.nodes()
        node_attrs = G._graph.nodes[1]
        assert node_attrs['label'] == 'first'
        assert node_attrs['weight'] == 1.0

    def test_custom_merged_attributes(self) -> None:
        """Test using custom merged attributes."""
        G = DirectedMultiGraph()
        G.add_node(1, label='first', weight=1.0)
        G.add_node(2, label='second', weight=2.0)
        G.add_node(3, label='third', weight=3.0)
        G.add_edge(1, 4)
        G.add_edge(2, 4)
        G.add_edge(3, 4)
        
        merged_attrs = {'label': 'merged', 'weight': 5.0, 'custom': 'value'}
        identify_vertices(G, [1, 2, 3], merged_attrs=merged_attrs)
        
        assert 1 in G.nodes()
        node_attrs = G._graph.nodes[1]
        assert node_attrs['label'] == 'merged'
        assert node_attrs['weight'] == 5.0
        assert node_attrs['custom'] == 'value'
        # Original attributes should be removed
        assert 'first' not in str(node_attrs.values())


class TestIdentifyVerticesParallelEdges:
    """Tests that identify_vertices may create parallel edges."""

    def test_parallel_edges_created(self) -> None:
        """Test that parallel edges can be created."""
        G = DirectedMultiGraph()
        G.add_edge(1, 3, weight=1.0)
        G.add_edge(2, 3, weight=2.0)
        
        identify_vertices(G, [1, 2])
        
        # Should have parallel edges from 1 to 3
        assert G._graph.number_of_edges(1, 3) == 2
        assert G.has_edge(1, 3)

    def test_parallel_edges_preserve_attributes(self) -> None:
        """Test that parallel edges preserve their attributes."""
        G = DirectedMultiGraph()
        key1 = G.add_edge(1, 3, weight=1.0, label='first')
        key2 = G.add_edge(2, 3, weight=2.0, label='second')
        
        identify_vertices(G, [1, 2])
        
        # Should have two parallel edges with different attributes
        assert G._graph.number_of_edges(1, 3) == 2
        edges_data = [G._graph[1][3][k] for k in sorted(G._graph[1][3].keys())]
        weights = [d.get('weight') for d in edges_data]
        assert 1.0 in weights
        assert 2.0 in weights


class TestIdentifyVerticesComplex:
    """Complex tests for identify_vertices."""

    def test_preserves_other_edges(self) -> None:
        """Test that edges not involving merged vertices are preserved."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        G.add_edge(4, 5)
        G.add_edge(5, 6)
        
        identify_vertices(G, [1, 2])
        
        # Edges involving 4, 5, 6 should be unchanged
        assert G.has_edge(4, 5)
        assert G.has_edge(5, 6)
        assert G.has_edge(1, 3)  # New edge from merge

    def test_multiple_merges_separate(self) -> None:
        """Test multiple separate vertex identifications."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        G.add_edge(4, 5)
        G.add_edge(5, 6)
        
        identify_vertices(G, [1, 2])
        identify_vertices(G, [4, 5])
        
        assert 1 in G.nodes()
        assert 2 not in G.nodes()
        assert 4 in G.nodes()
        assert 5 not in G.nodes()
        assert G.has_edge(1, 3)
        assert G.has_edge(4, 6)

