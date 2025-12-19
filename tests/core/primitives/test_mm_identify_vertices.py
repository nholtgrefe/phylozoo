"""
Tests for identify_vertices function in m_multigraph transformations.
"""

import pytest

from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
from phylozoo.core.primitives.m_multigraph.transformations import identify_vertices


class TestIdentifyVerticesBasic:
    """Basic tests for identify_vertices."""

    def test_empty_vertices_list_raises_error(self) -> None:
        """Test that empty vertices list raises ValueError."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        
        with pytest.raises(ValueError, match="Vertices list cannot be empty"):
            identify_vertices(G, [])

    def test_single_vertex_no_change(self) -> None:
        """Test that identifying a single vertex does nothing."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        
        initial_nodes = list(G.nodes())
        identify_vertices(G, [1])
        
        assert list(G.nodes()) == initial_nodes
        assert G.has_edge(1, 2)

    def test_vertex_not_in_graph_raises_error(self) -> None:
        """Test that identifying a non-existent vertex raises ValueError."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        
        with pytest.raises(ValueError, match="not found in graph"):
            identify_vertices(G, [1, 99])

    def test_two_vertices_basic_undirected(self) -> None:
        """Test identifying two vertices with undirected edges."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        
        identify_vertices(G, [1, 2])
        
        assert 2 not in G.nodes()
        assert 1 in G.nodes()
        assert G.has_edge(1, 3)

    def test_two_vertices_basic_directed(self) -> None:
        """Test identifying two vertices with directed edges."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2)
        G.add_directed_edge(2, 3)
        
        identify_vertices(G, [1, 2])
        
        assert 2 not in G.nodes()
        assert 1 in G.nodes()
        assert G.has_edge(1, 3)

    def test_multiple_vertices(self) -> None:
        """Test identifying multiple vertices."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_undirected_edge(3, 4)
        G.add_undirected_edge(4, 5)
        
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

    def test_no_self_loops_created_undirected(self) -> None:
        """Test that self-loops are not created with undirected edges."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 1)  # Same edge, undirected
        
        identify_vertices(G, [1, 2])
        
        # Should not have self-loop
        assert not G.has_edge(1, 1)
        assert 1 in G.nodes()
        assert 2 not in G.nodes()

    def test_no_self_loops_created_directed(self) -> None:
        """Test that self-loops are not created with directed edges."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2)
        G.add_directed_edge(2, 1)  # Bidirectional
        
        identify_vertices(G, [1, 2])
        
        # Should not have self-loop
        assert not G.has_edge(1, 1)
        assert 1 in G.nodes()
        assert 2 not in G.nodes()


class TestIdentifyVerticesMutualExclusivityError:
    """Tests for mutual exclusivity error cases."""

    def test_directed_and_undirected_error(self) -> None:
        """Test that creating both directed and undirected edges raises error."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 3)  # 1 -> 3 (directed)
        G.add_undirected_edge(2, 3)  # 2 - 3 (undirected)
        
        # After merging 1 and 2: would have both directed 1->3 and undirected 1-3
        with pytest.raises(ValueError, match="both directed and undirected"):
            identify_vertices(G, [1, 2])

    def test_directed_and_undirected_error_reverse(self) -> None:
        """Test mutual exclusivity error in reverse direction."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 3)  # 1 - 3 (undirected)
        G.add_directed_edge(2, 3)  # 2 -> 3 (directed)
        
        # After merging 1 and 2: would have both undirected 1-3 and directed 1->3
        with pytest.raises(ValueError, match="both directed and undirected"):
            identify_vertices(G, [1, 2])

    def test_mutual_exclusivity_error_complex(self) -> None:
        """Test mutual exclusivity error with multiple vertices."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 4)  # 1 -> 4
        G.add_undirected_edge(2, 4)  # 2 - 4
        G.add_directed_edge(3, 5)  # 3 -> 5
        
        # After merging [1, 2, 3]: would have both directed 1->4 and undirected 1-4
        with pytest.raises(ValueError, match="both directed and undirected"):
            identify_vertices(G, [1, 2, 3])


class TestIdentifyVerticesAttributes:
    """Tests for node attribute handling."""

    def test_preserve_first_vertex_attributes(self) -> None:
        """Test that first vertex's attributes are preserved by default."""
        G = MixedMultiGraph()
        G.add_node(1)
        G.add_node(2)
        G._undirected.nodes[1]['label'] = 'first'
        G._undirected.nodes[1]['weight'] = 1.0
        G._undirected.nodes[2]['label'] = 'second'
        G._undirected.nodes[2]['weight'] = 2.0
        G.add_undirected_edge(1, 2)
        
        identify_vertices(G, [1, 2])
        
        assert 1 in G.nodes()
        node_attrs = G._undirected.nodes[1]
        assert node_attrs['label'] == 'first'
        assert node_attrs['weight'] == 1.0

    def test_custom_merged_attributes(self) -> None:
        """Test using custom merged attributes."""
        G = MixedMultiGraph()
        G.add_node(1)
        G.add_node(2)
        G.add_node(3)
        G._undirected.nodes[1]['label'] = 'first'
        G._undirected.nodes[2]['label'] = 'second'
        G._undirected.nodes[3]['label'] = 'third'
        G.add_undirected_edge(1, 4)
        G.add_undirected_edge(2, 4)
        G.add_undirected_edge(3, 4)
        
        merged_attrs = {'label': 'merged', 'weight': 5.0, 'custom': 'value'}
        identify_vertices(G, [1, 2, 3], merged_attrs=merged_attrs)
        
        assert 1 in G.nodes()
        node_attrs = G._undirected.nodes[1]
        assert node_attrs['label'] == 'merged'
        assert node_attrs['weight'] == 5.0
        assert node_attrs['custom'] == 'value'


class TestIdentifyVerticesParallelEdges:
    """Tests that identify_vertices may create parallel edges."""

    def test_parallel_edges_created_undirected(self) -> None:
        """Test that parallel undirected edges can be created."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 3, weight=1.0)
        G.add_undirected_edge(2, 3, weight=2.0)
        
        identify_vertices(G, [1, 2])
        
        # Should have parallel undirected edges from 1 to 3
        assert G._undirected.number_of_edges(1, 3) == 2
        assert G.has_edge(1, 3)

    def test_parallel_edges_created_directed(self) -> None:
        """Test that parallel directed edges can be created."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 3, weight=1.0)
        G.add_directed_edge(2, 3, weight=2.0)
        
        identify_vertices(G, [1, 2])
        
        # Should have parallel directed edges from 1 to 3
        assert G._directed.number_of_edges(1, 3) == 2
        assert G.has_edge(1, 3)


class TestIdentifyVerticesComplex:
    """Complex tests for identify_vertices."""

    def test_preserves_other_edges(self) -> None:
        """Test that edges not involving merged vertices are preserved."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_directed_edge(4, 5)
        G.add_directed_edge(5, 6)
        
        identify_vertices(G, [1, 2])
        
        # Edges involving 4, 5, 6 should be unchanged
        assert G.has_edge(4, 5)
        assert G.has_edge(5, 6)
        assert G.has_edge(1, 3)  # New edge from merge

    def test_mixed_edge_types_preserved(self) -> None:
        """Test that different edge types are preserved correctly."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_directed_edge(1, 4)
        G.add_directed_edge(2, 5)
        
        identify_vertices(G, [1, 2])
        
        # Should have both undirected and directed edges from 1
        assert G.has_edge(1, 3)  # Undirected
        assert G.has_edge(1, 4)  # Directed
        assert G.has_edge(1, 5)  # Directed

