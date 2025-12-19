"""
Tests for identify_parallel_edge function in d_multigraph transformations.
"""

import pytest

from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
from phylozoo.core.primitives.d_multigraph.transformations import identify_parallel_edge


class TestIdentifyParallelEdgeBasic:
    """Basic tests for identify_parallel_edge."""

    def test_single_edge_no_change(self) -> None:
        """Test that identifying parallel edges on a single edge does nothing."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0)
        
        initial_edges = G.number_of_edges()
        identify_parallel_edge(G, 1, 2)
        
        assert G.number_of_edges() == initial_edges
        assert G.has_edge(1, 2)
        assert G._graph[1][2][0]['weight'] == 1.0

    def test_two_parallel_edges(self) -> None:
        """Test identifying two parallel edges."""
        G = DirectedMultiGraph()
        key1 = G.add_edge(1, 2, weight=1.0)
        key2 = G.add_edge(1, 2, weight=2.0)
        
        assert G._graph.number_of_edges(1, 2) == 2
        
        identify_parallel_edge(G, 1, 2)
        
        assert G._graph.number_of_edges(1, 2) == 1
        assert G.has_edge(1, 2)
        # First edge (lowest key) should be kept
        assert G._graph[1][2][key1]['weight'] == 2.0  # Merged: last value wins

    def test_three_parallel_edges(self) -> None:
        """Test identifying three parallel edges."""
        G = DirectedMultiGraph()
        key1 = G.add_edge(1, 2, weight=1.0, label='first')
        key2 = G.add_edge(1, 2, weight=2.0, label='second')
        key3 = G.add_edge(1, 2, weight=3.0, label='third')
        
        assert G._graph.number_of_edges(1, 2) == 3
        
        identify_parallel_edge(G, 1, 2)
        
        assert G._graph.number_of_edges(1, 2) == 1
        assert G.has_edge(1, 2)
        # First edge key should be kept with merged attributes
        edge_data = G._graph[1][2][key1]
        assert edge_data['weight'] == 3.0  # Last value wins
        assert edge_data['label'] == 'third'  # Last value wins

    def test_parallel_edges_with_different_attributes(self) -> None:
        """Test identifying parallel edges with different attribute sets."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0, color='red')
        G.add_edge(1, 2, size=5.0, color='blue')
        G.add_edge(1, 2, weight=3.0, label='test')
        
        identify_parallel_edge(G, 1, 2)
        
        assert G._graph.number_of_edges(1, 2) == 1
        edge_data = G._graph[1][2][0]
        # All attributes should be present, last values win
        assert edge_data['weight'] == 3.0
        assert edge_data['color'] == 'blue'
        assert edge_data['size'] == 5.0
        assert edge_data['label'] == 'test'


class TestIdentifyParallelEdgeWithMergedAttrs:
    """Tests for identify_parallel_edge with merged_attrs parameter."""

    def test_custom_merged_attributes(self) -> None:
        """Test using custom merged attributes."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0, label='first')
        G.add_edge(1, 2, weight=2.0, label='second')
        
        merged_attrs = {'weight': 5.0, 'merged': True, 'custom': 'value'}
        identify_parallel_edge(G, 1, 2, merged_attrs=merged_attrs)
        
        assert G._graph.number_of_edges(1, 2) == 1
        edge_data = G._graph[1][2][0]
        assert edge_data['weight'] == 5.0
        assert edge_data['merged'] is True
        assert edge_data['custom'] == 'value'
        # Original attributes should be replaced
        assert 'label' not in edge_data

    def test_empty_merged_attributes(self) -> None:
        """Test with empty merged attributes dict."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0, label='test')
        G.add_edge(1, 2, weight=2.0)
        
        identify_parallel_edge(G, 1, 2, merged_attrs={})
        
        assert G._graph.number_of_edges(1, 2) == 1
        edge_data = G._graph[1][2][0]
        # Should have no attributes (empty dict)
        assert len(edge_data) == 0 or edge_data == {}


class TestIdentifyParallelEdgeErrorCases:
    """Tests for error cases."""

    def test_node_not_in_graph_raises_error(self) -> None:
        """Test that identifying edges with non-existent node raises ValueError."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        
        with pytest.raises(ValueError, match="not found in graph"):
            identify_parallel_edge(G, 99, 2)
        
        with pytest.raises(ValueError, match="not found in graph"):
            identify_parallel_edge(G, 1, 99)

    def test_no_edges_between_nodes_raises_error(self) -> None:
        """Test that identifying edges between nodes with no edges raises ValueError."""
        G = DirectedMultiGraph()
        G.add_node(1)
        G.add_node(2)
        G.add_edge(1, 3)  # Edge exists, but not between 1 and 2
        
        with pytest.raises(ValueError, match="No edges exist between nodes"):
            identify_parallel_edge(G, 1, 2)


class TestIdentifyParallelEdgeComplex:
    """Complex tests for identify_parallel_edge."""

    def test_identify_parallel_edges_preserves_other_edges(self) -> None:
        """Test that identifying parallel edges doesn't affect other edges."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0)
        G.add_edge(1, 2, weight=2.0)
        G.add_edge(2, 3, weight=3.0)
        G.add_edge(3, 4, weight=4.0)
        G.add_edge(1, 4, weight=5.0)
        
        initial_other_edges = {
            (2, 3): G._graph[2][3][0].copy(),
            (3, 4): G._graph[3][4][0].copy(),
            (1, 4): G._graph[1][4][0].copy(),
        }
        
        identify_parallel_edge(G, 1, 2)
        
        # Other edges should be unchanged
        assert G.has_edge(2, 3)
        assert G.has_edge(3, 4)
        assert G.has_edge(1, 4)
        assert G._graph[2][3][0]['weight'] == 3.0
        assert G._graph[3][4][0]['weight'] == 4.0
        assert G._graph[1][4][0]['weight'] == 5.0

    def test_identify_parallel_edges_multiple_pairs(self) -> None:
        """Test identifying parallel edges in multiple edge pairs."""
        G = DirectedMultiGraph()
        # Multiple parallel edges between 1 and 2
        G.add_edge(1, 2, weight=1.0)
        G.add_edge(1, 2, weight=2.0)
        # Multiple parallel edges between 3 and 4
        G.add_edge(3, 4, weight=3.0)
        G.add_edge(3, 4, weight=4.0)
        G.add_edge(3, 4, weight=5.0)
        
        identify_parallel_edge(G, 1, 2)
        identify_parallel_edge(G, 3, 4)
        
        assert G._graph.number_of_edges(1, 2) == 1
        assert G._graph.number_of_edges(3, 4) == 1
        assert G.number_of_edges() == 2

    def test_identify_parallel_edges_with_none_values(self) -> None:
        """Test identifying parallel edges when some attributes are missing."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0)
        G.add_edge(1, 2, weight=2.0, label='test')
        
        identify_parallel_edge(G, 1, 2)
        
        assert G._graph.number_of_edges(1, 2) == 1
        edge_data = G._graph[1][2][0]
        assert edge_data['weight'] == 2.0
        assert edge_data['label'] == 'test'  # Missing attribute should be added from second edge

    def test_identify_parallel_edges_preserves_first_key(self) -> None:
        """Test that the first (lowest) key is preserved."""
        G = DirectedMultiGraph()
        key1 = G.add_edge(1, 2, weight=1.0)
        key2 = G.add_edge(1, 2, weight=2.0)
        key3 = G.add_edge(1, 2, weight=3.0)
        
        # Verify keys are sequential
        assert key1 < key2 < key3
        
        identify_parallel_edge(G, 1, 2)
        
        # First key should be kept
        assert G.has_edge(1, 2, key=key1)
        assert not G.has_edge(1, 2, key=key2)
        assert not G.has_edge(1, 2, key=key3)

