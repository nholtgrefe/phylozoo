"""
Tests for MixedMultiGraph isomorphism functions.
"""

import pytest

from phylozoo.core.primitives.m_multigraph import MixedMultiGraph, is_isomorphic


class TestIsIsomorphic:
    """Tests for is_isomorphic function."""
    
    def test_simple_isomorphic_undirected_only(self) -> None:
        """Test that simple isomorphic graphs with only undirected edges are identified correctly."""
        G1 = MixedMultiGraph(undirected_edges=[(1, 2), (2, 3)])
        G2 = MixedMultiGraph(undirected_edges=[(4, 5), (5, 6)])
        
        assert is_isomorphic(G1, G2) is True
    
    def test_simple_isomorphic_directed_only(self) -> None:
        """Test that simple isomorphic graphs with only directed edges are identified correctly."""
        G1 = MixedMultiGraph(directed_edges=[(1, 2), (2, 3)])
        G2 = MixedMultiGraph(directed_edges=[(4, 5), (5, 6)])
        
        assert is_isomorphic(G1, G2) is True
    
    def test_simple_isomorphic_mixed(self) -> None:
        """Test that simple isomorphic graphs with mixed edges are identified correctly."""
        G1 = MixedMultiGraph(
            undirected_edges=[(1, 2)],
            directed_edges=[(2, 3)]
        )
        G2 = MixedMultiGraph(
            undirected_edges=[(4, 5)],
            directed_edges=[(5, 6)]
        )
        
        assert is_isomorphic(G1, G2) is True
    
    def test_simple_non_isomorphic(self) -> None:
        """Test that non-isomorphic graphs are identified correctly."""
        G1 = MixedMultiGraph(undirected_edges=[(1, 2), (2, 3)])  # Path of length 2
        G2 = MixedMultiGraph(undirected_edges=[(4, 5), (4, 6), (5, 6)])  # Triangle (3 edges)
        
        assert is_isomorphic(G1, G2) is False
    
    def test_empty_graphs(self) -> None:
        """Test that empty graphs are isomorphic."""
        G1 = MixedMultiGraph()
        G2 = MixedMultiGraph()
        
        assert is_isomorphic(G1, G2) is True
    
    def test_single_node_graphs(self) -> None:
        """Test that single node graphs are isomorphic."""
        G1 = MixedMultiGraph()
        G1.add_node(1)
        G2 = MixedMultiGraph()
        G2.add_node(2)
        
        assert is_isomorphic(G1, G2) is True
    
    def test_parallel_undirected_edges_isomorphic(self) -> None:
        """Test that graphs with parallel undirected edges are handled correctly."""
        G1 = MixedMultiGraph(undirected_edges=[(1, 2), (1, 2)])
        G2 = MixedMultiGraph(undirected_edges=[(3, 4), (3, 4)])
        
        assert is_isomorphic(G1, G2) is True
    
    def test_parallel_directed_edges_isomorphic(self) -> None:
        """Test that graphs with parallel directed edges are handled correctly."""
        G1 = MixedMultiGraph(directed_edges=[(1, 2), (1, 2)])
        G2 = MixedMultiGraph(directed_edges=[(3, 4), (3, 4)])
        
        assert is_isomorphic(G1, G2) is True
    
    def test_parallel_edges_non_isomorphic(self) -> None:
        """Test that different numbers of parallel edges are not isomorphic."""
        G1 = MixedMultiGraph(undirected_edges=[(1, 2), (1, 2)])
        G2 = MixedMultiGraph(undirected_edges=[(3, 4)])
        
        assert is_isomorphic(G1, G2) is False
    
    def test_undirected_vs_directed_not_isomorphic(self) -> None:
        """Test that undirected edges don't match with directed edges."""
        G1 = MixedMultiGraph(undirected_edges=[(1, 2)])
        G2 = MixedMultiGraph(directed_edges=[(3, 4)])
        
        assert is_isomorphic(G1, G2) is False
    
    def test_self_loops_undirected(self) -> None:
        """Test that self-loops in undirected edges are handled correctly."""
        G1 = MixedMultiGraph(undirected_edges=[(1, 1), (1, 2)])
        G2 = MixedMultiGraph(undirected_edges=[(3, 3), (3, 4)])
        
        assert is_isomorphic(G1, G2) is True
    
    def test_self_loops_directed(self) -> None:
        """Test that self-loops in directed edges are handled correctly."""
        G1 = MixedMultiGraph(directed_edges=[(1, 1), (1, 2)])
        G2 = MixedMultiGraph(directed_edges=[(3, 3), (3, 4)])
        
        assert is_isomorphic(G1, G2) is True
    
    def test_with_matching_node_labels(self) -> None:
        """Test isomorphism with matching node labels."""
        G1 = MixedMultiGraph(undirected_edges=[(1, 2), (2, 3)])
        G1.add_node(1, label='A')
        G1.add_node(2, label='B')
        G1.add_node(3, label='C')
        
        G2 = MixedMultiGraph(undirected_edges=[(4, 5), (5, 6)])
        G2.add_node(4, label='A')
        G2.add_node(5, label='B')
        G2.add_node(6, label='C')
        
        assert is_isomorphic(G1, G2, node_attrs=['label']) is True
    
    def test_with_different_node_labels(self) -> None:
        """Test that graphs with different node labels are not isomorphic."""
        G1 = MixedMultiGraph(undirected_edges=[(1, 2), (2, 3)])
        G1.add_node(1, label='A')
        G1.add_node(2, label='B')
        G1.add_node(3, label='C')
        
        G2 = MixedMultiGraph(undirected_edges=[(4, 5), (5, 6)])
        G2.add_node(4, label='A')
        G2.add_node(5, label='B')
        G2.add_node(6, label='D')  # Different label
        
        assert is_isomorphic(G1, G2, node_attrs=['label']) is False
    
    def test_with_labels_ignored(self) -> None:
        """Test that graphs are isomorphic when labels are ignored."""
        G1 = MixedMultiGraph(undirected_edges=[(1, 2), (2, 3)])
        G1.add_node(1, label='A')
        G1.add_node(2, label='B')
        G1.add_node(3, label='C')
        
        G2 = MixedMultiGraph(undirected_edges=[(4, 5), (5, 6)])
        G2.add_node(4, label='X')
        G2.add_node(5, label='Y')
        G2.add_node(6, label='Z')
        
        # Should be isomorphic when labels are ignored
        assert is_isomorphic(G1, G2) is True
    
    def test_partial_node_labels(self) -> None:
        """Test isomorphism with partial node labels (some nodes labeled, some not)."""
        G1 = MixedMultiGraph(undirected_edges=[(1, 2), (2, 3)])
        G1.add_node(1, label='A')
        # Node 2 and 3 have no labels
        
        G2 = MixedMultiGraph(undirected_edges=[(4, 5), (5, 6)])
        G2.add_node(4, label='A')
        # Node 5 and 6 have no labels
        
        assert is_isomorphic(G1, G2, node_attrs=['label']) is True
    
    def test_multiple_node_attributes(self) -> None:
        """Test isomorphism with multiple node attributes."""
        G1 = MixedMultiGraph(undirected_edges=[(1, 2), (2, 3)])
        G1.add_node(1, label='A', type='leaf')
        G1.add_node(2, label='B', type='internal')
        G1.add_node(3, label='C', type='leaf')
        
        G2 = MixedMultiGraph(undirected_edges=[(4, 5), (5, 6)])
        G2.add_node(4, label='A', type='leaf')
        G2.add_node(5, label='B', type='internal')
        G2.add_node(6, label='C', type='leaf')
        
        assert is_isomorphic(G1, G2, node_attrs=['label', 'type']) is True
    
    def test_edge_attributes_undirected(self) -> None:
        """Test isomorphism with matching edge attributes on undirected edges."""
        G1 = MixedMultiGraph(undirected_edges=[
            {'u': 1, 'v': 2, 'weight': 1.0},
            {'u': 2, 'v': 3, 'weight': 2.0}
        ])
        G2 = MixedMultiGraph(undirected_edges=[
            {'u': 4, 'v': 5, 'weight': 1.0},
            {'u': 5, 'v': 6, 'weight': 2.0}
        ])
        
        assert is_isomorphic(G1, G2, edge_attrs=['weight']) is True
    
    def test_edge_attributes_directed(self) -> None:
        """Test isomorphism with matching edge attributes on directed edges."""
        G1 = MixedMultiGraph(directed_edges=[
            {'u': 1, 'v': 2, 'weight': 1.0},
            {'u': 2, 'v': 3, 'weight': 2.0}
        ])
        G2 = MixedMultiGraph(directed_edges=[
            {'u': 4, 'v': 5, 'weight': 1.0},
            {'u': 5, 'v': 6, 'weight': 2.0}
        ])
        
        assert is_isomorphic(G1, G2, edge_attrs=['weight']) is True
    
    def test_edge_attributes_mismatch(self) -> None:
        """Test that graphs with different edge attributes are not isomorphic."""
        G1 = MixedMultiGraph(undirected_edges=[
            {'u': 1, 'v': 2, 'weight': 1.0}
        ])
        G2 = MixedMultiGraph(undirected_edges=[
            {'u': 3, 'v': 4, 'weight': 2.0}  # Different weight
        ])
        
        assert is_isomorphic(G1, G2, edge_attrs=['weight']) is False
    
    def test_edge_attributes_ignored(self) -> None:
        """Test that graphs are isomorphic when edge attributes are ignored."""
        G1 = MixedMultiGraph(undirected_edges=[
            {'u': 1, 'v': 2, 'weight': 1.0}
        ])
        G2 = MixedMultiGraph(undirected_edges=[
            {'u': 3, 'v': 4, 'weight': 2.0}
        ])
        
        # Should be isomorphic when edge attributes are ignored
        assert is_isomorphic(G1, G2) is True
    
    def test_multiple_edge_attributes(self) -> None:
        """Test isomorphism with multiple edge attributes."""
        G1 = MixedMultiGraph(undirected_edges=[
            {'u': 1, 'v': 2, 'weight': 1.0, 'bootstrap': 0.95}
        ])
        G2 = MixedMultiGraph(undirected_edges=[
            {'u': 3, 'v': 4, 'weight': 1.0, 'bootstrap': 0.95}
        ])
        
        assert is_isomorphic(G1, G2, edge_attrs=['weight', 'bootstrap']) is True
    
    def test_graph_attributes_matching(self) -> None:
        """Test isomorphism with matching graph attributes."""
        G1 = MixedMultiGraph(attributes={'version': '1.0', 'source': 'test'})
        G2 = MixedMultiGraph(attributes={'version': '1.0', 'source': 'test'})
        
        assert is_isomorphic(G1, G2, graph_attrs=['version', 'source']) is True
    
    def test_graph_attributes_mismatch(self) -> None:
        """Test that graphs with different graph attributes are not isomorphic."""
        G1 = MixedMultiGraph(attributes={'version': '1.0'})
        G2 = MixedMultiGraph(attributes={'version': '2.0'})
        
        assert is_isomorphic(G1, G2, graph_attrs=['version']) is False
    
    def test_combined_attributes(self) -> None:
        """Test isomorphism with node, edge, and graph attributes combined."""
        G1 = MixedMultiGraph(
            undirected_edges=[{'u': 1, 'v': 2, 'weight': 1.0}],
            attributes={'version': '1.0'}
        )
        G1.add_node(1, label='A')
        G1.add_node(2, label='B')
        
        G2 = MixedMultiGraph(
            undirected_edges=[{'u': 3, 'v': 4, 'weight': 1.0}],
            attributes={'version': '1.0'}
        )
        G2.add_node(3, label='A')
        G2.add_node(4, label='B')
        
        assert is_isomorphic(
            G1, G2,
            node_attrs=['label'],
            edge_attrs=['weight'],
            graph_attrs=['version']
        ) is True
    
    def test_parallel_undirected_edges_with_attributes(self) -> None:
        """Test that parallel undirected edges with attributes are handled correctly."""
        G1 = MixedMultiGraph(undirected_edges=[
            {'u': 1, 'v': 2, 'weight': 1.0, 'key': 0},
            {'u': 1, 'v': 2, 'weight': 2.0, 'key': 1}
        ])
        G2 = MixedMultiGraph(undirected_edges=[
            {'u': 3, 'v': 4, 'weight': 1.0, 'key': 0},
            {'u': 3, 'v': 4, 'weight': 2.0, 'key': 1}
        ])
        
        assert is_isomorphic(G1, G2, edge_attrs=['weight']) is True
    
    def test_mixed_directed_undirected_complex(self) -> None:
        """Test isomorphism with complex mixed directed and undirected edges."""
        G1 = MixedMultiGraph(
            undirected_edges=[(1, 2), (2, 3)],
            directed_edges=[(3, 4), (4, 5), (5, 1)]
        )
        G2 = MixedMultiGraph(
            undirected_edges=[(6, 7), (7, 8)],
            directed_edges=[(8, 9), (9, 10), (10, 6)]
        )
        
        assert is_isomorphic(G1, G2) is True
    
    def test_mixed_directed_undirected_not_isomorphic(self) -> None:
        """Test that mixed graphs with different structures are not isomorphic."""
        G1 = MixedMultiGraph(
            undirected_edges=[(1, 2)],
            directed_edges=[(2, 3)]
        )
        G2 = MixedMultiGraph(
            undirected_edges=[(4, 5)],
            directed_edges=[(5, 6), (6, 7)]  # Extra edge
        )
        
        assert is_isomorphic(G1, G2) is False
    
    def test_complex_structure(self) -> None:
        """Test isomorphism with more complex graph structures."""
        G1 = MixedMultiGraph(undirected_edges=[
            (1, 2), (2, 3), (2, 4), (3, 5), (4, 5)
        ])
        G2 = MixedMultiGraph(undirected_edges=[
            (6, 7), (7, 8), (7, 9), (8, 10), (9, 10)
        ])
        
        assert is_isomorphic(G1, G2) is True
    
    def test_same_graph(self) -> None:
        """Test that a graph is isomorphic to itself."""
        G1 = MixedMultiGraph(
            undirected_edges=[(1, 2)],
            directed_edges=[(2, 3)]
        )
        G1.add_node(1, label='A')
        
        assert is_isomorphic(G1, G1) is True
        assert is_isomorphic(G1, G1, node_attrs=['label']) is True
    
    def test_bidirectional_undirected_edges(self) -> None:
        """Test that undirected edges are correctly converted to bidirectional."""
        # This test verifies that undirected edges become bidirectional in conversion
        G1 = MixedMultiGraph(undirected_edges=[(1, 2)])
        G2 = MixedMultiGraph(undirected_edges=[(3, 4)])
        
        # Should be isomorphic because both undirected edges become bidirectional
        assert is_isomorphic(G1, G2) is True
        
        # But undirected should NOT match with a single directed edge
        G3 = MixedMultiGraph(directed_edges=[(5, 6)])
        assert is_isomorphic(G1, G3) is False

