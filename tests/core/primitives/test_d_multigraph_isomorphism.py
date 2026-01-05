"""
Tests for DirectedMultiGraph isomorphism functions.
"""

import pytest

from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph, is_isomorphic


class TestIsIsomorphic:
    """Tests for is_isomorphic function."""
    
    def test_simple_isomorphic_graphs(self) -> None:
        """Test that simple isomorphic graphs are identified correctly."""
        G1 = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
        G2 = DirectedMultiGraph(edges=[(4, 5), (5, 6)])
        
        assert is_isomorphic(G1, G2) is True
    
    def test_simple_non_isomorphic_graphs(self) -> None:
        """Test that non-isomorphic graphs are identified correctly."""
        G1 = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
        G2 = DirectedMultiGraph(edges=[(4, 5), (6, 5)])  # Different structure
        
        assert is_isomorphic(G1, G2) is False
    
    def test_empty_graphs(self) -> None:
        """Test that empty graphs are isomorphic."""
        G1 = DirectedMultiGraph()
        G2 = DirectedMultiGraph()
        
        assert is_isomorphic(G1, G2) is True
    
    def test_single_node_graphs(self) -> None:
        """Test that single node graphs are isomorphic."""
        G1 = DirectedMultiGraph()
        G1.add_node(1)
        G2 = DirectedMultiGraph()
        G2.add_node(2)
        
        assert is_isomorphic(G1, G2) is True
    
    def test_parallel_edges_isomorphic(self) -> None:
        """Test that graphs with parallel edges are handled correctly."""
        G1 = DirectedMultiGraph(edges=[(1, 2), (1, 2)])
        G2 = DirectedMultiGraph(edges=[(3, 4), (3, 4)])
        
        assert is_isomorphic(G1, G2) is True
    
    def test_parallel_edges_non_isomorphic(self) -> None:
        """Test that different numbers of parallel edges are not isomorphic."""
        G1 = DirectedMultiGraph(edges=[(1, 2), (1, 2)])
        G2 = DirectedMultiGraph(edges=[(3, 4)])
        
        assert is_isomorphic(G1, G2) is False
    
    def test_self_loops(self) -> None:
        """Test that self-loops are handled correctly."""
        G1 = DirectedMultiGraph(edges=[(1, 1), (1, 2)])
        G2 = DirectedMultiGraph(edges=[(3, 3), (3, 4)])
        
        assert is_isomorphic(G1, G2) is True
    
    def test_with_matching_labels(self) -> None:
        """Test isomorphism with matching node labels."""
        G1 = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
        G1.add_node(1, label='A')
        G1.add_node(2, label='B')
        G1.add_node(3, label='C')
        
        G2 = DirectedMultiGraph(edges=[(4, 5), (5, 6)])
        G2.add_node(4, label='A')
        G2.add_node(5, label='B')
        G2.add_node(6, label='C')
        
        assert is_isomorphic(G1, G2, node_attrs=['label']) is True
    
    def test_with_different_labels(self) -> None:
        """Test that graphs with different labels are not isomorphic."""
        G1 = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
        G1.add_node(1, label='A')
        G1.add_node(2, label='B')
        G1.add_node(3, label='C')
        
        G2 = DirectedMultiGraph(edges=[(4, 5), (5, 6)])
        G2.add_node(4, label='A')
        G2.add_node(5, label='B')
        G2.add_node(6, label='D')  # Different label
        
        assert is_isomorphic(G1, G2, node_attrs=['label']) is False
    
    def test_with_labels_ignored(self) -> None:
        """Test that graphs are isomorphic when labels are ignored."""
        G1 = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
        G1.add_node(1, label='A')
        G1.add_node(2, label='B')
        G1.add_node(3, label='C')
        
        G2 = DirectedMultiGraph(edges=[(4, 5), (5, 6)])
        G2.add_node(4, label='X')
        G2.add_node(5, label='Y')
        G2.add_node(6, label='Z')
        
        # Should be isomorphic when labels are ignored
        assert is_isomorphic(G1, G2) is True
    
    def test_partial_labels(self) -> None:
        """Test isomorphism with partial labels (some nodes labeled, some not)."""
        G1 = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
        G1.add_node(1, label='A')
        # Node 2 and 3 have no labels
        
        G2 = DirectedMultiGraph(edges=[(4, 5), (5, 6)])
        G2.add_node(4, label='A')
        # Node 5 and 6 have no labels
        
        assert is_isomorphic(G1, G2, node_attrs=['label']) is True
    
    def test_partial_labels_mismatch(self) -> None:
        """Test that partial labels must match."""
        G1 = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
        G1.add_node(1, label='A')
        
        G2 = DirectedMultiGraph(edges=[(4, 5), (5, 6)])
        G2.add_node(4, label='B')  # Different label
        
        assert is_isomorphic(G1, G2, node_attrs=['label']) is False
    
    def test_custom_label_attribute(self) -> None:
        """Test isomorphism with custom label attribute name."""
        G1 = DirectedMultiGraph(edges=[(1, 2)])
        G1.add_node(1, taxon='A')
        G1.add_node(2, taxon='B')
        
        G2 = DirectedMultiGraph(edges=[(3, 4)])
        G2.add_node(3, taxon='A')
        G2.add_node(4, taxon='B')
        
        assert is_isomorphic(G1, G2, node_attrs=['taxon']) is True
    
    def test_multiple_node_attributes(self) -> None:
        """Test isomorphism with multiple node attributes."""
        G1 = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
        G1.add_node(1, label='A', type='leaf')
        G1.add_node(2, label='B', type='internal')
        G1.add_node(3, label='C', type='leaf')
        
        G2 = DirectedMultiGraph(edges=[(4, 5), (5, 6)])
        G2.add_node(4, label='A', type='leaf')
        G2.add_node(5, label='B', type='internal')
        G2.add_node(6, label='C', type='leaf')
        
        assert is_isomorphic(G1, G2, node_attrs=['label', 'type']) is True
    
    def test_multiple_node_attributes_mismatch(self) -> None:
        """Test that all node attributes must match."""
        G1 = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
        G1.add_node(1, label='A', type='leaf')
        G1.add_node(2, label='B', type='internal')
        
        G2 = DirectedMultiGraph(edges=[(4, 5), (5, 6)])
        G2.add_node(4, label='A', type='leaf')
        G2.add_node(5, label='B', type='hybrid')  # Different type
        
        assert is_isomorphic(G1, G2, node_attrs=['label', 'type']) is False
    
    def test_edge_attributes_matching(self) -> None:
        """Test isomorphism with matching edge attributes."""
        G1 = DirectedMultiGraph(edges=[
            {'u': 1, 'v': 2, 'weight': 1.0},
            {'u': 2, 'v': 3, 'weight': 2.0}
        ])
        G2 = DirectedMultiGraph(edges=[
            {'u': 4, 'v': 5, 'weight': 1.0},
            {'u': 5, 'v': 6, 'weight': 2.0}
        ])
        
        assert is_isomorphic(G1, G2, edge_attrs=['weight']) is True
    
    def test_edge_attributes_mismatch(self) -> None:
        """Test that graphs with different edge attributes are not isomorphic."""
        G1 = DirectedMultiGraph(edges=[
            {'u': 1, 'v': 2, 'weight': 1.0}
        ])
        G2 = DirectedMultiGraph(edges=[
            {'u': 3, 'v': 4, 'weight': 2.0}  # Different weight
        ])
        
        assert is_isomorphic(G1, G2, edge_attrs=['weight']) is False
    
    def test_edge_attributes_ignored(self) -> None:
        """Test that graphs are isomorphic when edge attributes are ignored."""
        G1 = DirectedMultiGraph(edges=[
            {'u': 1, 'v': 2, 'weight': 1.0}
        ])
        G2 = DirectedMultiGraph(edges=[
            {'u': 3, 'v': 4, 'weight': 2.0}
        ])
        
        # Should be isomorphic when edge attributes are ignored
        assert is_isomorphic(G1, G2) is True
    
    def test_multiple_edge_attributes(self) -> None:
        """Test isomorphism with multiple edge attributes."""
        G1 = DirectedMultiGraph(edges=[
            {'u': 1, 'v': 2, 'weight': 1.0, 'bootstrap': 0.95}
        ])
        G2 = DirectedMultiGraph(edges=[
            {'u': 3, 'v': 4, 'weight': 1.0, 'bootstrap': 0.95}
        ])
        
        assert is_isomorphic(G1, G2, edge_attrs=['weight', 'bootstrap']) is True
    
    def test_graph_attributes_matching(self) -> None:
        """Test isomorphism with matching graph attributes."""
        G1 = DirectedMultiGraph(attributes={'version': '1.0', 'source': 'test'})
        G2 = DirectedMultiGraph(attributes={'version': '1.0', 'source': 'test'})
        
        assert is_isomorphic(G1, G2, graph_attrs=['version', 'source']) is True
    
    def test_graph_attributes_mismatch(self) -> None:
        """Test that graphs with different graph attributes are not isomorphic."""
        G1 = DirectedMultiGraph(attributes={'version': '1.0'})
        G2 = DirectedMultiGraph(attributes={'version': '2.0'})
        
        assert is_isomorphic(G1, G2, graph_attrs=['version']) is False
    
    def test_graph_attributes_ignored(self) -> None:
        """Test that graphs are isomorphic when graph attributes are ignored."""
        G1 = DirectedMultiGraph(attributes={'version': '1.0'})
        G2 = DirectedMultiGraph(attributes={'version': '2.0'})
        
        # Should be isomorphic when graph attributes are ignored
        assert is_isomorphic(G1, G2) is True
    
    def test_combined_attributes(self) -> None:
        """Test isomorphism with node, edge, and graph attributes combined."""
        G1 = DirectedMultiGraph(
            edges=[{'u': 1, 'v': 2, 'weight': 1.0}],
            attributes={'version': '1.0'}
        )
        G1.add_node(1, label='A')
        G1.add_node(2, label='B')
        
        G2 = DirectedMultiGraph(
            edges=[{'u': 3, 'v': 4, 'weight': 1.0}],
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
    
    def test_combined_attributes_partial_mismatch(self) -> None:
        """Test that any mismatched attribute prevents isomorphism."""
        G1 = DirectedMultiGraph(
            edges=[{'u': 1, 'v': 2, 'weight': 1.0}],
            attributes={'version': '1.0'}
        )
        G1.add_node(1, label='A')
        
        G2 = DirectedMultiGraph(
            edges=[{'u': 3, 'v': 4, 'weight': 1.0}],
            attributes={'version': '2.0'}  # Different version
        )
        G2.add_node(3, label='A')
        
        assert is_isomorphic(
            G1, G2,
            node_attrs=['label'],
            edge_attrs=['weight'],
            graph_attrs=['version']
        ) is False
    
    def test_parallel_edges_with_attributes(self) -> None:
        """Test that parallel edges with attributes are handled correctly."""
        G1 = DirectedMultiGraph(edges=[
            {'u': 1, 'v': 2, 'weight': 1.0, 'key': 0},
            {'u': 1, 'v': 2, 'weight': 2.0, 'key': 1}
        ])
        G2 = DirectedMultiGraph(edges=[
            {'u': 3, 'v': 4, 'weight': 1.0, 'key': 0},
            {'u': 3, 'v': 4, 'weight': 2.0, 'key': 1}
        ])
        
        assert is_isomorphic(G1, G2, edge_attrs=['weight']) is True
    
    def test_parallel_edges_different_attributes(self) -> None:
        """Test that parallel edges must have matching attributes."""
        G1 = DirectedMultiGraph(edges=[
            {'u': 1, 'v': 2, 'weight': 1.0},
            {'u': 1, 'v': 2, 'weight': 2.0}
        ])
        G2 = DirectedMultiGraph(edges=[
            {'u': 3, 'v': 4, 'weight': 1.0},
            {'u': 3, 'v': 4, 'weight': 3.0}  # Different weight
        ])
        
        assert is_isomorphic(G1, G2, edge_attrs=['weight']) is False
    
    def test_complex_structure(self) -> None:
        """Test isomorphism with more complex graph structures."""
        G1 = DirectedMultiGraph(edges=[
            (1, 2), (2, 3), (2, 4), (3, 5), (4, 5)
        ])
        G2 = DirectedMultiGraph(edges=[
            (6, 7), (7, 8), (7, 9), (8, 10), (9, 10)
        ])
        
        assert is_isomorphic(G1, G2) is True
    
    def test_different_degrees(self) -> None:
        """Test that graphs with different node degrees are not isomorphic."""
        G1 = DirectedMultiGraph(edges=[(1, 2), (1, 3), (1, 4)])
        G2 = DirectedMultiGraph(edges=[(5, 6), (5, 7)])
        
        assert is_isomorphic(G1, G2) is False
    
    def test_same_graph(self) -> None:
        """Test that a graph is isomorphic to itself."""
        G1 = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
        G1.add_node(1, label='A')
        
        assert is_isomorphic(G1, G1) is True
        assert is_isomorphic(G1, G1, node_attrs=['label']) is True

