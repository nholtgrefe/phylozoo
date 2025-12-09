"""
Tests for orient_mixed_graph_from_root function.

This test suite provides comprehensive coverage of the orient_mixed_graph_from_root
function, including basic functionality, edge cases, error handling, and complex structures.
"""

import pytest
from typing import Dict, List, Set, Tuple

from phylozoo.core.primitives.m_multigraph import (
    MixedMultiGraph,
    orient_mixed_graph_from_root,
)
from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph


class TestOrientMixedGraphFromRoot:
    """Test cases for orient_mixed_graph_from_root function."""
    
    def test_simple_undirected_tree(self) -> None:
        """Test orienting a simple undirected tree."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(1, 3)
        G.add_undirected_edge(2, 4)
        G.add_undirected_edge(2, 5)
        
        dm = orient_mixed_graph_from_root(G, 1)
        
        # All edges should be directed away from root 1
        assert dm.has_edge(1, 2)
        assert dm.has_edge(1, 3)
        assert dm.has_edge(2, 4)
        assert dm.has_edge(2, 5)
        assert dm.number_of_edges() == G.number_of_edges()
        assert dm.number_of_nodes() == G.number_of_nodes()
    
    def test_mixed_directed_undirected(self) -> None:
        """Test orienting graph with both directed and undirected edges."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_undirected_edge(2, 4)
        G.add_directed_edge(3, 5)
        G.add_undirected_edge(4, 6)
        
        dm = orient_mixed_graph_from_root(G, 1)
        
        # Directed edges should be preserved
        assert dm.has_edge(1, 2)
        assert dm.has_edge(3, 5)
        # Undirected edges should be oriented
        assert dm.has_edge(2, 3)
        assert dm.has_edge(2, 4)
        assert dm.has_edge(4, 6)
        assert dm.number_of_edges() == G.number_of_edges()
    
    def test_parallel_undirected_edges(self) -> None:
        """Test orienting graph with parallel undirected edges."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(1, 2)  # Parallel
        G.add_undirected_edge(1, 2)  # Parallel
        G.add_undirected_edge(2, 3)
        
        dm = orient_mixed_graph_from_root(G, 1)
        
        # All parallel edges should be oriented
        assert dm.number_of_edges() == G.number_of_edges()
        # Check that all edges from 1 to 2 are present
        edges_1_2 = [e for e in dm.edges(keys=True) if e[0] == 1 and e[1] == 2]
        assert len(edges_1_2) == 3
        assert dm.has_edge(2, 3)
    
    def test_parallel_directed_edges(self) -> None:
        """Test orienting graph with parallel directed edges."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2)
        G.add_directed_edge(1, 2)  # Parallel
        G.add_undirected_edge(2, 3)
        G.add_undirected_edge(2, 4)
        
        dm = orient_mixed_graph_from_root(G, 1)
        
        # All edges should be preserved
        assert dm.number_of_edges() == G.number_of_edges()
        edges_1_2 = [e for e in dm.edges(keys=True) if e[0] == 1 and e[1] == 2]
        assert len(edges_1_2) == 2
        assert dm.has_edge(2, 3)
        assert dm.has_edge(2, 4)
    
    def test_hybrid_node_structure(self) -> None:
        """Test orienting graph with hybrid node (multiple incoming directed edges)."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(1, 3)
        G.add_directed_edge(2, 4)
        G.add_directed_edge(3, 4)  # Hybrid node 4
        G.add_undirected_edge(4, 5)
        G.add_undirected_edge(4, 6)
        
        dm = orient_mixed_graph_from_root(G, 1)
        
        # All edges should be oriented
        assert dm.has_edge(1, 2)
        assert dm.has_edge(1, 3)
        assert dm.has_edge(2, 4)
        assert dm.has_edge(3, 4)
        assert dm.has_edge(4, 5)
        assert dm.has_edge(4, 6)
        assert dm.number_of_edges() == G.number_of_edges()
    
    def test_large_tree_structure(self) -> None:
        """Test orienting a larger tree structure."""
        G = MixedMultiGraph()
        # Create a binary tree
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(1, 3)
        G.add_undirected_edge(2, 4)
        G.add_undirected_edge(2, 5)
        G.add_undirected_edge(3, 6)
        G.add_undirected_edge(3, 7)
        G.add_undirected_edge(4, 8)
        G.add_undirected_edge(5, 9)
        G.add_undirected_edge(6, 10)
        G.add_undirected_edge(7, 11)
        
        dm = orient_mixed_graph_from_root(G, 1)
        
        # Verify all edges are oriented away from root
        assert dm.number_of_edges() == G.number_of_edges()
        assert dm.number_of_nodes() == G.number_of_nodes()
        # Check some specific edges
        assert dm.has_edge(1, 2)
        assert dm.has_edge(1, 3)
        assert dm.has_edge(2, 4)
        assert dm.has_edge(2, 5)
    
    def test_root_not_in_graph_raises_error(self) -> None:
        """Test that root not in graph raises ValueError."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        
        with pytest.raises(ValueError, match="Root vertex.*not found"):
            orient_mixed_graph_from_root(G, 999)
    
    def test_directed_edge_pointing_wrong_way_raises_error(self) -> None:
        """Test that directed edge pointing towards BFS path raises error."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_directed_edge(4, 2)  # Points towards 2, which is in BFS path
        
        with pytest.raises(ValueError, match="points towards vertex"):
            orient_mixed_graph_from_root(G, 1)
    
    def test_directed_cycle_raises_error(self) -> None:
        """Test that directed cycle raises error."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_directed_edge(2, 3)
        G.add_directed_edge(3, 4)
        G.add_directed_edge(4, 2)  # Creates cycle 2->3->4->2
        
        with pytest.raises(ValueError, match="points towards vertex"):
            orient_mixed_graph_from_root(G, 1)
    
    def test_unreachable_undirected_edge_raises_error(self) -> None:
        """Test that unreachable undirected edge raises error."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(3, 4)  # Disconnected component
        
        with pytest.raises(ValueError, match="not reachable from root"):
            orient_mixed_graph_from_root(G, 1)
    
    def test_single_node_graph(self) -> None:
        """Test orienting graph with single node."""
        G = MixedMultiGraph()
        G.add_node(1)
        
        dm = orient_mixed_graph_from_root(G, 1)
        
        assert dm.number_of_nodes() == 1
        assert dm.number_of_edges() == 0
    
    def test_single_edge_graph(self) -> None:
        """Test orienting graph with single edge."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        
        dm = orient_mixed_graph_from_root(G, 1)
        
        assert dm.has_edge(1, 2)
        assert dm.number_of_edges() == 1
        assert dm.number_of_nodes() == 2
    
    def test_preserves_edge_attributes(self) -> None:
        """Test that edge attributes are preserved during orientation."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2, weight=1.5, label='test')
        G.add_undirected_edge(2, 3, weight=2.0)
        G.add_directed_edge(3, 4, weight=3.0)
        
        dm = orient_mixed_graph_from_root(G, 1)
        
        # Check attributes are preserved
        edge_data_1_2 = dict(dm._graph[1][2][0])
        assert edge_data_1_2.get('weight') == 1.5
        assert edge_data_1_2.get('label') == 'test'
        
        edge_data_2_3 = dict(dm._graph[2][3][0])
        assert edge_data_2_3.get('weight') == 2.0
        
        edge_data_3_4 = dict(dm._graph[3][4][0])
        assert edge_data_3_4.get('weight') == 3.0
    
    def test_preserves_node_attributes(self) -> None:
        """Test that node attributes are preserved during orientation."""
        G = MixedMultiGraph()
        G.add_node(1, label='root', color='red')
        G.add_node(2, label='leaf', color='blue')
        G.add_undirected_edge(1, 2)
        
        dm = orient_mixed_graph_from_root(G, 1)
        
        # Access node data using nodes() method
        node_data_1 = dict(dm.nodes(data=True))[1]
        node_data_2 = dict(dm.nodes(data=True))[2]
        assert node_data_1.get('label') == 'root'
        assert node_data_1.get('color') == 'red'
        assert node_data_2.get('label') == 'leaf'
        assert node_data_2.get('color') == 'blue'
    
    def test_complex_mixed_structure(self) -> None:
        """Test orienting a complex mixed structure."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(1, 3)
        G.add_directed_edge(2, 4)
        G.add_directed_edge(3, 4)  # Hybrid node
        G.add_undirected_edge(4, 5)
        G.add_undirected_edge(4, 6)
        G.add_directed_edge(5, 7)
        G.add_undirected_edge(6, 8)
        G.add_undirected_edge(7, 9)
        G.add_undirected_edge(8, 10)
        
        dm = orient_mixed_graph_from_root(G, 1)
        
        # Verify all edges are oriented
        assert dm.number_of_edges() == G.number_of_edges()
        assert dm.number_of_nodes() == G.number_of_nodes()
        # Check key edges
        assert dm.has_edge(1, 2)
        assert dm.has_edge(1, 3)
        assert dm.has_edge(2, 4)
        assert dm.has_edge(3, 4)
        assert dm.has_edge(4, 5)
        assert dm.has_edge(4, 6)
        assert dm.has_edge(5, 7)
        assert dm.has_edge(6, 8)
        assert dm.has_edge(7, 9)
        assert dm.has_edge(8, 10)
    
    def test_multiple_directed_paths(self) -> None:
        """Test graph with multiple directed paths from root."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2)
        G.add_directed_edge(1, 3)
        G.add_undirected_edge(2, 4)
        G.add_undirected_edge(2, 5)
        G.add_undirected_edge(3, 6)
        G.add_undirected_edge(3, 7)
        
        dm = orient_mixed_graph_from_root(G, 1)
        
        assert dm.has_edge(1, 2)
        assert dm.has_edge(1, 3)
        assert dm.has_edge(2, 4)
        assert dm.has_edge(2, 5)
        assert dm.has_edge(3, 6)
        assert dm.has_edge(3, 7)
        assert dm.number_of_edges() == G.number_of_edges()
    
    def test_all_directed_edges_already_oriented(self) -> None:
        """Test graph where all edges are already directed correctly."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2)
        G.add_directed_edge(1, 3)
        G.add_directed_edge(2, 4)
        G.add_directed_edge(3, 5)
        
        dm = orient_mixed_graph_from_root(G, 1)
        
        # Should preserve all directed edges
        assert dm.number_of_edges() == G.number_of_edges()
        assert dm.has_edge(1, 2)
        assert dm.has_edge(1, 3)
        assert dm.has_edge(2, 4)
        assert dm.has_edge(3, 5)
    
    def test_orientation_from_different_roots(self) -> None:
        """Test that orientation differs when using different roots."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_undirected_edge(2, 4)
        
        dm1 = orient_mixed_graph_from_root(G, 1)
        dm2 = orient_mixed_graph_from_root(G, 3)
        
        # Both should have all edges directed
        assert dm1.number_of_edges() == G.number_of_edges()
        assert dm2.number_of_edges() == G.number_of_edges()
        
        # But orientations should differ
        assert dm1.has_edge(1, 2)
        assert dm1.has_edge(2, 3)
        assert dm1.has_edge(2, 4)
        
        # From root 3, edges should point AWAY from 3 (not towards)
        assert dm2.has_edge(3, 2)  # 3->2 (away from root)
        assert dm2.has_edge(2, 1)  # 2->1 (away from root via 2)
        assert dm2.has_edge(2, 4)  # 2->4 (away from root via 2)
    
    def test_parallel_edges_with_attributes(self) -> None:
        """Test that parallel edges with different attributes are preserved."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2, weight=1.0, key=0)
        G.add_undirected_edge(1, 2, weight=2.0, key=1)
        G.add_undirected_edge(1, 2, weight=3.0, key=2)
        
        dm = orient_mixed_graph_from_root(G, 1)
        
        # All parallel edges should be oriented
        assert dm.number_of_edges() == 3
        edges_1_2 = list(dm._graph[1][2].values())
        weights = {e['weight'] for e in edges_1_2}
        assert weights == {1.0, 2.0, 3.0}
    
    def test_empty_graph(self) -> None:
        """Test orienting empty graph."""
        G = MixedMultiGraph()
        
        with pytest.raises(ValueError, match="not found"):
            orient_mixed_graph_from_root(G, 1)
    
    def test_star_graph(self) -> None:
        """Test orienting star graph (one central node connected to many leaves)."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(1, 3)
        G.add_undirected_edge(1, 4)
        G.add_undirected_edge(1, 5)
        G.add_undirected_edge(1, 6)
        
        dm = orient_mixed_graph_from_root(G, 1)
        
        # All edges should point away from center
        assert dm.number_of_edges() == 5
        assert dm.has_edge(1, 2)
        assert dm.has_edge(1, 3)
        assert dm.has_edge(1, 4)
        assert dm.has_edge(1, 5)
        assert dm.has_edge(1, 6)
    
    def test_path_graph(self) -> None:
        """Test orienting a path graph."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_undirected_edge(3, 4)
        G.add_undirected_edge(4, 5)
        
        dm = orient_mixed_graph_from_root(G, 1)
        
        # Should form a directed path
        assert dm.has_edge(1, 2)
        assert dm.has_edge(2, 3)
        assert dm.has_edge(3, 4)
        assert dm.has_edge(4, 5)
        assert dm.number_of_edges() == 4
    
    def test_path_graph_reverse_root(self) -> None:
        """Test orienting path graph from opposite end."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_undirected_edge(3, 4)
        G.add_undirected_edge(4, 5)
        
        dm = orient_mixed_graph_from_root(G, 5)
        
        # Should form a directed path AWAY from root 5
        assert dm.has_edge(5, 4)  # 5->4 (away from root)
        assert dm.has_edge(4, 3)  # 4->3 (away from root)
        assert dm.has_edge(3, 2)  # 3->2 (away from root)
        assert dm.has_edge(2, 1)  # 2->1 (away from root)
        assert dm.number_of_edges() == 4

