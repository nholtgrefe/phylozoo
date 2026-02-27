"""
Tests for orient_away_from_vertex function.

This test suite provides comprehensive coverage of the orient_away_from_vertex
function, including basic functionality, edge cases, error handling, and complex structures.
"""

import pytest
from typing import Dict, List, Set, Tuple

from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
from phylozoo.core.primitives.m_multigraph.transformations import orient_away_from_vertex
from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph


class TestOrientMixedGraphFromRoot:
    """Test cases for orient_away_from_vertex function."""
    
    def test_simple_undirected_tree(self) -> None:
        """Test orienting a simple undirected tree."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(1, 3)
        G.add_undirected_edge(2, 4)
        G.add_undirected_edge(2, 5)
        
        dm = orient_away_from_vertex(G, 1)
        
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
        
        dm = orient_away_from_vertex(G, 1)
        
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
        
        dm = orient_away_from_vertex(G, 1)
        
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
        
        dm = orient_away_from_vertex(G, 1)
        
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
        
        dm = orient_away_from_vertex(G, 1)
        
        # All edges should be oriented
        assert dm.has_edge(1, 2)
        assert dm.has_edge(1, 3)
        assert dm.has_edge(2, 4)
        assert dm.has_edge(3, 4)
        assert dm.has_edge(4, 5)
        assert dm.has_edge(4, 6)
        assert dm.number_of_edges() == G.number_of_edges()

    def test_directed_edge_into_vertex_reached_via_undirected(self) -> None:
        """
        Test orienting when a directed edge points into a vertex already on the BFS path.

        Root 6 is connected by undirected edges to 4 and 1. Vertex 1 is reached via 6-1.
        The directed edge 3->1 points into 1; this is valid (hybrid node) and should not raise.
        """
        G = MixedMultiGraph(
            directed_edges=[(4, 5, 0), (4, 5, 1), (2, 3, 0), (2, 3, 1), (3, 1, 0)],
            undirected_edges=[(4, 6, 0), (5, 2, 0), (1, 6, 0)],
        )
        dm = orient_away_from_vertex(G, 6)
        assert dm.number_of_edges() == G.number_of_edges()
        assert dm.number_of_nodes() == G.number_of_nodes()
        # Directed edge 3->1 preserved (hybrid at 1); undirected oriented away from root 6
        assert dm.has_edge(3, 1)
        assert dm.has_edge(6, 4)
        assert dm.has_edge(6, 1)
    
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
        
        dm = orient_away_from_vertex(G, 1)
        
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
            orient_away_from_vertex(G, 999)
    
    def test_unreachable_undirected_edge_raises_error(self) -> None:
        """Test that unreachable undirected edge raises error."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(3, 4)  # Disconnected component
        
        with pytest.raises(ValueError, match="not reachable from root"):
            orient_away_from_vertex(G, 1)
    
    def test_single_node_graph(self) -> None:
        """Test orienting graph with single node."""
        G = MixedMultiGraph()
        G.add_node(1)
        
        dm = orient_away_from_vertex(G, 1)
        
        assert dm.number_of_nodes() == 1
        assert dm.number_of_edges() == 0
    
    def test_single_edge_graph(self) -> None:
        """Test orienting graph with single edge."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        
        dm = orient_away_from_vertex(G, 1)
        
        assert dm.has_edge(1, 2)
        assert dm.number_of_edges() == 1
        assert dm.number_of_nodes() == 2
    
    def test_preserves_edge_attributes(self) -> None:
        """Test that edge attributes are preserved during orientation."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2, weight=1.5, label='test')
        G.add_undirected_edge(2, 3, weight=2.0)
        G.add_directed_edge(3, 4, weight=3.0)
        
        dm = orient_away_from_vertex(G, 1)
        
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
        
        dm = orient_away_from_vertex(G, 1)
        
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
        
        dm = orient_away_from_vertex(G, 1)
        
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
        
        dm = orient_away_from_vertex(G, 1)
        
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
        
        dm = orient_away_from_vertex(G, 1)
        
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
        
        dm1 = orient_away_from_vertex(G, 1)
        dm2 = orient_away_from_vertex(G, 3)
        
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
        
        dm = orient_away_from_vertex(G, 1)
        
        # All parallel edges should be oriented
        assert dm.number_of_edges() == 3
        edges_1_2 = list(dm._graph[1][2].values())
        weights = {e['weight'] for e in edges_1_2}
        assert weights == {1.0, 2.0, 3.0}
    
    def test_empty_graph(self) -> None:
        """Test orienting empty graph."""
        G = MixedMultiGraph()
        
        with pytest.raises(ValueError, match="not found"):
            orient_away_from_vertex(G, 1)
    
    def test_star_graph(self) -> None:
        """Test orienting star graph (one central node connected to many leaves)."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(1, 3)
        G.add_undirected_edge(1, 4)
        G.add_undirected_edge(1, 5)
        G.add_undirected_edge(1, 6)
        
        dm = orient_away_from_vertex(G, 1)
        
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
        
        dm = orient_away_from_vertex(G, 1)
        
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
        
        dm = orient_away_from_vertex(G, 5)
        
        # Should form a directed path AWAY from root 5
        assert dm.has_edge(5, 4)  # 5->4 (away from root)
        assert dm.has_edge(4, 3)  # 4->3 (away from root)
        assert dm.has_edge(3, 2)  # 3->2 (away from root)
        assert dm.has_edge(2, 1)  # 2->1 (away from root)
        assert dm.number_of_edges() == 4
    
    def test_no_duplicate_edges_bug_fix(self) -> None:
        """
        Test that duplicate edges are not created (bug fix verification).
        
        This test specifically verifies that the fix for duplicate edge processing
        works correctly. Previously, undirected edges could be processed twice
        (once from each endpoint), creating duplicate edges in the result.
        """
        G = MixedMultiGraph()
        # Create a structure where edges could be processed from multiple directions
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_undirected_edge(3, 4)
        G.add_undirected_edge(2, 5)  # Branch from node 2
        G.add_undirected_edge(3, 6)  # Branch from node 3
        
        dm = orient_away_from_vertex(G, 1)
        
        # Verify exact edge count matches (no duplicates)
        assert dm.number_of_edges() == G.number_of_edges()
        assert dm.number_of_edges() == 5
        
        # Verify each edge appears exactly once
        edges_list = list(dm.edges())
        edges_set = set(edges_list)
        assert len(edges_list) == len(edges_set), "Duplicate edges found!"
        
        # Verify all expected edges exist
        assert dm.has_edge(1, 2)
        assert dm.has_edge(2, 3)
        assert dm.has_edge(3, 4)
        assert dm.has_edge(2, 5)
        assert dm.has_edge(3, 6)
    
    def test_no_duplicate_edges_with_parallel_edges(self) -> None:
        """
        Test that parallel edges are not duplicated (bug fix verification).
        
        This test verifies that parallel undirected edges are correctly handled
        and not duplicated during orientation.
        """
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2, key=0, weight=1.0)
        G.add_undirected_edge(1, 2, key=1, weight=2.0)  # Parallel edge
        G.add_undirected_edge(2, 3, key=0, weight=3.0)
        
        dm = orient_away_from_vertex(G, 1)
        
        # Verify exact edge count (should have both parallel edges)
        assert dm.number_of_edges() == G.number_of_edges()
        assert dm.number_of_edges() == 3
        
        # Verify both parallel edges exist
        assert dm.has_edge(1, 2, key=0)
        assert dm.has_edge(1, 2, key=1)
        assert dm.has_edge(2, 3, key=0)
        
        # Verify no duplicates
        edges_list = list(dm.edges(keys=True))
        edges_set = set(edges_list)
        assert len(edges_list) == len(edges_set), "Duplicate edges found!"
    
    def test_node_attributes_preserved_comprehensive(self) -> None:
        """
        Test that node attributes are correctly copied (bug fix verification).
        
        This test verifies the fix for node attribute copying. Previously,
        attributes weren't being copied because the code tried to access
        a non-existent _node attribute.
        """
        G = MixedMultiGraph()
        # Add nodes with various attributes
        G.add_node(1, label='root', color='red', weight=10.0)
        G.add_node(2, label='internal', color='blue', weight=5.0)
        G.add_node(3, label='leaf1', color='green')
        G.add_node(4, label='leaf2', color='yellow', weight=1.0)
        
        # Add edges
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_undirected_edge(2, 4)
        
        dm = orient_away_from_vertex(G, 1)
        
        # Verify all nodes have their attributes
        node_data = dict(dm.nodes(data=True))
        
        assert node_data[1].get('label') == 'root'
        assert node_data[1].get('color') == 'red'
        assert node_data[1].get('weight') == 10.0
        
        assert node_data[2].get('label') == 'internal'
        assert node_data[2].get('color') == 'blue'
        assert node_data[2].get('weight') == 5.0
        
        assert node_data[3].get('label') == 'leaf1'
        assert node_data[3].get('color') == 'green'
        assert 'weight' not in node_data[3] or node_data[3].get('weight') is None
        
        assert node_data[4].get('label') == 'leaf2'
        assert node_data[4].get('color') == 'yellow'
        assert node_data[4].get('weight') == 1.0
    
    def test_node_attributes_from_both_graphs(self) -> None:
        """
        Test that node attributes from both _undirected and _directed graphs are preserved.
        
        This test verifies that when a node exists in both underlying graphs,
        attributes from both are correctly merged.
        """
        G = MixedMultiGraph()
        # Add node with attributes via add_node
        G.add_node(1, label='node1', attr1='value1')
        # Node will be in both _undirected and _directed after adding edges
        G.add_undirected_edge(1, 2)
        G.add_directed_edge(1, 3)
        
        # Add attributes to node via the underlying graphs (simulating real usage)
        G._undirected.nodes[1]['attr2'] = 'value2'
        G._directed.nodes[1]['attr3'] = 'value3'
        
        dm = orient_away_from_vertex(G, 1)
        
        node_data = dict(dm.nodes(data=True))
        
        # Verify all attributes are present
        assert node_data[1].get('label') == 'node1'
        assert node_data[1].get('attr1') == 'value1'
        assert node_data[1].get('attr2') == 'value2'
        assert node_data[1].get('attr3') == 'value3'
    
    def test_edge_count_always_matches(self) -> None:
        """
        Test that edge count always matches input graph (comprehensive bug fix verification).
        
        This test verifies that in various scenarios, the output graph always
        has exactly the same number of edges as the input graph, ensuring
        no duplicates and no missing edges.
        """
        test_cases = [
            # (description, graph_setup_function, valid_roots)
            ("Simple tree", lambda: self._create_simple_tree(), [1, 2, 3, 4, 5]),
            ("Star graph", lambda: self._create_star_graph(), [1, 2, 3, 4, 5]),
            ("Mixed directed/undirected", lambda: self._create_mixed_graph(), [1]),  # Only root 1 is valid
            ("Complex branching", lambda: self._create_complex_branching(), [1, 2, 3]),
        ]
        
        for desc, graph_func, valid_roots in test_cases:
            G = graph_func()
            original_edge_count = G.number_of_edges()
            
            # Test with valid roots only (some graphs have constraints on valid roots)
            for root in valid_roots:
                if root not in G.nodes():
                    continue
                    
                dm = orient_away_from_vertex(G, root)
                assert dm.number_of_edges() == original_edge_count, \
                    f"Edge count mismatch for {desc} with root {root}: " \
                    f"expected {original_edge_count}, got {dm.number_of_edges()}"
                
                # Also verify no duplicate edges
                edges_list = list(dm.edges(keys=True))
                edges_set = set(edges_list)
                assert len(edges_list) == len(edges_set), \
                    f"Duplicate edges found for {desc} with root {root}"
    
    def _create_simple_tree(self) -> MixedMultiGraph:
        """Helper to create a simple tree."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(1, 3)
        G.add_undirected_edge(2, 4)
        G.add_undirected_edge(2, 5)
        return G
    
    def _create_star_graph(self) -> MixedMultiGraph:
        """Helper to create a star graph."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(1, 3)
        G.add_undirected_edge(1, 4)
        G.add_undirected_edge(1, 5)
        return G
    
    def _create_mixed_graph(self) -> MixedMultiGraph:
        """Helper to create a mixed directed/undirected graph."""
        G = MixedMultiGraph()
        G.add_directed_edge(1, 2)
        G.add_undirected_edge(2, 3)
        G.add_undirected_edge(2, 4)
        G.add_directed_edge(3, 5)
        G.add_undirected_edge(4, 6)
        return G
    
    def _create_complex_branching(self) -> MixedMultiGraph:
        """Helper to create a complex branching structure."""
        G = MixedMultiGraph()
        G.add_undirected_edge(1, 2)
        G.add_undirected_edge(1, 3)
        G.add_undirected_edge(2, 4)
        G.add_undirected_edge(2, 5)
        G.add_undirected_edge(3, 6)
        G.add_undirected_edge(3, 7)
        G.add_undirected_edge(4, 8)
        G.add_undirected_edge(5, 9)
        return G

