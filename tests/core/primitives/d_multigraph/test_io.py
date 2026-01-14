"""
Tests for directed multi-graph I/O functions.

This test suite covers the IOMixin integration for DirectedMultiGraph, including
save, load, convert, and format conversion methods for DOT and edge-list formats.
"""

import os
import tempfile

import pytest

from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph


class TestDirectedMultiGraphDOTIO:
    """Test IOMixin methods for DirectedMultiGraph with DOT format."""
    
    def test_to_string(self) -> None:
        """Test converting graph to DOT string."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0)
        G.add_edge(2, 3, weight=2.0)
        
        dot_str = G.to_string()
        
        assert 'digraph' in dot_str
        assert '1 -> 2' in dot_str
        assert '2 -> 3' in dot_str
    
    def test_to_string_with_attributes(self) -> None:
        """Test DOT string with node and edge attributes."""
        G = DirectedMultiGraph(attributes={'name': 'test'})
        G.add_node(1, label='Node1', shape='box')
        G.add_node(2, label='Node2')
        G.add_edge(1, 2, weight=1.0, color='blue')
        
        dot_str = G.to_string()
        
        assert 'digraph' in dot_str
        assert 'name=test' in dot_str
        assert 'label=Node1' in dot_str or 'label="Node1"' in dot_str
        assert 'weight=1.0' in dot_str
    
    def test_to_string_parallel_edges(self) -> None:
        """Test DOT string with parallel edges."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0)
        G.add_edge(1, 2, weight=2.0)  # Parallel edge
        
        dot_str = G.to_string()
        
        assert 'digraph' in dot_str
        assert dot_str.count('1 -> 2') == 2  # Two edges
    
    def test_save_basic(self) -> None:
        """Test basic DOT file saving."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.dot')
            G.save(file_path, overwrite=True)
            
            assert os.path.exists(file_path)
            
            # Read and verify content
            with open(file_path, 'r') as f:
                content = f.read()
            
            assert 'digraph' in content
            assert '1 -> 2' in content
    
    def test_save_gv_extension(self) -> None:
        """Test that .gv extension is accepted."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.gv')
            G.save(file_path, overwrite=True)
            assert os.path.exists(file_path)
    
    def test_from_string(self) -> None:
        """Test parsing DOT string to graph."""
        dot_str = '''digraph {
    1 [label="Node1"];
    2 [label="Node2"];
    1 -> 2 [weight=1.0];
    2 -> 3 [weight=2.0];
}'''
        
        G = DirectedMultiGraph.from_string(dot_str, format='dot')
        
        assert G.number_of_nodes() == 3
        assert G.number_of_edges() == 2
        assert 1 in G
        assert 2 in G
        assert 3 in G
    
    def test_from_string_with_attributes(self) -> None:
        """Test parsing DOT string with attributes."""
        dot_str = '''digraph {
    name=test_graph;
    1 [label="Node1", shape=box, color=red];
    2 [label="Node2"];
    1 -> 2 [weight=1.0, label=edge1, color=blue];
}'''
        
        G = DirectedMultiGraph.from_string(dot_str, format='dot')
        
        assert G.number_of_nodes() == 2
        assert G.number_of_edges() == 1
        assert G._graph.graph.get('name') == 'test_graph'
        
        # Check node attributes
        node1_data = G._graph.nodes[1]
        assert node1_data.get('label') == 'Node1'
        assert node1_data.get('shape') == 'box'
        
        # Check edge attributes
        for u, v, key, data in G.edges_iter(keys=True, data=True):
            assert data.get('weight') == 1.0
            assert data.get('label') == 'edge1'
    
    def test_from_string_parallel_edges(self) -> None:
        """Test parsing DOT string with parallel edges."""
        dot_str = '''digraph {
    1 -> 2 [weight=1.0, key=0];
    1 -> 2 [weight=2.0, key=1];
}'''
        
        G = DirectedMultiGraph.from_string(dot_str, format='dot')
        
        assert G.number_of_nodes() == 2
        assert G.number_of_edges() == 2
        
        # Check both edges exist
        edges = list(G.edges_iter(keys=True, data=True))
        assert len(edges) == 2
        weights = {data.get('weight') for _, _, _, data in edges}
        assert weights == {1.0, 2.0}
    
    def test_load(self) -> None:
        """Test loading graph from DOT file."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0)
        G.add_edge(2, 3, weight=2.0)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.dot')
            G.save(file_path, overwrite=True)
            
            # Load it back
            loaded_G = DirectedMultiGraph.load(file_path)
            
            assert loaded_G.number_of_nodes() == G.number_of_nodes()
            assert loaded_G.number_of_edges() == G.number_of_edges()
    
    def test_round_trip(self) -> None:
        """Test round-trip conversion: save and load."""
        G = DirectedMultiGraph(attributes={'name': 'test'})
        G.add_node(1, label='Node1', shape='box')
        G.add_node(2, label='Node2')
        G.add_edge(1, 2, weight=1.0, color='blue')
        G.add_edge(2, 3, weight=2.0)
        G.add_edge(1, 2, weight=3.0)  # Parallel edge
        
        dot_str = G.to_string()
        G2 = DirectedMultiGraph.from_string(dot_str)
        
        assert G2.number_of_nodes() == G.number_of_nodes()
        assert G2.number_of_edges() == G.number_of_edges()
    
    def test_empty_graph(self) -> None:
        """Test I/O with empty graph."""
        G = DirectedMultiGraph()
        
        dot_str = G.to_string()
        assert 'digraph' in dot_str
        
        G2 = DirectedMultiGraph.from_string(dot_str)
        assert G2.number_of_nodes() == 0
        assert G2.number_of_edges() == 0
    
    def test_convert(self) -> None:
        """Test format conversion between files."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = os.path.join(tmpdir, 'input.dot')
            output_file = os.path.join(tmpdir, 'output.dot')
            
            G.save(input_file, overwrite=True)
            DirectedMultiGraph.convert(input_file, output_file, overwrite=True)
            
            assert os.path.exists(output_file)
            loaded = DirectedMultiGraph.load(output_file)
            assert loaded.number_of_edges() == G.number_of_edges()
    
    def test_from_string_malformed_raises_error(self) -> None:
        """Test that malformed DOT string raises PhyloZooParseError."""
        from phylozoo.utils.exceptions import PhyloZooParseError
        with pytest.raises(PhyloZooParseError, match="Could not find"):
            DirectedMultiGraph.from_string("invalid dot", format='dot')
    
    def test_string_node_ids(self) -> None:
        """Test I/O with string node IDs."""
        G = DirectedMultiGraph()
        G.add_edge('A', 'B', weight=1.0)
        G.add_edge('B', 'C', weight=2.0)
        
        dot_str = G.to_string()
        G2 = DirectedMultiGraph.from_string(dot_str)
        
        assert G2.number_of_nodes() == G.number_of_nodes()
        assert G2.number_of_edges() == G.number_of_edges()
        assert 'A' in G2
        assert 'B' in G2
        assert 'C' in G2


class TestDirectedMultiGraphEdgeListIO:
    """Test IOMixin methods for DirectedMultiGraph with edge-list format."""
    
    def test_to_string(self) -> None:
        """Test converting graph to edge-list string."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0)
        G.add_edge(2, 3, weight=2.0)
        
        el_str = G.to_string(format='edgelist')
        
        assert '1 2' in el_str
        assert '2 3' in el_str
    
    def test_to_string_with_attributes(self) -> None:
        """Test edge-list string with edge attributes."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0, color='blue')
        G.add_edge(2, 3, weight=2.0)
        
        el_str = G.to_string(format='edgelist')
        
        assert '1 2' in el_str
        assert 'weight=1.0' in el_str
        assert 'color=blue' in el_str
    
    def test_to_string_parallel_edges(self) -> None:
        """Test edge-list string with parallel edges."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0)
        G.add_edge(1, 2, weight=2.0)  # Parallel edge
        
        el_str = G.to_string(format='edgelist')
        
        assert '1 2 0' in el_str or '1 2' in el_str
        assert '1 2 1' in el_str or el_str.count('1 2') == 2
    
    def test_save_basic(self) -> None:
        """Test basic edge-list file saving."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2)
        G.add_edge(2, 3)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.el')
            G.save(file_path, overwrite=True)
            
            assert os.path.exists(file_path)
            
            # Read and verify content
            with open(file_path, 'r') as f:
                content = f.read()
            
            assert '1 2' in content
            assert '2 3' in content
    
    def test_from_string(self) -> None:
        """Test parsing edge-list string to graph."""
        el_str = '''1 2
2 3
3 4'''
        
        G = DirectedMultiGraph.from_string(el_str, format='edgelist')
        
        assert G.number_of_nodes() == 4
        assert G.number_of_edges() == 3
        assert G.has_edge(1, 2)
        assert G.has_edge(2, 3)
        assert G.has_edge(3, 4)
    
    def test_from_string_with_attributes(self) -> None:
        """Test parsing edge-list string with attributes."""
        el_str = '''1 2 weight=1.0 color=blue
2 3 weight=2.0
3 4 0 key1=value1'''
        
        G = DirectedMultiGraph.from_string(el_str, format='edgelist')
        
        assert G.number_of_nodes() == 4
        assert G.number_of_edges() == 3
        
        # Check edge attributes
        for u, v, key, data in G.edges_iter(keys=True, data=True):
            if u == 1 and v == 2:
                assert data.get('weight') == 1.0
                assert data.get('color') == 'blue'
    
    def test_from_string_with_keys(self) -> None:
        """Test parsing edge-list string with explicit keys."""
        el_str = '''1 2 0 weight=1.0
1 2 1 weight=2.0'''
        
        G = DirectedMultiGraph.from_string(el_str, format='edgelist')
        
        assert G.number_of_nodes() == 2
        assert G.number_of_edges() == 2
        
        # Check both edges exist with correct keys
        edges = list(G.edges_iter(keys=True, data=True))
        assert len(edges) == 2
        keys = {key for _, _, key, _ in edges}
        assert keys == {0, 1}
    
    def test_load(self) -> None:
        """Test loading graph from edge-list file."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0)
        G.add_edge(2, 3, weight=2.0)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.el')
            G.save(file_path, format='edgelist', overwrite=True)
            
            # Load it back
            loaded_G = DirectedMultiGraph.load(file_path)
            
            assert loaded_G.number_of_nodes() == G.number_of_nodes()
            assert loaded_G.number_of_edges() == G.number_of_edges()
    
    def test_round_trip(self) -> None:
        """Test round-trip conversion: save and load."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0, color='blue')
        G.add_edge(2, 3, weight=2.0)
        G.add_edge(1, 2, weight=3.0)  # Parallel edge
        
        el_str = G.to_string(format='edgelist')
        G2 = DirectedMultiGraph.from_string(el_str, format='edgelist')
        
        assert G2.number_of_nodes() == G.number_of_nodes()
        assert G2.number_of_edges() == G.number_of_edges()
    
    def test_empty_graph(self) -> None:
        """Test I/O with empty graph."""
        G = DirectedMultiGraph()
        
        el_str = G.to_string(format='edgelist')
        assert el_str.strip() == '' or el_str == '\n'
        
        G2 = DirectedMultiGraph.from_string(el_str, format='edgelist')
        assert G2.number_of_nodes() == 0
        assert G2.number_of_edges() == 0
    
    def test_from_string_invalid_line_raises_error(self) -> None:
        """Test that invalid edge line raises PhyloZooParseError."""
        from phylozoo.utils.exceptions import PhyloZooParseError
        with pytest.raises(PhyloZooParseError, match="Invalid edge line"):
            DirectedMultiGraph.from_string("1", format='edgelist')
    
    def test_string_node_ids(self) -> None:
        """Test I/O with string node IDs."""
        G = DirectedMultiGraph()
        G.add_edge('A', 'B', weight=1.0)
        G.add_edge('B', 'C', weight=2.0)
        
        el_str = G.to_string(format='edgelist')
        G2 = DirectedMultiGraph.from_string(el_str, format='edgelist')
        
        assert G2.number_of_nodes() == G.number_of_nodes()
        assert G2.number_of_edges() == G.number_of_edges()
        assert 'A' in G2
        assert 'B' in G2
        assert 'C' in G2


class TestFormatConversion:
    """Test format conversion between DOT and edge-list."""
    
    def test_convert_dot_to_edgelist(self) -> None:
        """Test converting DOT file to edge-list."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0)
        G.add_edge(2, 3, weight=2.0)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = os.path.join(tmpdir, 'input.dot')
            output_file = os.path.join(tmpdir, 'output.el')
            
            G.save(input_file, overwrite=True)
            DirectedMultiGraph.convert(input_file, output_file, overwrite=True)
            
            assert os.path.exists(output_file)
            loaded = DirectedMultiGraph.load(output_file)
            assert loaded.number_of_edges() == G.number_of_edges()
    
    def test_convert_edgelist_to_dot(self) -> None:
        """Test converting edge-list file to DOT."""
        G = DirectedMultiGraph()
        G.add_edge(1, 2, weight=1.0)
        G.add_edge(2, 3, weight=2.0)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = os.path.join(tmpdir, 'input.el')
            output_file = os.path.join(tmpdir, 'output.dot')
            
            G.save(input_file, format='edgelist', overwrite=True)
            DirectedMultiGraph.convert(input_file, output_file, overwrite=True)
            
            assert os.path.exists(output_file)
            loaded = DirectedMultiGraph.load(output_file)
            assert loaded.number_of_edges() == G.number_of_edges()
    
    def test_convert_string(self) -> None:
        """Test string format conversion."""
        dot_str = '''digraph {
    1 -> 2 [weight=1.0];
    2 -> 3 [weight=2.0];
}'''
        
        # Convert DOT to edge-list
        el_str = DirectedMultiGraph.convert_string(dot_str, 'dot', 'edgelist')
        assert '1 2' in el_str
        assert '2 3' in el_str
        
        # Convert back to DOT
        dot_str2 = DirectedMultiGraph.convert_string(el_str, 'edgelist', 'dot')
        assert 'digraph' in dot_str2
        assert '1 -> 2' in dot_str2

