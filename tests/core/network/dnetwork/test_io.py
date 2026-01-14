"""
Tests for DirectedPhyNetwork I/O functions.

This test suite covers the IOMixin integration for DirectedPhyNetwork, including
save, load, convert, and format conversion methods for eNewick format.
"""

import os
import tempfile

import pytest

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.network.dnetwork._enewick import ENewickParseError


class TestDirectedPhyNetworkENewickIO:
    """Test IOMixin methods for DirectedPhyNetwork with eNewick format."""
    
    def test_to_string(self) -> None:
        """Test converting network to eNewick string."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        enewick_str = net.to_string()
        
        assert ';' in enewick_str
        assert 'A' in enewick_str
        assert 'B' in enewick_str
    
    def test_to_string_with_branch_lengths(self) -> None:
        """Test eNewick string with branch lengths."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'branch_length': 0.5},
                {'u': 3, 'v': 2, 'branch_length': 0.3}
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        enewick_str = net.to_string()
        
        assert ':0.5' in enewick_str
        assert ':0.3' in enewick_str
    
    def test_from_string(self) -> None:
        """Test parsing eNewick string to network."""
        enewick_str = "((A,B),C);"
        
        net = DirectedPhyNetwork.from_string(enewick_str, format='enewick')
        
        assert net.number_of_nodes() >= 3
        assert 'A' in net.taxa or 'A' in [net.get_label(n) for n in net.nodes()]
        assert 'B' in net.taxa or 'B' in [net.get_label(n) for n in net.nodes()]
        assert 'C' in net.taxa or 'C' in [net.get_label(n) for n in net.nodes()]
    
    def test_from_string_with_branch_lengths(self) -> None:
        """Test parsing eNewick string with branch lengths."""
        enewick_str = "((A:0.5,B:0.3):0.1,C:0.2);"
        
        net = DirectedPhyNetwork.from_string(enewick_str, format='enewick')
        
        assert net.number_of_nodes() >= 3
        # Check that branch lengths are preserved
        root = net.root_node
        has_branch_lengths = False
        for u, v, k, d in net._graph.edges(keys=True, data=True):
            if 'branch_length' in d:
                has_branch_lengths = True
                break
        assert has_branch_lengths
    
    def test_save_basic(self) -> None:
        """Test basic eNewick file saving."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.enewick')
            net.save(file_path, overwrite=True)
            
            assert os.path.exists(file_path)
            
            # Read and verify content
            with open(file_path, 'r') as f:
                content = f.read()
            
            assert ';' in content
            assert 'A' in content or 'B' in content
    
    def test_save_enw_extension(self) -> None:
        """Test that .enw extension is accepted."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.enw')
            net.save(file_path, overwrite=True)
            assert os.path.exists(file_path)
    
    def test_load_basic(self) -> None:
        """Test loading network from eNewick file."""
        enewick_content = "((A,B),C);"
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.enewick')
            with open(file_path, 'w') as f:
                f.write(enewick_content)
            
            net = DirectedPhyNetwork.load(file_path)
            
            assert net.number_of_nodes() >= 3
    
    def test_round_trip(self) -> None:
        """Test round-trip conversion (to_string -> from_string)."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'branch_length': 0.5},
                {'u': 3, 'v': 2, 'branch_length': 0.3}
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        enewick_str = net.to_string()
        net2 = DirectedPhyNetwork.from_string(enewick_str, format='enewick')
        
        assert net2.number_of_nodes() == net.number_of_nodes()
        assert net2.number_of_edges() == net.number_of_edges()
        assert net2.taxa == net.taxa
    
    def test_round_trip_file(self) -> None:
        """Test round-trip conversion via file."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.enewick')
            net.save(file_path, overwrite=True)
            
            net2 = DirectedPhyNetwork.load(file_path)
            
            assert net2.number_of_nodes() == net.number_of_nodes()
            assert net2.number_of_edges() == net.number_of_edges()
    
    def test_round_trip_valid_hybrid(self) -> None:
        """Test round-trip with valid hybrid node structure."""
        # Create a valid hybrid network (hybrid has out-degree 1)
        net = DirectedPhyNetwork(
            edges=[
                (7, 5), (7, 6),  # Root to tree nodes
                (5, 4), (5, 8),  # Tree node 5 -> hybrid 4 and leaf 8
                (6, 4), (6, 9),  # Tree node 6 -> hybrid 4 and leaf 9
                (4, 10),  # Hybrid 4 -> tree node 10 (out-degree 1)
                (10, 1), (10, 2)  # Tree node 10 -> leaves
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (8, {'label': 'C'}), (9, {'label': 'D'})]
        )
        
        enewick_str = net.to_string()
        net2 = DirectedPhyNetwork.from_string(enewick_str, format='enewick')
        
        assert net2.number_of_nodes() == net.number_of_nodes()
        assert net2.number_of_edges() == net.number_of_edges()
        assert len(net2.hybrid_nodes) == len(net.hybrid_nodes)
        assert net2.taxa == net.taxa
    
    def test_invalid_format(self) -> None:
        """Test error handling for invalid format."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        from phylozoo.utils.exceptions import PhyloZooFormatError
        with pytest.raises(PhyloZooFormatError, match="Format 'invalid' not supported"):
            net.to_string(format='invalid')
    
    def test_malformed_string(self) -> None:
        """Test error handling for malformed eNewick string."""
        malformed = "not a valid enewick"
        
        with pytest.raises(ENewickParseError):
            DirectedPhyNetwork.from_string(malformed, format='enewick')
    
    def test_missing_semicolon(self) -> None:
        """Test error handling for missing semicolon."""
        malformed = "((A,B),C)"  # Missing semicolon
        
        with pytest.raises(ENewickParseError, match="must end with ';'"):
            DirectedPhyNetwork.from_string(malformed, format='enewick')
    
    def test_convert(self) -> None:
        """Test format conversion between files."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, 'input.enewick')
            output_path = os.path.join(tmpdir, 'output.enewick')
            
            net.save(input_path, overwrite=True)
            DirectedPhyNetwork.convert(input_path, output_path, overwrite=True)
            
            net2 = DirectedPhyNetwork.load(output_path)
            assert net2.number_of_nodes() == net.number_of_nodes()
            assert net2.number_of_edges() == net.number_of_edges()
    
    def test_convert_string(self) -> None:
        """Test format conversion between strings."""
        enewick_str = "((A,B),C);"
        
        # Convert to same format (should be identity)
        result = DirectedPhyNetwork.convert_string(
            enewick_str, 'enewick', 'enewick'
        )
        
        net1 = DirectedPhyNetwork.from_string(enewick_str, format='enewick')
        net2 = DirectedPhyNetwork.from_string(result, format='enewick')
        
        assert net2.number_of_nodes() == net1.number_of_nodes()
        assert net2.number_of_edges() == net1.number_of_edges()
    
    def test_auto_detect_format(self) -> None:
        """Test automatic format detection from file extension."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.enewick')
            net.save(file_path, overwrite=True)
            
            # Load without specifying format (should auto-detect)
            net2 = DirectedPhyNetwork.load(file_path)
            assert net2.number_of_edges() == 2
    
    def test_explicit_format_override(self) -> None:
        """Test explicit format override even with extension."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, 'test.txt')
            net.save(file_path, format='enewick', overwrite=True)
            
            # Load with explicit format
            net2 = DirectedPhyNetwork.load(file_path, format='enewick')
            assert net2.number_of_edges() == 2
    
    def test_quoted_labels(self) -> None:
        """Test parsing eNewick with quoted labels."""
        enewick_str = "('Species A','Species B');"
        
        net = DirectedPhyNetwork.from_string(enewick_str, format='enewick')
        
        assert net.number_of_nodes() >= 2
        # Check that labels are preserved
        all_nodes = list(net._graph.nodes())
        labels = {net.get_label(n) for n in all_nodes if net.get_label(n)}
        assert 'Species A' in labels or 'Species B' in labels
    
    def test_internal_node_labels(self) -> None:
        """Test parsing eNewick with internal node labels."""
        enewick_str = "((A,B)internal,C);"
        
        net = DirectedPhyNetwork.from_string(enewick_str, format='enewick')
        
        # Check that internal node has label
        has_internal_label = False
        all_nodes = list(net._graph.nodes())
        for node in all_nodes:
            if node not in net.leaves and net.get_label(node) == 'internal':
                has_internal_label = True
                break
        assert has_internal_label
    
    def test_single_node_network(self) -> None:
        """Test single-node network I/O."""
        with pytest.warns(UserWarning, match="Single-node network"):
            net = DirectedPhyNetwork(nodes=[(1, {'label': 'A'})])
        
        enewick_str = net.to_string()
        assert enewick_str.endswith(';')
        
        net2 = DirectedPhyNetwork.from_string(enewick_str, format='enewick')
        assert net2.number_of_nodes() == 1


class TestDirectedPhyNetworkDOTIO:
    """Test IOMixin methods for DirectedPhyNetwork with DOT format."""
    
    def test_to_string_basic(self) -> None:
        """Test converting network to DOT string."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        dot_str = net.to_string(format='dot')
        
        assert 'digraph' in dot_str
        assert 'A' in dot_str
        assert 'B' in dot_str
        assert '->' in dot_str
    
    def test_to_string_uses_labels_as_node_names(self) -> None:
        """Test that labels are used as node names in DOT."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'Species_A'}), (2, {'label': 'Species_B'})]
        )
        
        dot_str = net.to_string(format='dot')
        
        # Labels should be used as node names
        assert 'Species_A' in dot_str
        assert 'Species_B' in dot_str
        # Node IDs should not appear for labeled nodes
        assert '1 ->' not in dot_str
        assert '2 ->' not in dot_str
        # But unlabeled nodes (like root) should use their ID
        assert '3' in dot_str or '3 ->' in dot_str
    
    def test_to_string_with_edge_attributes(self) -> None:
        """Test DOT string with edge attributes."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'branch_length': 0.5, 'bootstrap': 0.95},
                {'u': 3, 'v': 2, 'branch_length': 0.3}
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        dot_str = net.to_string(format='dot')
        
        assert 'branch_length=0.5' in dot_str
        assert 'bootstrap=0.95' in dot_str
        assert 'branch_length=0.3' in dot_str
    
    def test_to_string_with_internal_node_label(self) -> None:
        """Test DOT string with internal node label."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'root'})]
        )
        
        dot_str = net.to_string(format='dot')
        
        # Internal node label should be used as node name
        assert 'root' in dot_str
        assert 'root ->' in dot_str
        # Node ID should not appear
        assert '3 ->' not in dot_str
    
    def test_to_string_with_graph_attributes(self) -> None:
        """Test DOT string with graph attributes."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})],
            attributes={'name': 'TestNetwork', 'version': '1.0'}
        )
        
        dot_str = net.to_string(format='dot')
        
        assert 'name=' in dot_str or 'version=' in dot_str
    
    def test_from_string_basic(self) -> None:
        """Test parsing DOT string to network."""
        dot_str = '''digraph {
    A;
    B;
    root -> A;
    root -> B;
}'''
        
        net = DirectedPhyNetwork.from_string(dot_str, format='dot')
        
        assert net.number_of_nodes() == 3
        assert net.number_of_edges() == 2
        # Labels should be preserved
        assert 'A' in net.taxa or any(net.get_label(n) == 'A' for n in net.leaves)
        assert 'B' in net.taxa or any(net.get_label(n) == 'B' for n in net.leaves)
    
    def test_from_string_with_edge_attributes(self) -> None:
        """Test parsing DOT string with edge attributes."""
        dot_str = '''digraph {
    A;
    B;
    root -> A [branch_length=0.5, bootstrap=0.95];
    root -> B [branch_length=0.3];
}'''
        
        net = DirectedPhyNetwork.from_string(dot_str, format='dot')
        
        assert net.number_of_nodes() == 3
        # Check that edge attributes are preserved
        root = net.root_node
        leaves = list(net.leaves)
        assert len(leaves) == 2
        
        # Find edges and check attributes
        for leaf in leaves:
            bl = net.get_branch_length(root, leaf)
            if bl is not None:
                assert bl in [0.5, 0.3]
    
    def test_round_trip_basic(self) -> None:
        """Test round-trip conversion for basic network."""
        net1 = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        dot_str = net1.to_string(format='dot')
        net2 = DirectedPhyNetwork.from_string(dot_str, format='dot')
        
        assert net1.number_of_nodes() == net2.number_of_nodes()
        assert net1.number_of_edges() == net2.number_of_edges()
        assert net1.taxa == net2.taxa
    
    def test_round_trip_with_attributes(self) -> None:
        """Test round-trip conversion with edge attributes."""
        net1 = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'branch_length': 0.5, 'bootstrap': 0.95},
                {'u': 3, 'v': 2, 'branch_length': 0.3}
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        dot_str = net1.to_string(format='dot')
        net2 = DirectedPhyNetwork.from_string(dot_str, format='dot')
        
        assert net1.number_of_nodes() == net2.number_of_nodes()
        assert net1.number_of_edges() == net2.number_of_edges()
        # Check that edge attributes are preserved
        root1 = net1.root_node
        root2 = net2.root_node
        leaves1 = list(net1.leaves)
        leaves2 = list(net2.leaves)
        
        # Find matching leaves by label
        for leaf1 in leaves1:
            label1 = net1.get_label(leaf1)
            for leaf2 in leaves2:
                if net2.get_label(leaf2) == label1:
                    bl1 = net1.get_branch_length(root1, leaf1)
                    bl2 = net2.get_branch_length(root2, leaf2)
                    assert bl1 == bl2
    
    def test_round_trip_with_internal_label(self) -> None:
        """Test round-trip conversion with internal node label."""
        net1 = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'root'})]
        )
        
        dot_str = net1.to_string(format='dot')
        net2 = DirectedPhyNetwork.from_string(dot_str, format='dot')
        
        assert net1.number_of_nodes() == net2.number_of_nodes()
        assert net1.get_label(net1.root_node) == net2.get_label(net2.root_node)
    
    def test_save_and_load_dot(self) -> None:
        """Test saving and loading DOT format from file."""
        net1 = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.dot', delete=False) as f:
            temp_path = f.name
        
        try:
            net1.save(temp_path, format='dot', overwrite=True)
            
            net2 = DirectedPhyNetwork.load(temp_path, format='dot')
            
            assert net1.number_of_nodes() == net2.number_of_nodes()
            assert net1.number_of_edges() == net2.number_of_edges()
            assert net1.taxa == net2.taxa
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)
    
    def test_auto_detect_dot_format(self) -> None:
        """Test auto-detection of DOT format from file extension."""
        net1 = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.dot', delete=False) as f:
            temp_path = f.name
        
        try:
            net1.save(temp_path, overwrite=True)  # Should auto-detect .dot extension
            
            net2 = DirectedPhyNetwork.load(temp_path)  # Should auto-detect .dot extension
            
            assert net1.number_of_nodes() == net2.number_of_nodes()
            assert net1.taxa == net2.taxa
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)
    
    def test_auto_detect_gv_format(self) -> None:
        """Test auto-detection of DOT format from .gv extension."""
        net1 = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gv', delete=False) as f:
            temp_path = f.name
        
        try:
            net1.save(temp_path, overwrite=True)  # Should auto-detect .gv extension
            
            net2 = DirectedPhyNetwork.load(temp_path)  # Should auto-detect .gv extension
            
            assert net1.number_of_nodes() == net2.number_of_nodes()
            assert net1.taxa == net2.taxa
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)
    
    def test_nodes_without_labels(self) -> None:
        """Test DOT conversion for nodes without labels."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), 2]  # Node 2 has no label
        )
        
        dot_str = net.to_string(format='dot')
        
        # Node with label should use label
        assert 'A' in dot_str
        # Node without label should use ID
        assert '2' in dot_str or '2 ->' in dot_str or '-> 2' in dot_str
    
    def test_parallel_edges(self) -> None:
        """Test DOT conversion with parallel edges."""
        # Create a valid network with parallel edges (tree node with parallel edges to hybrid)
        net = DirectedPhyNetwork(
            edges=[
                (5, 3),  # Root to tree node
                (3, 4, 0), (3, 4, 1),  # Two parallel edges from tree node to hybrid
                (3, 2),  # Tree node also has another child
                (4, 1)  # Hybrid to leaf
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        dot_str = net.to_string(format='dot')
        
        # Should have key attribute for parallel edges or show multiple edges
        assert 'key=' in dot_str or dot_str.count('->') >= 4
    
    def test_convert_between_formats(self) -> None:
        """Test converting between DOT and eNewick formats."""
        net1 = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        
        # Convert to DOT
        dot_str = net1.to_string(format='dot')
        net2 = DirectedPhyNetwork.from_string(dot_str, format='dot')
        
        # Convert back to eNewick
        enewick_str = net2.to_string(format='enewick')
        net3 = DirectedPhyNetwork.from_string(enewick_str, format='enewick')
        
        assert net1.taxa == net3.taxa
        assert net1.number_of_nodes() == net3.number_of_nodes()
    
    def test_malformed_dot_string(self) -> None:
        """Test error handling for malformed DOT string."""
        from phylozoo.utils.exceptions import PhyloZooParseError
        malformed = "not a valid digraph"
        
        with pytest.raises(PhyloZooParseError, match="Could not find digraph"):
            DirectedPhyNetwork.from_string(malformed, format='dot')
    
    def test_empty_dot_string(self) -> None:
        """Test error handling for empty DOT string."""
        empty = "digraph {}"
        
        # Empty network may be valid or raise an error depending on validation
        # Let's check what actually happens
        try:
            net = DirectedPhyNetwork.from_string(empty, format='dot')
            # If it succeeds, network should be empty
            assert net.number_of_nodes() == 0
        except ValueError:
            # If it fails, that's also acceptable
            pass
    
    def test_dot_with_quoted_labels(self) -> None:
        """Test DOT conversion with quoted labels containing special characters."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'Species A'}), (2, {'label': 'Species-B'})]
        )
        
        dot_str = net.to_string(format='dot')
        
        # Quoted labels should be properly escaped
        assert 'Species A' in dot_str or '"Species A"' in dot_str
        assert 'Species-B' in dot_str

