"""
Tests for binary_resolution transformation.

This test suite covers binary resolution of non-binary networks, including
high out-degree nodes, high in-degree nodes (hybrids), gamma preservation,
branch length handling, and edge cases.
"""

import pytest

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.transformations import binary_resolution
from phylozoo.core.network.dnetwork.classifications import is_binary, has_parallel_edges


class TestBinaryResolution:
    """Test cases for binary_resolution function."""
    
    def test_already_binary_network(self) -> None:
        """Test that binary networks are returned as copies."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        result = binary_resolution(net)
        assert is_binary(result)
        assert result.taxa == net.taxa
        assert result is not net  # Should be a copy
    
    def test_high_outdegree_node(self) -> None:
        """Test resolution of node with out-degree > 2."""
        net = DirectedPhyNetwork(
            edges=[(5, 1), (5, 2), (5, 3), (5, 4)],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (4, {'label': 'D'})
            ]
        )
        result = binary_resolution(net)
        assert is_binary(result)
        assert result.taxa == net.taxa
    
    
    def test_with_branch_lengths(self) -> None:
        """Test that new edges get branch_length=0.0 when original has branch lengths."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 5, 'v': 1, 'branch_length': 0.5},
                {'u': 5, 'v': 2, 'branch_length': 0.3},
                {'u': 5, 'v': 3, 'branch_length': 0.2}
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'})
            ]
        )
        result = binary_resolution(net)
        assert is_binary(result)
        
        # Node 5 is removed, find the caterpillar node that connects to 1
        cat_nodes = [n for n in result._graph.nodes() if n not in [1, 2, 3]]
        for cat in cat_nodes:
            if result._graph.has_edge(cat, 1):
                # Original edge 5->1 had branch_length 0.5, should be preserved
                assert result.get_branch_length(cat, 1) == 0.5
        
        # Check that all edges have branch_length (either original or 0.0)
        for u, v, key, data in result._graph.edges(keys=True, data=True):
            assert 'branch_length' in data
    
    def test_attributes_removed_except_branch_length_gamma(self) -> None:
        """Test that all attributes except branch_length and gamma are removed."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 5, 'v': 1, 'branch_length': 0.5, 'bootstrap': 0.9, 'custom': 'value'},
                {'u': 5, 'v': 2, 'branch_length': 0.3, 'bootstrap': 0.8},
                {'u': 5, 'v': 3, 'branch_length': 0.2}
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'})
            ]
        )
        result = binary_resolution(net)
        
        # Check that bootstrap and custom attributes are removed
        for u, v, key, data in result._graph.edges(keys=True, data=True):
            assert 'bootstrap' not in data
            assert 'custom' not in data
            # branch_length should be preserved
            if 'branch_length' in data:
                assert isinstance(data['branch_length'], (int, float))

    def test_without_branch_lengths(self) -> None:
        """Test that new edges have no branch_length when original has none."""
        net = DirectedPhyNetwork(
            edges=[(5, 1), (5, 2), (5, 3)],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'})
            ]
        )
        result = binary_resolution(net)
        assert is_binary(result)
        
        # Check that edges have no branch_length
        for u, v, key, data in result._graph.edges(keys=True, data=True):
            assert 'branch_length' not in data or data.get('branch_length') is None
    
    def test_parallel_edges_raises_error(self) -> None:
        """Test that networks with parallel edges raise ValueError."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 5),  # Root
                (5, 4, 0), (5, 4, 1),  # Parallel edges
                (4, 1)
            ],
            nodes=[
                (1, {'label': 'A'}),
            ]
        )
        with pytest.raises(ValueError, match="parallel edges"):
            binary_resolution(net)
    
    def test_preserves_labels(self) -> None:
        """Test that node labels are preserved."""
        net = DirectedPhyNetwork(
            edges=[(5, 1), (5, 2), (5, 3)],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (5, {'label': 'Internal'})
            ]
        )
        result = binary_resolution(net)
        assert result.get_label(1) == 'A'
        assert result.get_label(2) == 'B'
        assert result.get_label(3) == 'C'
        # Node 5 is removed, so its label is gone
    
    def test_single_node_network(self) -> None:
        """Test that single-node networks are handled correctly."""
        net = DirectedPhyNetwork(
            nodes=[(1, {'label': 'A'})]
        )
        result = binary_resolution(net)
        assert is_binary(result)
        assert result.taxa == net.taxa
    
    def test_empty_network(self) -> None:
        """Test that empty networks are handled correctly."""
        net = DirectedPhyNetwork(edges=[], nodes=[])
        result = binary_resolution(net)
        assert is_binary(result)
        assert result.number_of_nodes() == 0
    
    def test_fixture_network(self) -> None:
        """Test binary resolution on a fixture network."""
        from tests.fixtures.directed_networks import LEVEL_2_DNETWORK_DIAMOND_HYBRID
        
        net = LEVEL_2_DNETWORK_DIAMOND_HYBRID
        result = binary_resolution(net)
        assert is_binary(result)
        assert result.taxa == net.taxa
        assert not has_parallel_edges(result)
        
        # Check that high-degree nodes are removed
        # Node 10 has out-degree 3, node 11 has out-degree 3, node 6 has in-degree 3
        assert 10 not in result._graph.nodes()
        assert 11 not in result._graph.nodes()
        assert 6 not in result._graph.nodes()
