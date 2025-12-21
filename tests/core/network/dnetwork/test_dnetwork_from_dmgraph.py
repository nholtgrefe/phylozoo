"""
Tests for dnetwork_from_dmgraph function.
"""

import pytest

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.io import dnetwork_from_dmgraph
from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
from phylozoo.utils.validation import no_validation


class TestDnetworkFromDmgraph:
    """Test conversion from DirectedMultiGraph to DirectedPhyNetwork."""

    def test_empty_graph(self) -> None:
        """Convert empty graph."""
        dmgraph = DirectedMultiGraph()
        result = dnetwork_from_dmgraph(dmgraph)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_nodes() == 0
        assert result.number_of_edges() == 0

    def test_single_node_no_attributes(self) -> None:
        """Convert graph with single node without attributes."""
        dmgraph = DirectedMultiGraph()
        dmgraph.add_node(1)
        result = dnetwork_from_dmgraph(dmgraph)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_nodes() == 1
        assert result.number_of_edges() == 0

    def test_single_node_with_label(self) -> None:
        """Convert graph with single node with label."""
        dmgraph = DirectedMultiGraph()
        dmgraph.add_node(1, label='A')
        result = dnetwork_from_dmgraph(dmgraph)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_nodes() == 1
        assert result.get_label(1) == 'A'

    def test_simple_tree(self) -> None:
        """Convert simple tree structure."""
        dmgraph = DirectedMultiGraph(edges=[(3, 1), (3, 2)])
        dmgraph._graph.nodes[1]['label'] = 'A'
        dmgraph._graph.nodes[2]['label'] = 'B'
        result = dnetwork_from_dmgraph(dmgraph)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_nodes() == 3
        assert result.number_of_edges() == 2
        assert result.get_label(1) == 'A'
        assert result.get_label(2) == 'B'
        assert result.root_node == 3

    def test_edges_with_attributes(self) -> None:
        """Convert graph with edges that have attributes."""
        dmgraph = DirectedMultiGraph(edges=[
            {'u': 1, 'v': 2, 'branch_length': 0.5, 'bootstrap': 0.9}
        ])
        result = dnetwork_from_dmgraph(dmgraph)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_edges() == 1
        # Check that edge attributes are preserved
        edges = list(result._graph.edges(data=True))
        assert len(edges) == 1
        u, v, data = edges[0]
        assert data.get('branch_length') == 0.5
        assert data.get('bootstrap') == 0.9

    def test_parallel_edges(self) -> None:
        """Convert graph with parallel edges."""
        dmgraph = DirectedMultiGraph(edges=[
            {'u': 1, 'v': 2, 'key': 0, 'branch_length': 0.5},
            {'u': 1, 'v': 2, 'key': 1, 'branch_length': 0.5}  # Same branch_length for parallel edges
        ])
        with no_validation():
            result = dnetwork_from_dmgraph(dmgraph)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_edges() == 2

    def test_network_validity(self) -> None:
        """Ensure resulting network is valid."""
        dmgraph = DirectedMultiGraph(edges=[(3, 1), (3, 2)])
        dmgraph._graph.nodes[1]['label'] = 'A'
        dmgraph._graph.nodes[2]['label'] = 'B'
        result = dnetwork_from_dmgraph(dmgraph)
        # Should not raise an exception
        result.validate()

    def test_node_order_preserved(self) -> None:
        """Test that node ordering is preserved."""
        dmgraph = DirectedMultiGraph(edges=[(1, 2), (2, 3), (2, 4)])
        dmgraph._graph.nodes[3]['label'] = 'A'
        dmgraph._graph.nodes[4]['label'] = 'B'
        result = dnetwork_from_dmgraph(dmgraph)
        assert result.get_label(3) == 'A'
        assert result.get_label(4) == 'B'
