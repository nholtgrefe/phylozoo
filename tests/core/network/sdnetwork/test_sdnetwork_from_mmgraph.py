"""
Tests for sdnetwork_from_mmgraph function.
"""

import pytest

from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork, MixedPhyNetwork
from phylozoo.core.network.sdnetwork.io import sdnetwork_from_mmgraph
from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
from phylozoo.utils.validation import no_validation


class TestSdnetworkFromMmgraph:
    """Test conversion from MixedMultiGraph to SemiDirectedPhyNetwork/MixedPhyNetwork."""

    def test_empty_graph_semi_directed(self) -> None:
        """Convert empty graph to semi-directed network."""
        mmgraph = MixedMultiGraph()
        result = sdnetwork_from_mmgraph(mmgraph, network_type='semi-directed')
        assert isinstance(result, SemiDirectedPhyNetwork)
        assert result.number_of_nodes() == 0

    def test_empty_graph_mixed(self) -> None:
        """Convert empty graph to mixed network."""
        mmgraph = MixedMultiGraph()
        result = sdnetwork_from_mmgraph(mmgraph, network_type='mixed')
        assert isinstance(result, MixedPhyNetwork)
        assert result.number_of_nodes() == 0

    def test_default_network_type(self) -> None:
        """Test default network type is semi-directed."""
        mmgraph = MixedMultiGraph()
        result = sdnetwork_from_mmgraph(mmgraph)
        assert isinstance(result, SemiDirectedPhyNetwork)

    def test_single_node_no_attributes(self) -> None:
        """Convert graph with single node without attributes."""
        mmgraph = MixedMultiGraph()
        mmgraph.add_node(1)
        result = sdnetwork_from_mmgraph(mmgraph, network_type='semi-directed')
        assert isinstance(result, SemiDirectedPhyNetwork)
        assert result.number_of_nodes() == 1

    def test_single_node_with_label(self) -> None:
        """Convert graph with single node with label."""
        mmgraph = MixedMultiGraph()
        mmgraph.add_node(1, label='A')
        result = sdnetwork_from_mmgraph(mmgraph, network_type='semi-directed')
        assert isinstance(result, SemiDirectedPhyNetwork)
        assert result.get_label(1) == 'A'

    def test_directed_edges_only(self) -> None:
        """Convert graph with only directed edges."""
        mmgraph = MixedMultiGraph(directed_edges=[(1, 2), (2, 3)])
        mmgraph.add_node(3, label='A')
        with no_validation():
            result = sdnetwork_from_mmgraph(mmgraph, network_type='mixed')
        assert isinstance(result, MixedPhyNetwork)
        assert result.get_label(3) == 'A'

    def test_undirected_edges_only(self) -> None:
        """Convert graph with only undirected edges."""
        mmgraph = MixedMultiGraph(undirected_edges=[(1, 2), (2, 3)])
        mmgraph.add_node(3, label='A')
        with no_validation():
            result = sdnetwork_from_mmgraph(mmgraph, network_type='mixed')
        assert isinstance(result, MixedPhyNetwork)
        assert result.get_label(3) == 'A'

    def test_mixed_edges(self) -> None:
        """Convert graph with both directed and undirected edges."""
        mmgraph = MixedMultiGraph(
            directed_edges=[(1, 2)],
            undirected_edges=[(2, 3)]
        )
        mmgraph.add_node(3, label='A')
        with no_validation():
            result = sdnetwork_from_mmgraph(mmgraph, network_type='mixed')
        assert isinstance(result, MixedPhyNetwork)
        assert result.get_label(3) == 'A'

    def test_node_order_preserved(self) -> None:
        """Test that node ordering is preserved."""
        mmgraph = MixedMultiGraph(undirected_edges=[(1, 2), (2, 3), (2, 4)])
        mmgraph._undirected.nodes[3]['label'] = 'A'
        mmgraph._undirected.nodes[4]['label'] = 'B'
        result = sdnetwork_from_mmgraph(mmgraph, network_type='semi-directed')
        assert result.get_label(3) == 'A'
        assert result.get_label(4) == 'B'
