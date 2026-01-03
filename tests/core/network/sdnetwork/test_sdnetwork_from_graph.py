"""
Tests for sdnetwork_from_graph function.
"""

import pytest
import networkx as nx

from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork, MixedPhyNetwork
from phylozoo.core.network.sdnetwork.conversions import sdnetwork_from_graph
from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
from phylozoo.utils.validation import no_validation


class TestSdnetworkFromGraph:
    """Test conversion from various graph types to SemiDirectedPhyNetwork/MixedPhyNetwork."""

    def test_empty_graph_semi_directed_from_mmgraph(self) -> None:
        """Convert empty MixedMultiGraph to semi-directed network."""
        mmgraph = MixedMultiGraph()
        result = sdnetwork_from_graph(mmgraph, network_type='semi-directed')
        assert isinstance(result, SemiDirectedPhyNetwork)
        assert result.number_of_nodes() == 0

    def test_empty_graph_semi_directed_from_nx_graph(self) -> None:
        """Convert empty NetworkX Graph to semi-directed network."""
        graph = nx.Graph()
        result = sdnetwork_from_graph(graph, network_type='semi-directed')
        assert isinstance(result, SemiDirectedPhyNetwork)
        assert result.number_of_nodes() == 0

    def test_empty_graph_mixed_from_mmgraph(self) -> None:
        """Convert empty MixedMultiGraph to mixed network."""
        mmgraph = MixedMultiGraph()
        result = sdnetwork_from_graph(mmgraph, network_type='mixed')
        assert isinstance(result, MixedPhyNetwork)
        assert result.number_of_nodes() == 0

    def test_default_network_type(self) -> None:
        """Test default network type is semi-directed."""
        mmgraph = MixedMultiGraph()
        result = sdnetwork_from_graph(mmgraph)
        assert isinstance(result, SemiDirectedPhyNetwork)

    def test_single_node_no_attributes_from_mmgraph(self) -> None:
        """Convert MixedMultiGraph with single node without attributes."""
        mmgraph = MixedMultiGraph()
        mmgraph.add_node(1)
        result = sdnetwork_from_graph(mmgraph, network_type='semi-directed')
        assert isinstance(result, SemiDirectedPhyNetwork)
        assert result.number_of_nodes() == 1

    def test_single_node_with_label_from_mmgraph(self) -> None:
        """Convert MixedMultiGraph with single node with label."""
        mmgraph = MixedMultiGraph()
        mmgraph.add_node(1, label='A')
        result = sdnetwork_from_graph(mmgraph, network_type='semi-directed')
        assert isinstance(result, SemiDirectedPhyNetwork)
        assert result.get_label(1) == 'A'

    def test_single_node_with_label_from_nx_graph(self) -> None:
        """Convert NetworkX Graph with single node with label."""
        graph = nx.Graph()
        graph.add_node(1)
        graph.nodes[1]['label'] = 'A'
        result = sdnetwork_from_graph(graph, network_type='semi-directed')
        assert isinstance(result, SemiDirectedPhyNetwork)
        assert result.get_label(1) == 'A'

    def test_directed_edges_only_from_mmgraph(self) -> None:
        """Convert MixedMultiGraph with only directed edges."""
        mmgraph = MixedMultiGraph(directed_edges=[(1, 2), (2, 3)])
        mmgraph.add_node(3, label='A')
        with no_validation():
            result = sdnetwork_from_graph(mmgraph, network_type='mixed')
        assert isinstance(result, MixedPhyNetwork)
        assert result.get_label(3) == 'A'

    def test_undirected_edges_only_from_mmgraph(self) -> None:
        """Convert MixedMultiGraph with only undirected edges."""
        mmgraph = MixedMultiGraph(undirected_edges=[(0, 1), (0, 2), (0, 3)])
        mmgraph.add_node(1, label='A')
        mmgraph.add_node(2, label='B')
        mmgraph.add_node(3, label='C')
        result = sdnetwork_from_graph(mmgraph, network_type='semi-directed')
        assert isinstance(result, SemiDirectedPhyNetwork)
        assert result.get_label(1) == 'A'

    def test_undirected_edges_only_from_nx_graph(self) -> None:
        """Convert NetworkX Graph with only undirected edges."""
        graph = nx.Graph()
        graph.add_edge(0, 1)
        graph.add_edge(0, 2)
        graph.add_edge(0, 3)
        graph.nodes[1]['label'] = 'A'
        graph.nodes[2]['label'] = 'B'
        graph.nodes[3]['label'] = 'C'
        result = sdnetwork_from_graph(graph, network_type='semi-directed')
        assert isinstance(result, SemiDirectedPhyNetwork)
        assert result.get_label(1) == 'A'

    def test_mixed_edges_from_mmgraph(self) -> None:
        """Convert MixedMultiGraph with both directed and undirected edges."""
        mmgraph = MixedMultiGraph(
            directed_edges=[(1, 2)],
            undirected_edges=[(2, 3)]
        )
        mmgraph.add_node(3, label='A')
        with no_validation():
            result = sdnetwork_from_graph(mmgraph, network_type='mixed')
        assert isinstance(result, MixedPhyNetwork)
        assert result.get_label(3) == 'A'

    def test_edges_with_attributes_from_nx_graph(self) -> None:
        """Convert NetworkX Graph with edges that have attributes."""
        graph = nx.Graph()
        graph.add_edge(0, 1, branch_length=0.5, bootstrap=0.9)
        graph.nodes[1]['label'] = 'A'
        result = sdnetwork_from_graph(graph, network_type='semi-directed')
        assert isinstance(result, SemiDirectedPhyNetwork)
        assert result.get_branch_length(0, 1) == 0.5
        assert result.get_bootstrap(0, 1) == 0.9

    def test_graph_attributes_preserved_from_mmgraph(self) -> None:
        """Test that graph attributes are preserved from MixedMultiGraph."""
        mmgraph = MixedMultiGraph(
            undirected_edges=[(0, 1)],
            attributes={'source': 'test', 'version': '1.0'}
        )
        mmgraph.add_node(1, label='A')
        result = sdnetwork_from_graph(mmgraph, network_type='semi-directed')
        assert result.get_network_attribute('source') == 'test'
        assert result.get_network_attribute('version') == '1.0'

    def test_graph_attributes_preserved_from_nx_graph(self) -> None:
        """Test that graph attributes are preserved from NetworkX Graph."""
        graph = nx.Graph()
        graph.add_edge(0, 1)
        graph.nodes[1]['label'] = 'A'
        graph.graph['source'] = 'test'
        graph.graph['version'] = '1.0'
        result = sdnetwork_from_graph(graph, network_type='semi-directed')
        assert result.get_network_attribute('source') == 'test'
        assert result.get_network_attribute('version') == '1.0'

    def test_parallel_edges_from_nx_multigraph(self) -> None:
        """Convert NetworkX MultiGraph with parallel edges."""
        graph = nx.MultiGraph()
        graph.add_edge(0, 1, branch_length=0.5)
        graph.add_edge(0, 1, branch_length=0.3)
        graph.nodes[1]['label'] = 'A'
        with no_validation():
            result = sdnetwork_from_graph(graph, network_type='semi-directed')
        assert isinstance(result, SemiDirectedPhyNetwork)
        assert result.number_of_edges() == 2

    def test_simple_tree_from_nx_multigraph(self) -> None:
        """Convert simple tree structure from NetworkX MultiGraph."""
        graph = nx.MultiGraph()
        graph.add_edge(0, 1)
        graph.add_edge(0, 2)
        graph.add_edge(0, 3)
        graph.nodes[1]['label'] = 'A'
        graph.nodes[2]['label'] = 'B'
        graph.nodes[3]['label'] = 'C'
        result = sdnetwork_from_graph(graph, network_type='semi-directed')
        assert isinstance(result, SemiDirectedPhyNetwork)
        assert result.number_of_nodes() == 4
        assert result.number_of_edges() == 3
        assert result.get_label(1) == 'A'
        assert result.get_label(2) == 'B'
        assert result.get_label(3) == 'C'

    def test_edges_with_attributes_from_nx_multigraph(self) -> None:
        """Convert NetworkX MultiGraph with edges that have attributes."""
        graph = nx.MultiGraph()
        graph.add_edge(0, 1, branch_length=0.5, bootstrap=0.9)
        graph.nodes[1]['label'] = 'A'
        result = sdnetwork_from_graph(graph, network_type='semi-directed')
        assert isinstance(result, SemiDirectedPhyNetwork)
        assert result.get_branch_length(0, 1) == 0.5
        assert result.get_bootstrap(0, 1) == 0.9

    def test_graph_attributes_preserved_from_nx_multigraph(self) -> None:
        """Test that graph attributes are preserved from NetworkX MultiGraph."""
        graph = nx.MultiGraph()
        graph.add_edge(0, 1)
        graph.nodes[1]['label'] = 'A'
        graph.graph['source'] = 'test'
        graph.graph['version'] = '1.0'
        result = sdnetwork_from_graph(graph, network_type='semi-directed')
        assert result.get_network_attribute('source') == 'test'
        assert result.get_network_attribute('version') == '1.0'

    def test_empty_graph_from_nx_multigraph(self) -> None:
        """Convert empty NetworkX MultiGraph."""
        graph = nx.MultiGraph()
        result = sdnetwork_from_graph(graph, network_type='semi-directed')
        assert isinstance(result, SemiDirectedPhyNetwork)
        assert result.number_of_nodes() == 0
        assert result.number_of_edges() == 0

    def test_node_order_preserved(self) -> None:
        """Test that node ordering is preserved."""
        mmgraph = MixedMultiGraph(undirected_edges=[(0, 1), (0, 2), (0, 3)])
        mmgraph._undirected.nodes[1]['label'] = 'A'
        mmgraph._undirected.nodes[2]['label'] = 'B'
        mmgraph._undirected.nodes[3]['label'] = 'C'
        result = sdnetwork_from_graph(mmgraph, network_type='semi-directed')
        assert result.get_label(1) == 'A'
        assert result.get_label(2) == 'B'
        assert result.get_label(3) == 'C'

    def test_invalid_graph_type_raises_error(self) -> None:
        """Test that invalid graph type raises TypeError."""
        with pytest.raises(TypeError):
            sdnetwork_from_graph("not a graph")  # type: ignore
