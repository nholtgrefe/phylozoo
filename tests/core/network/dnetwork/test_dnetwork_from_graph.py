"""
Tests for dnetwork_from_graph function.
"""

import pytest
import networkx as nx

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.conversions import dnetwork_from_graph
from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
from phylozoo.utils.validation import no_validation


class TestDnetworkFromGraph:
    """Test conversion from various graph types to DirectedPhyNetwork."""

    def test_empty_graph_from_dmgraph(self) -> None:
        """Convert empty DirectedMultiGraph."""
        dmgraph = DirectedMultiGraph()
        result = dnetwork_from_graph(dmgraph)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_nodes() == 0
        assert result.number_of_edges() == 0

    def test_empty_graph_from_nx_digraph(self) -> None:
        """Convert empty NetworkX DiGraph."""
        graph = nx.DiGraph()
        result = dnetwork_from_graph(graph)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_nodes() == 0
        assert result.number_of_edges() == 0

    def test_single_node_no_attributes_from_dmgraph(self) -> None:
        """Convert DirectedMultiGraph with single node without attributes."""
        dmgraph = DirectedMultiGraph()
        dmgraph.add_node(1)
        result = dnetwork_from_graph(dmgraph)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_nodes() == 1
        assert result.number_of_edges() == 0

    def test_single_node_with_label_from_dmgraph(self) -> None:
        """Convert DirectedMultiGraph with single node with label."""
        dmgraph = DirectedMultiGraph()
        dmgraph.add_node(1, label='A')
        result = dnetwork_from_graph(dmgraph)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_nodes() == 1
        assert result.get_label(1) == 'A'

    def test_single_node_with_label_from_nx_digraph(self) -> None:
        """Convert NetworkX DiGraph with single node with label."""
        graph = nx.DiGraph()
        graph.add_node(1)
        graph.nodes[1]['label'] = 'A'
        result = dnetwork_from_graph(graph)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_nodes() == 1
        assert result.get_label(1) == 'A'

    def test_simple_tree_from_dmgraph(self) -> None:
        """Convert simple tree structure from DirectedMultiGraph."""
        dmgraph = DirectedMultiGraph(edges=[(3, 1), (3, 2)])
        dmgraph._graph.nodes[1]['label'] = 'A'
        dmgraph._graph.nodes[2]['label'] = 'B'
        result = dnetwork_from_graph(dmgraph)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_nodes() == 3
        assert result.number_of_edges() == 2
        assert result.get_label(1) == 'A'
        assert result.get_label(2) == 'B'
        assert result.root_node == 3

    def test_simple_tree_from_nx_digraph(self) -> None:
        """Convert simple tree structure from NetworkX DiGraph."""
        graph = nx.DiGraph()
        graph.add_edge(3, 1)
        graph.add_edge(3, 2)
        graph.nodes[1]['label'] = 'A'
        graph.nodes[2]['label'] = 'B'
        result = dnetwork_from_graph(graph)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_nodes() == 3
        assert result.number_of_edges() == 2
        assert result.get_label(1) == 'A'
        assert result.get_label(2) == 'B'
        assert result.root_node == 3

    def test_edges_with_attributes_from_dmgraph(self) -> None:
        """Convert DirectedMultiGraph with edges that have attributes."""
        dmgraph = DirectedMultiGraph(edges=[
            {'u': 3, 'v': 1, 'branch_length': 0.5, 'bootstrap': 0.9}
        ])
        dmgraph._graph.nodes[1]['label'] = 'A'
        result = dnetwork_from_graph(dmgraph)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_edges() == 1
        assert result.get_branch_length(3, 1) == 0.5
        assert result.get_bootstrap(3, 1) == 0.9

    def test_edges_with_attributes_from_nx_digraph(self) -> None:
        """Convert NetworkX DiGraph with edges that have attributes."""
        graph = nx.DiGraph()
        graph.add_edge(3, 1, branch_length=0.5, bootstrap=0.9)
        graph.nodes[1]['label'] = 'A'
        result = dnetwork_from_graph(graph)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_edges() == 1
        assert result.get_branch_length(3, 1) == 0.5
        assert result.get_bootstrap(3, 1) == 0.9

    def test_graph_attributes_preserved_from_dmgraph(self) -> None:
        """Test that graph attributes are preserved from DirectedMultiGraph."""
        dmgraph = DirectedMultiGraph(
            edges=[(3, 1)],
            attributes={'source': 'test', 'version': '1.0'}
        )
        dmgraph._graph.nodes[1]['label'] = 'A'
        result = dnetwork_from_graph(dmgraph)
        assert result.get_network_attribute('source') == 'test'
        assert result.get_network_attribute('version') == '1.0'

    def test_graph_attributes_preserved_from_nx_digraph(self) -> None:
        """Test that graph attributes are preserved from NetworkX DiGraph."""
        graph = nx.DiGraph()
        graph.add_edge(3, 1)
        graph.nodes[1]['label'] = 'A'
        graph.graph['source'] = 'test'
        graph.graph['version'] = '1.0'
        result = dnetwork_from_graph(graph)
        assert result.get_network_attribute('source') == 'test'
        assert result.get_network_attribute('version') == '1.0'

    def test_node_attributes_preserved(self) -> None:
        """Test that node attributes are preserved."""
        graph = nx.DiGraph()
        graph.add_edge(3, 1)
        graph.nodes[1]['label'] = 'A'
        graph.nodes[1]['custom_attr'] = 'custom_value'
        result = dnetwork_from_graph(graph)
        # Custom attributes are preserved in the graph structure
        assert result.get_label(1) == 'A'

    def test_parallel_edges_from_dmgraph(self) -> None:
        """Convert DirectedMultiGraph with parallel edges."""
        dmgraph = DirectedMultiGraph(edges=[
            {'u': 1, 'v': 2, 'key': 0, 'branch_length': 0.5},
            {'u': 1, 'v': 2, 'key': 1, 'branch_length': 0.5}  # Same branch_length for parallel edges
        ])
        with no_validation():
            result = dnetwork_from_graph(dmgraph)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_edges() == 2

    def test_parallel_edges_from_nx_multidigraph(self) -> None:
        """Convert NetworkX MultiDiGraph with parallel edges."""
        graph = nx.MultiDiGraph()
        graph.add_edge(1, 2, branch_length=0.5)
        graph.add_edge(1, 2, branch_length=0.5)
        with no_validation():
            result = dnetwork_from_graph(graph)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_edges() == 2

    def test_simple_tree_from_nx_multidigraph(self) -> None:
        """Convert simple tree structure from NetworkX MultiDiGraph."""
        graph = nx.MultiDiGraph()
        graph.add_edge(3, 1)
        graph.add_edge(3, 2)
        graph.nodes[1]['label'] = 'A'
        graph.nodes[2]['label'] = 'B'
        result = dnetwork_from_graph(graph)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_nodes() == 3
        assert result.number_of_edges() == 2
        assert result.get_label(1) == 'A'
        assert result.get_label(2) == 'B'
        assert result.root_node == 3

    def test_edges_with_attributes_from_nx_multidigraph(self) -> None:
        """Convert NetworkX MultiDiGraph with edges that have attributes."""
        graph = nx.MultiDiGraph()
        graph.add_edge(3, 1, branch_length=0.5, bootstrap=0.9)
        graph.nodes[1]['label'] = 'A'
        result = dnetwork_from_graph(graph)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_edges() == 1
        assert result.get_branch_length(3, 1) == 0.5
        assert result.get_bootstrap(3, 1) == 0.9

    def test_graph_attributes_preserved_from_nx_multidigraph(self) -> None:
        """Test that graph attributes are preserved from NetworkX MultiDiGraph."""
        graph = nx.MultiDiGraph()
        graph.add_edge(3, 1)
        graph.nodes[1]['label'] = 'A'
        graph.graph['source'] = 'test'
        graph.graph['version'] = '1.0'
        result = dnetwork_from_graph(graph)
        assert result.get_network_attribute('source') == 'test'
        assert result.get_network_attribute('version') == '1.0'

    def test_network_validity(self) -> None:
        """Ensure resulting network is valid."""
        dmgraph = DirectedMultiGraph(edges=[(3, 1), (3, 2)])
        dmgraph._graph.nodes[1]['label'] = 'A'
        dmgraph._graph.nodes[2]['label'] = 'B'
        result = dnetwork_from_graph(dmgraph)
        # Should not raise an exception
        result.validate()

    def test_invalid_network_raises_error(self) -> None:
        """Test that invalid network structure raises ValueError."""
        # Create a graph that violates DirectedPhyNetwork constraints
        graph = nx.DiGraph()
        graph.add_edge(1, 2)
        graph.add_edge(2, 1)  # Creates a cycle, invalid for DAG
        graph.nodes[1]['label'] = 'A'
        graph.nodes[2]['label'] = 'B'
        with pytest.raises(ValueError):
            dnetwork_from_graph(graph)

    def test_node_order_preserved(self) -> None:
        """Test that node ordering is preserved."""
        dmgraph = DirectedMultiGraph(edges=[(1, 2), (2, 3), (2, 4)])
        dmgraph._graph.nodes[3]['label'] = 'A'
        dmgraph._graph.nodes[4]['label'] = 'B'
        result = dnetwork_from_graph(dmgraph)
        assert result.get_label(3) == 'A'
        assert result.get_label(4) == 'B'

    def test_invalid_graph_type_raises_error(self) -> None:
        """Test that invalid graph type raises TypeError."""
        with pytest.raises(TypeError):
            dnetwork_from_graph("not a graph")  # type: ignore
