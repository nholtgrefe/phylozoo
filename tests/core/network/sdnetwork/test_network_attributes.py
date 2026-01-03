"""
Tests for network attributes in MixedPhyNetwork and SemiDirectedPhyNetwork.
"""

import pytest

from phylozoo.core.network.sdnetwork import MixedPhyNetwork, SemiDirectedPhyNetwork


class TestMixedPhyNetworkAttributes:
    """Tests for network attributes in MixedPhyNetwork."""

    def test_attributes_initialization(self):
        """Test that attributes can be set during initialization."""
        attributes = {'source': 'file.nex', 'version': '1.0', 'created': '2024-01-01'}
        net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})],
            attributes=attributes
        )

        assert net.get_network_attribute() == attributes

    def test_get_network_attribute(self):
        """Test getting a network attribute."""
        attributes = {'source': 'file.nex', 'version': '1.0'}
        net = MixedPhyNetwork(
            undirected_edges=[(3, 1)],
            nodes=[(1, {'label': 'A'})],
            attributes=attributes
        )

        assert net.get_network_attribute('source') == 'file.nex'
        assert net.get_network_attribute('version') == '1.0'
        assert net.get_network_attribute('nonexistent') is None
        # Test getting all attributes
        assert net.get_network_attribute() == attributes

    def test_no_attributes(self):
        """Test that network without attributes has empty dict."""
        net = MixedPhyNetwork(
            undirected_edges=[(3, 1)],
            nodes=[(1, {'label': 'A'})]
        )

        assert net.get_network_attribute() == {}
        assert net.get_network_attribute('any') is None

    def test_attributes_preserved_on_copy(self):
        """Test that attributes are preserved when copying."""
        attributes = {'source': 'file.nex', 'version': '1.0'}
        net = MixedPhyNetwork(
            undirected_edges=[(3, 1)],
            nodes=[(1, {'label': 'A'})],
            attributes=attributes
        )
        net_copy = net.copy()

        assert net_copy.get_network_attribute() == attributes
        assert net_copy.get_network_attribute('source') == 'file.nex'

    def test_attributes_from_graph(self):
        """Test that attributes are preserved when creating network from graph."""
        from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
        from phylozoo.core.network.sdnetwork.conversions import sdnetwork_from_graph

        graph = MixedMultiGraph(
            undirected_edges=[(3, 1)],
            attributes={'source': 'file.nex', 'version': '1.0'}
        )
        graph.add_node(1, label='A')

        net = sdnetwork_from_graph(graph, network_type='mixed')

        assert net.get_network_attribute() == {'source': 'file.nex', 'version': '1.0'}
        assert net.get_network_attribute('source') == 'file.nex'


class TestSemiDirectedPhyNetworkAttributes:
    """Tests for network attributes in SemiDirectedPhyNetwork."""

    def test_attributes_initialization(self):
        """Test that attributes can be set during initialization."""
        attributes = {'source': 'file.nex', 'version': '1.0', 'created': '2024-01-01'}
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})],
            attributes=attributes
        )

        assert net.get_network_attribute() == attributes

    def test_get_network_attribute(self):
        """Test getting a network attribute."""
        attributes = {'source': 'file.nex', 'version': '1.0'}
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1)],
            nodes=[(1, {'label': 'A'})],
            attributes=attributes
        )

        assert net.get_network_attribute('source') == 'file.nex'
        assert net.get_network_attribute('version') == '1.0'
        assert net.get_network_attribute('nonexistent') is None
        # Test getting all attributes
        assert net.get_network_attribute() == attributes

    def test_attributes_inherited_from_mixed(self):
        """Test that SemiDirectedPhyNetwork inherits attribute methods from MixedPhyNetwork."""
        attributes = {'source': 'file.nex', 'version': '1.0'}
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1)],
            nodes=[(1, {'label': 'A'})],
            attributes=attributes
        )

        # Should have the same methods as MixedPhyNetwork
        assert hasattr(net, 'get_network_attribute')
        assert net.get_network_attribute() == attributes

    def test_attributes_from_graph(self):
        """Test that attributes are preserved when creating network from graph."""
        from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
        from phylozoo.core.network.sdnetwork.conversions import sdnetwork_from_graph

        graph = MixedMultiGraph(
            undirected_edges=[(3, 1)],
            attributes={'source': 'file.nex', 'version': '1.0'}
        )
        graph.add_node(1, label='A')

        net = sdnetwork_from_graph(graph, network_type='semi-directed')

        assert net.get_network_attribute() == {'source': 'file.nex', 'version': '1.0'}
        assert net.get_network_attribute('source') == 'file.nex'

