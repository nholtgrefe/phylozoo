"""
Tests for network attributes in DirectedPhyNetwork.
"""

import pytest

from phylozoo.core.network import DirectedPhyNetwork


class TestDirectedPhyNetworkAttributes:
    """Tests for network attributes in DirectedPhyNetwork."""

    def test_attributes_initialization(self):
        """Test that attributes can be set during initialization."""
        attributes = {'source': 'file.nex', 'version': '1.0', 'created': '2024-01-01'}
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})],
            attributes=attributes
        )

        assert net.get_network_attribute() == attributes

    def test_get_network_attribute(self):
        """Test getting a network attribute."""
        attributes = {'source': 'file.nex', 'version': '1.0'}
        net = DirectedPhyNetwork(
            edges=[(3, 1)],
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
        net = DirectedPhyNetwork(
            edges=[(3, 1)],
            nodes=[(1, {'label': 'A'})]
        )

        assert net.get_network_attribute() == {}
        assert net.get_network_attribute('any') is None

    def test_attributes_preserved_on_copy(self):
        """Test that attributes are preserved when copying."""
        attributes = {'source': 'file.nex', 'version': '1.0'}
        net = DirectedPhyNetwork(
            edges=[(3, 1)],
            nodes=[(1, {'label': 'A'})],
            attributes=attributes
        )
        net_copy = net.copy()

        assert net_copy.get_network_attribute() == attributes
        assert net_copy.get_network_attribute('source') == 'file.nex'

    def test_attributes_from_graph(self):
        """Test that attributes are preserved when creating network from graph."""
        from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
        from phylozoo.core.network.dnetwork.conversions import dnetwork_from_graph

        graph = DirectedMultiGraph(
            edges=[(3, 1)],
            attributes={'source': 'file.nex', 'version': '1.0'}
        )
        graph.add_node(1, label='A')

        net = dnetwork_from_graph(graph)

        assert net.get_network_attribute() == {'source': 'file.nex', 'version': '1.0'}
        assert net.get_network_attribute('source') == 'file.nex'

