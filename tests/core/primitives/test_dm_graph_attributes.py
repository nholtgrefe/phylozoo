"""
Tests for graph attributes in DirectedMultiGraph.
"""

import pytest

from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph


class TestDirectedMultiGraphAttributes:
    """Tests for graph attributes in DirectedMultiGraph."""

    def test_attributes_initialization(self):
        """Test that attributes can be set during initialization."""
        attributes = {'source': 'file.nex', 'version': '1.0', 'created': '2024-01-01'}
        graph = DirectedMultiGraph(edges=[(1, 2)], attributes=attributes)

        assert graph._graph.graph == attributes

    def test_get_graph_attribute(self):
        """Test getting a graph attribute via direct access."""
        attributes = {'source': 'file.nex', 'version': '1.0'}
        graph = DirectedMultiGraph(edges=[(1, 2)], attributes=attributes)

        assert graph._graph.graph.get('source') == 'file.nex'
        assert graph._graph.graph.get('version') == '1.0'
        assert graph._graph.graph.get('nonexistent') is None
        assert graph._graph.graph.get('nonexistent', 'default') == 'default'

    def test_graph_attributes_property(self):
        """Test accessing graph attributes directly."""
        attributes = {'source': 'file.nex', 'version': '1.0'}
        graph = DirectedMultiGraph(edges=[(1, 2)], attributes=attributes)

        attrs = graph._graph.graph.copy()
        assert attrs == attributes
        # Should be a copy, not the same object
        assert attrs is not graph._graph.graph

    def test_no_attributes(self):
        """Test that graph without attributes has empty dict."""
        graph = DirectedMultiGraph(edges=[(1, 2)])

        assert graph._graph.graph == {}
        assert graph._graph.graph.get('any') is None

    def test_set_graph_attribute(self):
        """Test setting a graph attribute in all underlying graphs."""
        graph = DirectedMultiGraph(edges=[(1, 2)])
        graph.set_graph_attribute('probability', 0.5)
        
        assert graph._graph.graph.get('probability') == 0.5
        assert graph._combined.graph.get('probability') == 0.5
        
        # Test updating the attribute
        graph.set_graph_attribute('probability', 0.8)
        assert graph._graph.graph.get('probability') == 0.8
        assert graph._combined.graph.get('probability') == 0.8

    def test_attributes_preserved_on_copy(self):
        """Test that attributes are preserved when copying."""
        attributes = {'source': 'file.nex', 'version': '1.0'}
        graph = DirectedMultiGraph(edges=[(1, 2)], attributes=attributes)
        graph_copy = graph.copy()

        assert graph_copy._graph.graph == attributes
        assert graph_copy._graph.graph.get('source') == 'file.nex'

    def test_warn_on_keyword_attribute(self):
        """Test that using Python keywords as attribute keys warns."""
        with pytest.warns(UserWarning, match="Graph attribute key .* is a Python keyword") as record:
            graph = DirectedMultiGraph(
                edges=[(1, 2)],
                attributes={'for': 'value', 'if': 'other'}
            )
        # Should have warned for both 'for' and 'if'
        assert len(record) == 2
        # Should still work despite warning
        assert graph._graph.graph.get('for') == 'value'
        assert graph._graph.graph.get('if') == 'other'

