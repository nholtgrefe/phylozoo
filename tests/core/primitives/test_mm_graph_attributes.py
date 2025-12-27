"""
Tests for graph attributes in MixedMultiGraph.
"""

import pytest

from phylozoo.core.primitives.m_multigraph import MixedMultiGraph


class TestMixedMultiGraphAttributes:
    """Tests for graph attributes in MixedMultiGraph."""

    def test_attributes_initialization(self):
        """Test that attributes can be set during initialization."""
        attributes = {'source': 'file.nex', 'version': '1.0', 'created': '2024-01-01'}
        graph = MixedMultiGraph(undirected_edges=[(1, 2)], attributes=attributes)

        assert graph._directed.graph == attributes

    def test_get_graph_attribute(self):
        """Test getting a graph attribute via direct access."""
        attributes = {'source': 'file.nex', 'version': '1.0'}
        graph = MixedMultiGraph(undirected_edges=[(1, 2)], attributes=attributes)

        assert graph._directed.graph.get('source') == 'file.nex'
        assert graph._directed.graph.get('version') == '1.0'
        assert graph._directed.graph.get('nonexistent') is None
        assert graph._directed.graph.get('nonexistent', 'default') == 'default'

    def test_graph_attributes_property(self):
        """Test accessing graph attributes directly."""
        attributes = {'source': 'file.nex', 'version': '1.0'}
        graph = MixedMultiGraph(undirected_edges=[(1, 2)], attributes=attributes)

        attrs = graph._directed.graph.copy()
        assert attrs == attributes
        # Should be a copy, not the same object
        assert attrs is not graph._directed.graph

    def test_no_attributes(self):
        """Test that graph without attributes has empty dict."""
        graph = MixedMultiGraph(undirected_edges=[(1, 2)])

        assert graph._directed.graph == {}
        assert graph._directed.graph.get('any') is None

    def test_attributes_preserved_on_copy(self):
        """Test that attributes are preserved when copying."""
        attributes = {'source': 'file.nex', 'version': '1.0'}
        graph = MixedMultiGraph(undirected_edges=[(1, 2)], attributes=attributes)
        graph_copy = graph.copy()

        assert graph_copy._directed.graph == attributes
        assert graph_copy._directed.graph.get('source') == 'file.nex'

    def test_attributes_on_both_graphs(self):
        """Test that attributes are set on both directed and undirected graphs."""
        attributes = {'source': 'file.nex', 'version': '1.0'}
        graph = MixedMultiGraph(undirected_edges=[(1, 2)], attributes=attributes)

        # Both should have the same attributes
        assert graph._directed.graph == attributes
        assert graph._undirected.graph == attributes

    def test_warn_on_keyword_attribute(self):
        """Test that using Python keywords as attribute keys warns."""
        with pytest.warns(UserWarning, match="Graph attribute key .* is a Python keyword") as record:
            graph = MixedMultiGraph(
                undirected_edges=[(1, 2)],
                attributes={'for': 'value', 'if': 'other'}
            )
        # Should have warned for both 'for' and 'if'
        assert len(record) == 2
        # Should still work despite warning
        assert graph._directed.graph.get('for') == 'value'
        assert graph._directed.graph.get('if') == 'other'

