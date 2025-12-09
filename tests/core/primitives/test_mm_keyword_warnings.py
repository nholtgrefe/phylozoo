"""
Tests for keyword identifier warnings in MixedMultiGraph.
"""

import pytest

from phylozoo.core.primitives.m_multigraph.mm_graph import MixedMultiGraph


def test_keyword_node_id_warns_mm() -> None:
    """Using a Python keyword as node id should warn."""
    with pytest.warns(UserWarning, match="Python keyword"):
        g = MixedMultiGraph()
        g.add_undirected_edge("for", "x")


def test_keyword_edge_key_warns_mm() -> None:
    """Using a Python keyword as an edge key should warn."""
    with pytest.warns(UserWarning, match="Python keyword"):
        g = MixedMultiGraph()
        g.add_directed_edge("a", "b", key="class")

