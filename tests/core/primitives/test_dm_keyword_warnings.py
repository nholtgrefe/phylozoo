"""
Tests for keyword identifier warnings in DirectedMultiGraph.
"""

import pytest

from phylozoo.core.primitives.d_multigraph.dm_graph import DirectedMultiGraph


def test_keyword_node_id_warns_dm() -> None:
    """Using a Python keyword as node id should warn."""
    with pytest.warns(UserWarning, match="Python keyword"):
        g = DirectedMultiGraph()
        g.add_edge("for", "x")


def test_keyword_edge_key_warns_dm() -> None:
    """Using a Python keyword as an edge key should warn."""
    with pytest.warns(UserWarning, match="Python keyword"):
        g = DirectedMultiGraph()
        g.add_edge("a", "b", key="class")

