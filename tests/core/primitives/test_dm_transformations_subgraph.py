from __future__ import annotations

from typing import Any

import pytest

from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph
from phylozoo.core.primitives.d_multigraph.transformations import subgraph


def test_subgraph_basic() -> None:
    G = DirectedMultiGraph()
    # add node attributes
    G.add_node(1, label='A')
    G.add_node(2, label='B')
    G.add_node(3, label='C')

    k1 = G.add_edge(1, 2, weight=1.0)
    k2 = G.add_edge(2, 3, weight=2.0)

    H = subgraph(G, [1, 2])

    assert 1 in H
    assert 2 in H
    assert 3 not in H

    # Node attributes preserved
    assert H._graph.nodes[1]['label'] == 'A'
    assert H._graph.nodes[2]['label'] == 'B'

    edges = list(H.edges_iter(keys=True, data=True))
    assert edges == [(1, 2, k1, {'weight': 1.0})]


def test_subgraph_parallel_edges_and_keys() -> None:
    G = DirectedMultiGraph()
    G.add_node('u', label='U')
    G.add_node('v', label='V')

    k1 = G.add_edge('u', 'v', label='first')
    k2 = G.add_edge('u', 'v', label='second')

    H = subgraph(G, ['u', 'v'])

    # Both parallel edges preserved
    edges = sorted(list(H.edges_iter(keys=True, data=True)), key=lambda x: x[2])
    assert len(edges) == 2
    assert edges[0][2] == k1
    assert edges[1][2] == k2
    # Attributes preserved
    labels = {d['label'] for (_, _, _, d) in edges}
    assert labels == {'first', 'second'}


def test_subgraph_empty_nodes_returns_empty_graph() -> None:
    G = DirectedMultiGraph()
    G.add_edge(1, 2)
    H = subgraph(G, [])
    assert len(H) == 0
    assert list(H.edges_iter(keys=True, data=True)) == []


def test_subgraph_missing_node_raises() -> None:
    G = DirectedMultiGraph()
    G.add_edge(1, 2)
    with pytest.raises(ValueError):
        subgraph(G, [1, 3])
