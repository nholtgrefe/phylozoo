from __future__ import annotations

from typing import Any

import pytest

from phylozoo.core.primitives.m_multigraph.base import MixedMultiGraph
from phylozoo.core.primitives.m_multigraph.transformations import subgraph


def test_subgraph_basic() -> None:
    G = MixedMultiGraph()
    # add node attributes
    G.add_node(1, label='A')
    G.add_node(2, label='B')
    G.add_node(3, label='C')

    k1 = G.add_undirected_edge(1, 2, weight=1.0)
    k2 = G.add_directed_edge(2, 3, weight=2.0)

    H = subgraph(G, [1, 2])

    assert 1 in H
    assert 2 in H
    assert 3 not in H

    # Node attributes preserved
    # Attributes are stored on underlying graphs; check undirected for node 1 and directed/undirected for 2
    assert H._undirected.nodes[1]['label'] == 'A'
    assert H._undirected.nodes[2]['label'] == 'B' or H._directed.nodes[2].get('label', 'B') == 'B'

    undirected_edges = list(H.undirected_edges_iter(keys=True, data=True))
    assert undirected_edges == [(1, 2, k1, {'weight': 1.0})]


def test_subgraph_parallel_edges_and_keys() -> None:
    G = MixedMultiGraph()
    G.add_node('a', label='A')
    G.add_node('b', label='B')

    ku1 = G.add_undirected_edge('a', 'b', label='first')
    ku2 = G.add_undirected_edge('a', 'b', label='second')

    G2 = MixedMultiGraph()
    G2.add_node('x', label='X')
    G2.add_node('y', label='Y')

    kd1 = G.add_directed_edge('x', 'y', label='dfirst')
    kd2 = G.add_directed_edge('x', 'y', label='dsecond')

    H = subgraph(G, ['a', 'b'])

    # Both parallel undirected edges preserved
    u_edges = sorted(list(H.undirected_edges_iter(keys=True, data=True)), key=lambda x: x[2])
    assert len(u_edges) == 2
    assert u_edges[0][2] == ku1
    assert u_edges[1][2] == ku2
    labels = {d['label'] for (_, _, _, d) in u_edges}
    assert labels == {'first', 'second'}

    # Also ensure directed parallel edges preserved when selecting their nodes
    H2 = subgraph(G, ['x', 'y'])
    d_edges = sorted(list(H2.directed_edges_iter(keys=True, data=True)), key=lambda x: x[2])
    assert len(d_edges) == 2
    assert d_edges[0][2] == kd1
    assert d_edges[1][2] == kd2
    dlabels = {d['label'] for (_, _, _, d) in d_edges}
    assert dlabels == {'dfirst', 'dsecond'}


def test_subgraph_empty_nodes_returns_empty_graph() -> None:
    G = MixedMultiGraph()
    G.add_directed_edge(1, 2)
    H = subgraph(G, [])
    assert len(H) == 0
    assert list(H.edges_iter(keys=True, data=True)) == []


def test_subgraph_missing_node_raises() -> None:
    G = MixedMultiGraph()
    G.add_directed_edge(1, 2)
    with pytest.raises(ValueError):
        subgraph(G, [1, 3])
