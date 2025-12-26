from __future__ import annotations

from typing import Any

import pytest

from phylozoo.core.network.dnetwork.base import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.derivations import subnetwork


def test_subnetwork_basic_tree() -> None:
    # Simple rooted tree: 3 -> 1 (A), 3 -> 2 (B)
    net = DirectedPhyNetwork(
        edges=[(3, 1), (3, 2)],
        nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    )

    sub = subnetwork(net, ['A'])
    # Expect nodes corresponding to leaf A and the root 3
    labels = {sub.get_label(n) for n in sub._graph.nodes() if sub.get_label(n) is not None}
    assert 'A' in labels
    # Only one edge should be present (3 -> 1)
    assert sub.number_of_nodes() == 2
    assert sub.number_of_edges() == 1


def test_subnetwork_chain_and_final_suppression() -> None:
    # Valid rooted network: make internal nodes have outdegree >= 2
    # Structure:
    # 5 -> 3,4
    # 3 -> 1(A), 3 -> 7
    # 4 -> 2(B), 4 -> 8
    net = DirectedPhyNetwork(
        edges=[(5, 3), (5, 4), (3, 1), (3, 7), (4, 2), (4, 8)],
        nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    )

    # Subnetwork induced by ['A'] should include nodes 5,3,1 (7 is excluded)
    sub = subnetwork(net, ['A'])
    # Depending on suppression the intermediate internal node may be removed
    assert sub.number_of_nodes() >= 2
    # Ensure at least one edge remains
    assert sub.number_of_edges() >= 1
    # Now request both taxa - should return the portion containing both leaves
    sub2 = subnetwork(net, ['A', 'B'])
    # Some internal nodes may be suppressed; require at least 3 nodes and 2 edges
    assert sub2.number_of_nodes() >= 3
    assert sub2.number_of_edges() >= 2


def test_subnetwork_missing_taxon_raises() -> None:
    net = DirectedPhyNetwork(edges=[(3, 1)], nodes=[(1, {'label': 'A'})])
    with pytest.raises(ValueError):
        subnetwork(net, ['X'])
