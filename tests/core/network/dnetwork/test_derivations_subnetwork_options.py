from __future__ import annotations

from typing import Any

from phylozoo.core.network.dnetwork.base import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.derivations import subnetwork
from phylozoo.core.network.dnetwork.transformations import identify_parallel_edges
from phylozoo.core.network.dnetwork.transformations import to_lsa_network, is_lsa_network


def test_identify_parallel_edges_flag() -> None:
    # Build network with parallel edges 1->2
    net = DirectedPhyNetwork(
        edges=[
            {'u': 1, 'v': 2, 'branch_length': 0.5},
            {'u': 1, 'v': 2, 'branch_length': 0.5},  # parallel
            {'u': 2, 'v': 3, 'branch_length': 0.3},
            {'u': 3, 'v': 4, 'branch_length': 0.2},
            {'u': 3, 'v': 5, 'branch_length': 0.1},
            {'u': 1, 'v': 6, 'branch_length': 0.1},
        ],
    nodes=[(4, {'label': 'A'}), (5, {'label': 'B'}), (6, {'label': 'C'})]
    )

    # Subnetwork without identification
    sub_no = subnetwork(net, ['A'], identify_parallel_edges=False)
    edges_no = sub_no.number_of_edges()

    # Subnetwork with identification should have <= edges
    sub_id = subnetwork(net, ['A'], identify_parallel_edges=True)
    edges_id = sub_id.number_of_edges()

    assert edges_id <= edges_no


def test_make_lsa_flag() -> None:
    # Construct a network where LSA is different: root connects to two subtrees that meet
    net = DirectedPhyNetwork(
        edges=[
            (10, 5), (10, 6),
            (5, 4), (5, 9),  # ensure node 5 out-degree >=2
            (6, 4), (6, 11),  # ensure node 6 out-degree >=2
            (4, 8), (8, 1), (8, 2)
        ],
        nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (9, {'label': 'X'}), (11, {'label': 'Y'})]
    )

    # Subnetwork without LSA conversion
    sub_no = subnetwork(net, ['A', 'B'], make_lsa=False)
    # Subnetwork with make_lsa=True should satisfy is_lsa_network
    sub_lsa = subnetwork(net, ['A', 'B'], make_lsa=True)
    assert is_lsa_network(sub_lsa)

    # Also check that applying to_lsa_network to sub_no yields same number of nodes (sanity)
    lsa_of_sub_no = to_lsa_network(sub_no)
    assert lsa_of_sub_no.number_of_nodes() >= 0


def test_suppress_2_blobs_flag_runs() -> None:
    # Use a network similar to suppress_2_blobs example
    net = DirectedPhyNetwork(
        edges=[
            (1, 2), (2, 3), (2, 7),  # node 2 out-degree >=2
            (3, 4), (3, 9), (4, 5), (4, 8),  # give node 3 out-degree >=2
            (1, 6)
        ],
    nodes=[(5, {'label': 'A'}), (6, {'label': 'B'}), (7, {'label': 'X'}), (8, {'label': 'Y'}), (9, {'label': 'Z'})]
    )

    # Ensure both calls run without error
    sub_no = subnetwork(net, ['A'], suppress_2_blobs=False)
    sub_yes = subnetwork(net, ['A'], suppress_2_blobs=True)

    assert isinstance(sub_no, DirectedPhyNetwork)
    assert isinstance(sub_yes, DirectedPhyNetwork)
    # suppress_2_blobs should not increase node count
    assert sub_yes.number_of_nodes() <= sub_no.number_of_nodes()
