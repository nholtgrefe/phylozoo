from __future__ import annotations

from typing import Any

from phylozoo.core.network.dnetwork.derivations import subnetwork
from phylozoo.core.network.dnetwork.classifications import is_lsa_network

from tests.fixtures.directed_networks import (
    DTREE_MEDIUM_BINARY,
    DTREE_LARGE_BINARY,
    LEVEL_1_DNETWORK_SINGLE_HYBRID,
    LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE,
    LEVEL_1_DNETWORK_PARALLEL_EDGES,
    LEVEL_1_DNETWORK_PARALLEL_EDGES_HYBRID,
    LEVEL_2_DNETWORK_MANY_PARALLEL_EDGES,
    LEVEL_1_DNETWORK_NON_LSA,
    LEVEL_2_DNETWORK_NON_LSA,
    LEVEL_2_DNETWORK_PARALLEL_HYBRIDS,
)


def test_subnetwork_identify_parallel_edges_fixture() -> None:
    """Using a fixture with parallel edges, ensure identify_parallel_edges reduces or preserves edge count."""
    net = LEVEL_1_DNETWORK_PARALLEL_EDGES

    # pick a leaf label from the fixture
    taxa = ['A']

    sub_no = subnetwork(net, taxa, identify_parallel_edges=False)
    sub_id = subnetwork(net, taxa, identify_parallel_edges=True)

    assert sub_id.number_of_edges() <= sub_no.number_of_edges()


def test_subnetwork_make_lsa_fixture() -> None:
    """Using a non-LSA fixture, make_lsa should produce an LSA network."""
    net = LEVEL_1_DNETWORK_NON_LSA
    taxa = ['A', 'B']

    sub_no = subnetwork(net, taxa, make_lsa=False)
    sub_yes = subnetwork(net, taxa, make_lsa=True)

    assert isinstance(sub_no, type(net))
    assert isinstance(sub_yes, type(net))
    assert is_lsa_network(sub_yes)


def test_subnetwork_suppress_2_blobs_fixture() -> None:
    """Using a fixture with a top 2-blob, suppress_2_blobs should not increase nodes and should run."""
    net = LEVEL_1_DNETWORK_NON_LSA
    taxa = ['A']

    sub_no = subnetwork(net, taxa, suppress_2_blobs=False)
    sub_yes = subnetwork(net, taxa, suppress_2_blobs=True)

    assert isinstance(sub_no, type(net))
    assert isinstance(sub_yes, type(net))
    assert sub_yes.number_of_nodes() <= sub_no.number_of_nodes()


import pytest


@pytest.mark.parametrize(
    "network",
    [
        DTREE_MEDIUM_BINARY,
        DTREE_LARGE_BINARY,
        LEVEL_1_DNETWORK_SINGLE_HYBRID,
        LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE,
        LEVEL_1_DNETWORK_PARALLEL_EDGES,
        LEVEL_1_DNETWORK_PARALLEL_EDGES_HYBRID,
        LEVEL_2_DNETWORK_MANY_PARALLEL_EDGES,
        LEVEL_1_DNETWORK_NON_LSA,
        LEVEL_2_DNETWORK_NON_LSA,
        LEVEL_2_DNETWORK_PARALLEL_HYBRIDS,
    ],
)
def test_subnetwork_options_across_many_fixtures(network) -> None:
    """Run a few option combinations across many fixture networks to catch regressions."""
    # Choose 1-2 taxa deterministically from the network's taxa set
    taxa_list = sorted(network.taxa)
    if not taxa_list:
        pytest.skip("Fixture has no taxa")
    taxa = taxa_list[:2] if len(taxa_list) >= 2 else taxa_list

    # identify_parallel_edges should not increase the number of edges
    sub_no = subnetwork(network, taxa, identify_parallel_edges=False)
    sub_id = subnetwork(network, taxa, identify_parallel_edges=True)
    assert sub_id.number_of_edges() <= sub_no.number_of_edges()

    # make_lsa should produce an LSA-network
    sub_lsa = subnetwork(network, taxa, make_lsa=True)
    assert is_lsa_network(sub_lsa)

    # suppress_2_blobs should not increase node count
    sub_no_blob = subnetwork(network, taxa, suppress_2_blobs=False)
    sub_yes_blob = subnetwork(network, taxa, suppress_2_blobs=True)
    assert sub_yes_blob.number_of_nodes() <= sub_no_blob.number_of_nodes()
