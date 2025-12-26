import pytest
from typing import List

from phylozoo.core.network.dnetwork.derivations import subnetwork

from tests.fixtures.directed_networks import (
    DTREE_LARGE_BINARY,
    LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE,
    LEVEL_1_DNETWORK_PARALLEL_EDGES_HYBRID,
    LEVEL_2_DNETWORK_NON_LSA,
)

FIXTURES = [
    ("large_tree", DTREE_LARGE_BINARY),
    ("two_hybrids", LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE),
    ("parallel_edges_hybrid", LEVEL_1_DNETWORK_PARALLEL_EDGES_HYBRID),
    ("level2_non_lsa", LEVEL_2_DNETWORK_NON_LSA),
]


@pytest.mark.parametrize("name,network", FIXTURES)
def test_subnetwork_builds_and_contains_requested_taxa(name: str, network) -> None:
    """Ensure subnetwork builds for a selection of fixtures and contains the requested taxa.

    This is a smoke test that verifies the function doesn't raise and that
    the returned subnetwork preserves the requested leaf labels.
    """
    taxa_list: List[str] = sorted(network.taxa)
    assert taxa_list, "fixture must expose some taxa"

    # choose up to 3 leaves for a compact induced subnetwork
    taxa = taxa_list[:3]

    sub = subnetwork(network, taxa)

    # subnetwork must contain all requested taxa (labels)
    assert set(taxa).issubset(set(sub.taxa)), f"Requested taxa {taxa} not found in subnetwork for {name}"

    # must have at least one node and be a valid network object
    assert sub.number_of_nodes() >= 1
    assert sub.number_of_edges() >= 0
