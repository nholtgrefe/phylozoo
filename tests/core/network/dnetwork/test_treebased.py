"""
Tests for tree-based classification.
"""

from __future__ import annotations

import itertools

import pytest

from phylozoo.core.network import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.classifications import has_parallel_edges, is_binary, is_treebased
from phylozoo.core.network.dnetwork.features import omnians

from tests.fixtures.directed_networks import (
    DTREE_EMPTY,
    DTREE_SINGLE_NODE,
    DTREE_SMALL_BINARY,
    LEVEL_1_DNETWORK_SINGLE_HYBRID_BINARY,
    LEVEL_1_DNETWORK_TWO_BLOBS,
    LEVEL_2_DNETWORK_NESTED_HYBRIDS,
    LEVEL_2_DNETWORK_NON_TREEBASED,
    LEVEL_2_DNETWORK_SINGLE_BLOB,
    LEVEL_3_DNETWORK_CHAIN_HYBRIDS,
)


def _brute_omnian_treebased(network: DirectedPhyNetwork) -> bool:
    """
    Check the omnian corollary by brute force: for all S ⊆ U, |children(S)| >= |S|.

    Parameters
    ----------
    network : DirectedPhyNetwork
        Binary directed phylogenetic network without parallel edges.

    Returns
    -------
    bool
        True iff for all S ⊆ U (U = omnians), |children(S)| >= |S|.
    """
    U = sorted(omnians(network))
    children = {u: set(network.children(u)) for u in U}

    for r in range(1, len(U) + 1):
        for subset in itertools.combinations(U, r):
            neigh: set[object] = set()
            for u in subset:
                neigh |= children[u]
            if len(neigh) < len(subset):
                return False

    return True


class TestIsTreebased:
    """Tests for `is_treebased`."""

    def test_treebased_matches_bruteforce_on_small_binary_network(self) -> None:
        """Check `is_treebased` matches brute-force omnian condition on a small binary network."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 11),
                (10, 12),
                (11, 20),
                (11, 21),
                (12, 20),
                (12, 21),
                (20, 1),
                (21, 2),
            ],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"})],
        )
        assert is_binary(net)
        assert not has_parallel_edges(net)
        assert omnians(net) == {11, 12}
        expected = _brute_omnian_treebased(net)
        assert is_treebased(net) == expected

    @pytest.mark.parametrize(
        "edges,nodes",
        [
            ([(3, 1), (3, 2)], [(1, {"label": "A"}), (2, {"label": "B"})]),
            ([], []),
            ([], [(1, {"label": "A"})]),
        ],
    )
    def test_treebased_trivial_cases(
        self, edges: list[tuple[int, int]], nodes: list[tuple[int, dict]]
    ) -> None:
        """Trivial cases (tree, empty, single-node) are tree-based."""
        net = DirectedPhyNetwork(edges=edges, nodes=nodes)
        if net.number_of_nodes() <= 1:
            assert is_treebased(net) is True
        else:
            assert is_binary(net)
            assert is_treebased(net) is True

    @pytest.mark.parametrize(
        "fixture",
        [
            DTREE_EMPTY,
            DTREE_SINGLE_NODE,
            DTREE_SMALL_BINARY,
            LEVEL_1_DNETWORK_SINGLE_HYBRID_BINARY,
            LEVEL_1_DNETWORK_TWO_BLOBS,
            LEVEL_2_DNETWORK_NESTED_HYBRIDS,
            LEVEL_2_DNETWORK_SINGLE_BLOB,
            LEVEL_3_DNETWORK_CHAIN_HYBRIDS,
        ],
    )
    def test_treebased_fixture_networks(self, fixture: DirectedPhyNetwork) -> None:
        """Fixture networks that are binary and tree-based return True."""
        if fixture.number_of_nodes() == 0:
            assert is_treebased(fixture) is True
            return
        if not is_binary(fixture) or has_parallel_edges(fixture):
            pytest.skip("Fixture is non-binary or has parallel edges")
        assert is_treebased(fixture) is True

    def test_treebased_returns_false_for_non_treebased_network(self) -> None:
        """A network that violates the omnian condition is not tree-based."""
        net = LEVEL_2_DNETWORK_NON_TREEBASED
        assert is_binary(net)
        assert not has_parallel_edges(net)
        assert omnians(net) == {11, 12, 13, 14, 52, 53}
        assert _brute_omnian_treebased(net) is False
        assert is_treebased(net) is False

    def test_treebased_raises_for_non_binary(self) -> None:
        """Raises PhyloZooNotImplementedError for non-binary networks."""
        from phylozoo.utils.exceptions import PhyloZooNotImplementedError

        # Root has out-degree 3 (non-binary)
        net = DirectedPhyNetwork(
            edges=[
                (10, 5),
                (10, 6),
                (10, 7),
                (5, 1),
                (5, 4),
                (6, 2),
                (6, 4),
                (7, 3),
                (7, 4),
                (4, 8),
            ],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (3, {"label": "C"}), (8, {"label": "D"})],
        )
        assert not is_binary(net)
        with pytest.raises(PhyloZooNotImplementedError, match="non-binary"):
            is_treebased(net)

    def test_treebased_raises_for_parallel_edges(self) -> None:
        """Raises PhyloZooNotImplementedError for networks with parallel edges."""
        from phylozoo.utils.exceptions import PhyloZooNotImplementedError

        net = DirectedPhyNetwork(
            edges=[(10, 5), (10, 5), (5, 1)],
            nodes=[(1, {"label": "A"})],
        )
        assert has_parallel_edges(net)
        with pytest.raises(PhyloZooNotImplementedError, match="parallel"):
            is_treebased(net)

