"""
Tests for the ``to_sd_network`` conversion utility.
"""

from typing import Dict, Set, Tuple

import pytest

from phylozoo.core.network import DirectedPhyNetwork, SemiDirectedPhyNetwork
from phylozoo.core.network.dnetwork.operations import to_sd_network


class TestToSDNetwork:
    """Tests for converting directed phylogenetic networks to semi-directed networks."""

    def test_root_suppression_merges_branch_lengths(self) -> None:
        """
        A root with degree 2 is suppressed; branch lengths on its edges are merged (summed).
        """
        dnet = DirectedPhyNetwork(
            edges=[
                {"u": 3, "v": 1, "branch_length": 0.5},
                {"u": 3, "v": 2, "branch_length": 0.7},
            ],
            taxa={1: "A", 2: "B"},
        )

        sdnet = to_sd_network(dnet)
        assert isinstance(sdnet, SemiDirectedPhyNetwork)

        # Root (3) should be suppressed; only leaves remain with one undirected edge.
        assert 3 not in sdnet._graph.nodes()
        undirected_data = sdnet._graph._undirected.get_edge_data(1, 2)
        assert undirected_data is not None
        # Expect a single edge with merged branch_length = 0.5 + 0.7
        key = next(iter(undirected_data))
        assert pytest.approx(undirected_data[key].get("branch_length", 0.0)) == 1.2
        # Taxa preserved
        assert sdnet.taxa == {"A", "B"}

    def test_hybrid_edges_remain_directed_after_lsa_and_suppression(self) -> None:
        """
        Hybrid edges stay directed and tree edges become undirected after LSA conversion and suppression.
        """
        # Network where LSA is below the root; root will be suppressed after undirecting tree edges.
        dnet = DirectedPhyNetwork(
            edges=[
                (7, 5),
                (7, 6),
                {"u": 5, "v": 4, "gamma": 0.55},
                {"u": 6, "v": 4, "gamma": 0.45},
                (5, 8),
                (6, 9),
                (4, 10),
                (10, 1),
                (10, 2),
            ],
            taxa={1: "A", 2: "B", 8: "C", 9: "D"},
        )

        sdnet = to_sd_network(dnet)

        # Directed (hybrid) edges should be exactly into node 4
        directed_edges: Set[Tuple[int, int]] = {
            (u, v) for u, v, _ in sdnet._graph._directed.edges(keys=True)
        }
        assert directed_edges == {(5, 4), (6, 4)}

        # Undirected tree edges should connect remaining tree/leaf nodes (root 7 suppressed)
        undirected_edges: Set[Tuple[int, int]] = {
            (min(u, v), max(u, v)) for u, v in sdnet.tree_edges
        }
        expected_undirected = {
            (min(4, 10), max(4, 10)),
            (min(5, 6), max(5, 6)),
            (min(1, 10), max(1, 10)),
            (min(2, 10), max(2, 10)),
            (min(8, 5), max(8, 5)),
            (min(9, 6), max(9, 6)),
        }
        assert undirected_edges == expected_undirected

        # Taxa preserved
        assert sdnet.taxa == {"A", "B", "C", "D"}

    def test_branch_lengths_preserved_below_lsa(self) -> None:
        """
        After converting to the LSA network, edges below the LSA retain their attributes.
        """
        dnet = DirectedPhyNetwork(
            edges=[
                (7, 5),
                (7, 6),
                {"u": 5, "v": 4, "gamma": 0.6},
                {"u": 6, "v": 4, "gamma": 0.4},
                {"u": 4, "v": 10, "branch_length": 0.3},
                {"u": 10, "v": 1, "branch_length": 0.4},
                {"u": 10, "v": 2, "branch_length": 0.5},
                (5, 11),  # ensure tree node 5 has out-degree >= 2
                (6, 12),  # ensure tree node 6 has out-degree >= 2
            ],
            taxa={1: "A", 2: "B", 11: "C", 12: "D"},
        )

        sdnet = to_sd_network(dnet)

        # Edge (10,1) should be undirected with preserved branch_length
        data_10_1 = sdnet._graph._undirected.get_edge_data(10, 1)
        assert data_10_1 is not None
        key_10_1 = next(iter(data_10_1))
        assert pytest.approx(data_10_1[key_10_1].get("branch_length", 0.0)) == 0.4

        # Edge (10,2) branch_length preserved
        data_10_2 = sdnet._graph._undirected.get_edge_data(10, 2)
        assert data_10_2 is not None
        key_10_2 = next(iter(data_10_2))
        assert pytest.approx(data_10_2[key_10_2].get("branch_length", 0.0)) == 0.5

        # Hybrid edges remain directed
        directed_edges: Set[Tuple[int, int]] = {
            (u, v) for u, v, _ in sdnet._graph._directed.edges(keys=True)
        }
        assert directed_edges == {(5, 4), (6, 4)}



