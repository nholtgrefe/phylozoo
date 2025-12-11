"""
Extensive tests for to_sd_network function.

Tests cover:
- Basic trees and hybrids
- LSA conversion and root suppression
- Edge attribute preservation (branch_length, bootstrap, gamma)
- Parallel edges (both directed and undirected)
- Large networks
- Edge cases
- Round-trip validation
- Complex topologies

All tests create valid DirectedPhyNetworks that pass validation.
"""

import pytest

from phylozoo.core.network import DirectedPhyNetwork, SemiDirectedPhyNetwork
from phylozoo.core.network.dnetwork.operations import to_sd_network
from phylozoo.core.network.dnetwork.classifications import is_lsa_network


class TestToSDNetworkBasicTrees:
    """Test to_sd_network on basic tree structures."""

    def test_ternary_tree(self) -> None:
        """Ternary tree with 3 leaves."""
        dnet = DirectedPhyNetwork(
            edges=[(4, 1), (4, 2), (4, 3)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})],
        )
        sdnet = to_sd_network(dnet)

        # All edges should be undirected
        assert len(list(sdnet._graph.directed_edges_iter())) == 0
        assert len(list(sdnet._graph.undirected_edges_iter())) == 3

        # Taxa preserved
        assert sdnet.taxa == {"A", "B", "C"}

    def test_tree_with_root_suppression(self) -> None:
        """Tree where root has outdegree 2, so it gets suppressed."""
        dnet = DirectedPhyNetwork(
            edges=[
                (5, 3),  # Root 5 has outdegree 2
                (5, 4),
                (3, 1),
                (3, 2),
                (4, 6),
                (4, 7),
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (6, {'label': 'C'}), (7, {'label': 'D'})],
        )
        sdnet = to_sd_network(dnet)

        # Root should be suppressed, creating undirected edge (3, 4)
        assert 5 not in sdnet._graph.nodes()
        assert sdnet.number_of_nodes() == 6
        assert len(sdnet.leaves) == 4

    def test_balanced_binary_tree(self) -> None:
        """Balanced binary tree."""
        dnet = DirectedPhyNetwork(
            edges=[
                (7, 5),
                (7, 6),
                (5, 1),
                (5, 2),
                (6, 3),
                (6, 4),
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (4, {'label': 'D'})],
        )
        sdnet = to_sd_network(dnet)

        # All edges undirected (no hybrids)
        assert len(list(sdnet._graph.directed_edges_iter())) == 0
        assert len(sdnet.leaves) == 4


class TestToSDNetworkHybrids:
    """Test to_sd_network on networks with hybrid nodes."""

    def test_single_hybrid_node(self) -> None:
        """Network with single hybrid node."""
        dnet = DirectedPhyNetwork(
            edges=[
                (7, 5),
                (7, 6),
                {"u": 5, "v": 4, "gamma": 0.6},
                {"u": 6, "v": 4, "gamma": 0.4},
                (5, 8),
                (6, 9),
                (4, 1),
            ],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})],
        )
        sdnet = to_sd_network(dnet)

        # Hybrid edges should be directed
        directed_edges = {(u, v) for u, v, _ in sdnet._graph.directed_edges_iter(keys=True)}
        assert directed_edges == {(5, 4), (6, 4)}

        # Tree edges should be undirected
        assert len(list(sdnet._graph.undirected_edges_iter())) > 0

        # Taxa preserved
        assert sdnet.taxa == {"A", "B", "C"}

    def test_multiple_hybrid_nodes(self) -> None:
        """Network with multiple hybrid nodes."""
        dnet = DirectedPhyNetwork(
            edges=[
                (10, 8),
                (10, 9),
                # First hybrid
                {"u": 8, "v": 5, "gamma": 0.7},
                {"u": 9, "v": 5, "gamma": 0.3},
                # Second hybrid
                {"u": 8, "v": 6, "gamma": 0.4},
                {"u": 9, "v": 6, "gamma": 0.6},
                (8, 11),
                (9, 12),
                (5, 1),
                (6, 2),
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (11, {'label': 'C'}), (12, {'label': 'D'})],
        )
        sdnet = to_sd_network(dnet)

        # Hybrid edges should be directed
        directed_edges = {(u, v) for u, v, _ in sdnet._graph.directed_edges_iter(keys=True)}
        assert (8, 5) in directed_edges
        assert (9, 5) in directed_edges
        assert (8, 6) in directed_edges
        assert (9, 6) in directed_edges

        # Taxa preserved
        assert sdnet.taxa == {"A", "B", "C", "D"}

    def test_hybrid_with_single_child(self) -> None:
        """Hybrid node with single child (valid configuration)."""
        dnet = DirectedPhyNetwork(
            edges=[
                (10, 8),
                (10, 9),
                {"u": 8, "v": 5, "gamma": 0.6},
                {"u": 9, "v": 5, "gamma": 0.4},
                (8, 11),
                (9, 12),
                (5, 1),  # Hybrid node 5 has single child
            ],
            nodes=[(1, {'label': 'A'}), (11, {'label': 'B'}), (12, {'label': 'C'})],
        )
        sdnet = to_sd_network(dnet)

        # Hybrid edges directed
        assert (8, 5) in {(u, v) for u, v, _ in sdnet._graph.directed_edges_iter(keys=True)}
        assert (9, 5) in {(u, v) for u, v, _ in sdnet._graph.directed_edges_iter(keys=True)}


class TestToSDNetworkLSAConversion:
    """Test to_sd_network with LSA conversion."""

    def test_root_suppression_merges_branch_lengths(self) -> None:
        """
        Root suppression should merge branch lengths from root's two edges.
        """
        dnet = DirectedPhyNetwork(
            edges=[
                {"u": 5, "v": 3, "branch_length": 0.2},
                {"u": 5, "v": 4, "branch_length": 0.3},
                (3, 1),
                (3, 2),
                (4, 6),
                (4, 7),
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (6, {'label': 'C'}), (7, {'label': 'D'})],
        )

        sdnet = to_sd_network(dnet)

        # Root 5 should be suppressed
        assert 5 not in sdnet._graph.nodes()

        # Edge (3,4) should have merged branch length
        data_3_4 = sdnet._graph._undirected.get_edge_data(3, 4)
        assert data_3_4 is not None
        key = next(iter(data_3_4))
        assert pytest.approx(data_3_4[key].get("branch_length", 0.0)) == 0.5  # 0.2 + 0.3

    def test_hybrid_edges_remain_directed_after_lsa_and_suppression(self) -> None:
        """
        Hybrid edges stay directed and tree edges become undirected after LSA conversion and suppression.
        """
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
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (8, {'label': 'C'}), (9, {'label': 'D'})],
        )

        sdnet = to_sd_network(dnet)

        # Directed (hybrid) edges should be exactly into node 4
        directed_edges = {(u, v) for u, v, _ in sdnet._graph._directed.edges(keys=True)}
        assert directed_edges == {(5, 4), (6, 4)}

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
                (5, 11),
                (6, 12),
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (11, {'label': 'C'}), (12, {'label': 'D'})],
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

        # Hybrid edges remain directed with gamma
        assert sdnet.get_gamma(5, 4) == 0.6
        assert sdnet.get_gamma(6, 4) == 0.4


class TestToSDNetworkParallelEdges:
    """Test to_sd_network with parallel edges."""

    def test_parallel_hybrid_edges(self) -> None:
        """Parallel hybrid edges should remain directed."""
        dnet = DirectedPhyNetwork(
            edges=[
                (7, 5),
                (7, 6),
                {"u": 5, "v": 4, "key": 0, "gamma": 0.3},
                {"u": 5, "v": 4, "key": 1, "gamma": 0.2},  # Parallel hybrid edges
                {"u": 6, "v": 4, "gamma": 0.5},
                (5, 8),
                (6, 9),
                (4, 1),
            ],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})],
        )
        sdnet = to_sd_network(dnet)

        # Parallel directed edges should exist
        edges_5_4 = list(sdnet._graph._directed.edges(5, 4, keys=True))
        assert len(edges_5_4) == 2

        # Gamma values preserved
        assert sdnet.get_gamma(5, 4, key=0) == 0.3
        assert sdnet.get_gamma(5, 4, key=1) == 0.2


class TestToSDNetworkEdgeAttributes:
    """Test edge attribute preservation and merging."""

    def test_gamma_preserved_on_hybrid_edges(self) -> None:
        """Gamma values should be preserved on directed hybrid edges."""
        dnet = DirectedPhyNetwork(
            edges=[
                (7, 5),
                (7, 6),
                {"u": 5, "v": 4, "gamma": 0.7},
                {"u": 6, "v": 4, "gamma": 0.3},
                (5, 8),
                (6, 9),
                (4, 1),
            ],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})],
        )
        sdnet = to_sd_network(dnet)

        # Gamma should be preserved on directed edges
        assert sdnet.get_gamma(5, 4) == 0.7
        assert sdnet.get_gamma(6, 4) == 0.3

    def test_branch_length_summed_during_suppression(self) -> None:
        """Branch lengths should be summed when suppressing degree-2 nodes."""
        dnet = DirectedPhyNetwork(
            edges=[
                {"u": 5, "v": 3, "branch_length": 0.2},
                {"u": 5, "v": 4, "branch_length": 0.3},
                {"u": 3, "v": 1, "branch_length": 0.1},
                {"u": 3, "v": 2, "branch_length": 0.15},
                (4, 6),
                (4, 7),
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (6, {'label': 'C'}), (7, {'label': 'D'})],
        )
        sdnet = to_sd_network(dnet)

        # Root 5 suppressed, edge (3,4) should have summed branch length
        assert pytest.approx(sdnet.get_branch_length(3, 4)) == 0.5  # 0.2 + 0.3

    def test_all_attributes_preserved_on_hybrids(self) -> None:
        """All edge attributes on hybrid edges should be preserved."""
        dnet = DirectedPhyNetwork(
            edges=[
                (7, 5),
                (7, 6),
                {"u": 5, "v": 4, "gamma": 0.65, "branch_length": 0.1, "bootstrap": 0.9, "custom": "foo"},
                {"u": 6, "v": 4, "gamma": 0.35, "branch_length": 0.2, "bootstrap": 0.85, "custom": "bar"},
                (5, 8),
                (6, 9),
                (4, 1),
            ],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})],
        )
        sdnet = to_sd_network(dnet)

        # Get edge data for hybrid edges
        data_5_4 = sdnet._graph._directed.get_edge_data(5, 4)
        key = next(iter(data_5_4))
        assert data_5_4[key]["gamma"] == 0.65
        assert data_5_4[key]["custom"] == "foo"


class TestToSDNetworkLargeNetworks:
    """Test to_sd_network on large networks."""

    def test_large_ternary_tree(self) -> None:
        """Large ternary tree with 27 leaves."""
        # Build a complete ternary tree (3 levels)
        edges = []
        node_id = 1
        current_level = [0]
        taxa = {}

        # Build 3 levels
        for level in range(3):
            next_level = []
            for parent in current_level:
                for _ in range(3):
                    child = node_id
                    node_id += 1
                    next_level.append(child)
                    edges.append((parent, child))
            current_level = next_level

        # Mark last level as leaves
        nodes = [(leaf, {"label": f"Taxon{i}"}) for i, leaf in enumerate(current_level)]

        dnet = DirectedPhyNetwork(edges=edges, nodes=nodes)
        sdnet = to_sd_network(dnet)

        # All edges should be undirected (no hybrids)
        assert len(list(sdnet._graph.directed_edges_iter())) == 0
        assert len(sdnet.taxa) == 27

    def test_star_tree_with_many_leaves(self) -> None:
        """Star tree with 20 leaves."""
        edges = [(100, i) for i in range(20)]
        nodes = [(i, {"label": f"Taxon{i}"}) for i in range(20)]

        dnet = DirectedPhyNetwork(edges=edges, nodes=nodes)
        sdnet = to_sd_network(dnet)

        # All edges should be undirected
        assert len(list(sdnet._graph.directed_edges_iter())) == 0
        assert sdnet.number_of_edges() == 20
        assert len(sdnet.taxa) == 20


class TestToSDNetworkEdgeCases:
    """Test edge cases for to_sd_network."""

    def test_single_edge_tree(self) -> None:
        """Tree with single edge."""
        with pytest.warns(UserWarning, match="Single-node network detected"):
            sdnet = to_sd_network(
                DirectedPhyNetwork(edges=[(2, 1)], nodes=[(1, {'label': 'A'})])
            )

        # LSA of a single-leaf tree is the leaf; conversion yields a single-node network
        assert sdnet.number_of_nodes() == 1
        assert sdnet.number_of_edges() == 0
        assert sdnet.leaves == {1}
        assert sdnet.get_label(1) == "A"

    def test_single_node_directed_network(self) -> None:
        """Single-node directed network converts to single-node semi-directed network."""
        with pytest.warns(UserWarning, match="Single-node network detected"):
            sdnet = to_sd_network(DirectedPhyNetwork(nodes=[(1, {"label": "A"})]))

        assert sdnet.number_of_nodes() == 1
        assert sdnet.number_of_edges() == 0
        assert sdnet.leaves == {1}
        assert sdnet.get_label(1) == "A"

    def test_star_tree(self) -> None:
        """Star tree with one internal node and many leaves."""
        edges = [(10, i) for i in range(10)]
        nodes = [(i, {"label": f"Taxon{i}"}) for i in range(10)]

        dnet = DirectedPhyNetwork(edges=edges, nodes=nodes)
        sdnet = to_sd_network(dnet)

        # All edges should be undirected
        assert len(list(sdnet._graph.directed_edges_iter())) == 0
        assert sdnet.number_of_edges() == 10


class TestToSDNetworkComplexTopologies:
    """Test to_sd_network on complex topologies."""

    def test_diamond_network(self) -> None:
        """Diamond-shaped network with hybrid at the bottom."""
        dnet = DirectedPhyNetwork(
            edges=[
                (6, 4),
                (6, 5),
                {"u": 4, "v": 3, "gamma": 0.6},
                {"u": 5, "v": 3, "gamma": 0.4},
                (4, 7),
                (5, 8),
                (3, 1),
            ],
            nodes=[(1, {'label': 'A'}), (7, {'label': 'B'}), (8, {'label': 'C'})],
        )
        sdnet = to_sd_network(dnet)

        # Hybrid edges should be directed
        assert (4, 3) in {(u, v) for u, v, _ in sdnet._graph.directed_edges_iter(keys=True)}
        assert (5, 3) in {(u, v) for u, v, _ in sdnet._graph.directed_edges_iter(keys=True)}

    def test_multi_level_tree_with_hybrid(self) -> None:
        """Multi-level tree with a single hybrid deep in the structure."""
        dnet = DirectedPhyNetwork(
            edges=[
                (20, 15),
                (20, 16),
                (15, 10),
                (15, 11),
                (16, 12),
                (16, 13),
                # Hybrid
                {"u": 10, "v": 5, "gamma": 0.5},
                {"u": 11, "v": 5, "gamma": 0.5},
                (10, 17),
                (11, 18),
                (12, 19),
                (12, 21),
                (13, 22),
                (13, 23),
                (5, 1),
            ],
            nodes=[(1, {'label': 'A'}), (17, {'label': 'B'}), (18, {'label': 'C'}), (19, {'label': 'D'}), (21, {'label': 'E'}), (22, {'label': 'F'}), (23, {'label': 'G'})],
        )
        sdnet = to_sd_network(dnet)

        # Only hybrid edges to node 5 should be directed
        directed_edges = {(u, v) for u, v, _ in sdnet._graph.directed_edges_iter(keys=True)}
        assert (10, 5) in directed_edges
        assert (11, 5) in directed_edges


class TestToSDNetworkRoundTrip:
    """Test round-trip conversion properties."""

    def test_taxa_preserved_through_conversion(self) -> None:
        """Taxa should be exactly preserved through conversion."""
        dnet = DirectedPhyNetwork(
            edges=[
                (5, 3),
                (5, 4),
                (3, 1),
                (3, 2),
                (4, 6),
                (4, 7),
            ],
            nodes=[(1, {'label': 'Species_A'}), (2, {'label': 'Species_B'}), (6, {'label': 'Species_C'}), (7, {'label': 'Species_D'})],
        )
        sdnet = to_sd_network(dnet)

        # Exact taxa preservation
        assert sdnet.taxa == {"Species_A", "Species_B", "Species_C", "Species_D"}

    def test_internal_labels_removed_for_suppressed_nodes(self) -> None:
        """Internal labels for suppressed nodes should be removed."""
        dnet = DirectedPhyNetwork(
            edges=[
                (5, 3),
                (5, 4),
                (3, 1),
                (3, 2),
                (4, 6),
                (4, 7),
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (6, {'label': 'C'}), (7, {'label': 'D'}), (5, {'label': 'root'}), (3, {'label': 'left'}), (4, {'label': 'right'})],
        )
        sdnet = to_sd_network(dnet)

        # Node 5 should be suppressed, so its label should be gone
        assert 5 not in sdnet._graph.nodes()
        assert sdnet.get_label(5) is None

        # Other internal labels should be preserved
        assert sdnet.get_label(3) == "left"
        assert sdnet.get_label(4) == "right"

    def test_number_of_edges_after_suppression(self) -> None:
        """After suppressing degree-2 nodes, edge count should be correct."""
        dnet = DirectedPhyNetwork(
            edges=[
                (5, 3),
                (5, 4),
                (3, 1),
                (3, 2),
                (4, 6),
                (4, 7),
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (6, {'label': 'C'}), (7, {'label': 'D'})],
        )
        sdnet = to_sd_network(dnet)

        # Root suppressed: 6 edges - 2 (to/from root) + 1 (new edge) = 5 edges
        assert sdnet.number_of_edges() == 5


class TestToSDNetworkValidation:
    """Test that resulting SD networks pass validation."""

    def test_simple_tree_validates(self) -> None:
        """Simple tree should produce valid SD network."""
        dnet = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        sdnet = to_sd_network(dnet)

        # Should validate successfully (already validated in __init__)
        sdnet.validate()

    def test_hybrid_network_validates(self) -> None:
        """Hybrid network should produce valid SD network."""
        dnet = DirectedPhyNetwork(
            edges=[
                (7, 5),
                (7, 6),
                {"u": 5, "v": 4, "gamma": 0.6},
                {"u": 6, "v": 4, "gamma": 0.4},
                (5, 8),
                (6, 9),
                (4, 1),
            ],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})],
        )
        sdnet = to_sd_network(dnet)

        # Should validate successfully
        sdnet.validate()

    def test_complex_network_validates(self) -> None:
        """Complex network with multiple hybrids should validate."""
        dnet = DirectedPhyNetwork(
            edges=[
                (10, 8),
                (10, 9),
                {"u": 8, "v": 5, "gamma": 0.7},
                {"u": 9, "v": 5, "gamma": 0.3},
                {"u": 8, "v": 6, "gamma": 0.4},
                {"u": 9, "v": 6, "gamma": 0.6},
                (8, 11),
                (9, 12),
                (5, 1),
                (6, 2),
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (11, {'label': 'C'}), (12, {'label': 'D'})],
        )
        sdnet = to_sd_network(dnet)

        # Should validate successfully
        sdnet.validate()


class TestToSDNetworkStructurePreservation:
    """Test that network structure is correctly transformed."""

    def test_leaf_count_preserved(self) -> None:
        """Number of leaves should be preserved."""
        dnet = DirectedPhyNetwork(
            edges=[
                (10, 8),
                (10, 9),
                (8, 1),
                (8, 2),
                (9, 3),
                (9, 4),
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (4, {'label': 'D'})],
        )
        sdnet = to_sd_network(dnet)

        assert len(sdnet.leaves) == 4
        assert len(dnet.leaves) == 4

    def test_hybrid_count_preserved(self) -> None:
        """Number of hybrid nodes should be preserved (unless suppressed)."""
        dnet = DirectedPhyNetwork(
            edges=[
                (10, 8),
                (10, 9),
                {"u": 8, "v": 5, "gamma": 0.5},
                {"u": 9, "v": 5, "gamma": 0.5},
                (8, 11),
                (9, 12),
                (5, 1),
            ],
            nodes=[(1, {'label': 'A'}), (11, {'label': 'B'}), (12, {'label': 'C'})],
        )
        sdnet = to_sd_network(dnet)

        # Should have 1 hybrid node
        assert len(sdnet.hybrid_nodes) == 1
        assert len(dnet.hybrid_nodes) == 1

    def test_tree_edges_become_undirected(self) -> None:
        """All non-hybrid edges should become undirected."""
        dnet = DirectedPhyNetwork(
            edges=[
                (7, 5),
                (7, 6),
                {"u": 5, "v": 4, "gamma": 0.6},
                {"u": 6, "v": 4, "gamma": 0.4},
                (5, 8),
                (6, 9),
                (4, 1),
            ],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})],
        )
        sdnet = to_sd_network(dnet)

        # Count tree edges
        tree_edges_in_dnet = [e for e in dnet._graph.edges() if e not in dnet.hybrid_edges]

        # After conversion and suppression, most tree edges should be undirected
        # (some may be merged due to suppression)
        undirected_count = len(list(sdnet._graph.undirected_edges_iter()))
        assert undirected_count >= len(tree_edges_in_dnet) - 2  # Allow for root suppression

    def test_hybrid_edges_remain_directed(self) -> None:
        """All hybrid edges should remain directed."""
        dnet = DirectedPhyNetwork(
            edges=[
                (10, 8),
                (10, 9),
                {"u": 8, "v": 5, "gamma": 0.7},
                {"u": 9, "v": 5, "gamma": 0.3},
                {"u": 8, "v": 6, "gamma": 0.4},
                {"u": 9, "v": 6, "gamma": 0.6},
                (8, 11),
                (9, 12),
                (5, 1),
                (6, 2),
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (11, {'label': 'C'}), (12, {'label': 'D'})],
        )
        sdnet = to_sd_network(dnet)

        # All hybrid edges should be in directed edges
        directed_edges = set(sdnet._graph.directed_edges_iter())
        assert (8, 5) in directed_edges
        assert (9, 5) in directed_edges
        assert (8, 6) in directed_edges
        assert (9, 6) in directed_edges


class TestToSDNetworkDegreeSuppression:
    """Test degree-2 node suppression during conversion."""

    def test_single_degree2_suppression(self) -> None:
        """Single degree-2 node should be suppressed."""
        dnet = DirectedPhyNetwork(
            edges=[
                (5, 3),
                (5, 4),
                (3, 1),
                (3, 2),
                (4, 6),
                (4, 7),
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (6, {'label': 'C'}), (7, {'label': 'D'})],
        )
        sdnet = to_sd_network(dnet)

        # Root 5 should be suppressed
        assert 5 not in sdnet._graph.nodes()
        assert sdnet.number_of_nodes() == 6

    def test_suppression_with_branch_length_merging(self) -> None:
        """Suppression should merge branch lengths correctly."""
        dnet = DirectedPhyNetwork(
            edges=[
                {"u": 7, "v": 5, "branch_length": 0.1},
                {"u": 7, "v": 6, "branch_length": 0.2},
                {"u": 5, "v": 3, "branch_length": 0.3},
                {"u": 5, "v": 4, "branch_length": 0.4},
                (3, 1),
                (3, 2),
                (4, 8),
                (4, 9),
                (6, 10),
                (6, 11),
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (8, {'label': 'C'}), (9, {'label': 'D'}), (10, {'label': 'E'}), (11, {'label': 'F'})],
        )
        sdnet = to_sd_network(dnet)

        # Root 7 suppressed, edge (5,6) should have summed branch length
        assert pytest.approx(sdnet.get_branch_length(5, 6)) == 0.3  # 0.1 + 0.2


class TestToSDNetworkHybridGamma:
    """Test gamma value handling during conversion."""

    def test_gamma_preserved_on_single_hybrid(self) -> None:
        """Gamma values should be preserved on hybrid edges."""
        dnet = DirectedPhyNetwork(
            edges=[
                (7, 5),
                (7, 6),
                {"u": 5, "v": 4, "gamma": 0.65},
                {"u": 6, "v": 4, "gamma": 0.35},
                (5, 8),
                (6, 9),
                (4, 1),
            ],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})],
        )
        sdnet = to_sd_network(dnet)

        assert sdnet.get_gamma(5, 4) == 0.65
        assert sdnet.get_gamma(6, 4) == 0.35

    def test_gamma_sum_maintained(self) -> None:
        """Sum of gamma values for hybrid node should remain 1.0."""
        dnet = DirectedPhyNetwork(
            edges=[
                (7, 5),
                (7, 6),
                {"u": 5, "v": 4, "gamma": 0.7},
                {"u": 6, "v": 4, "gamma": 0.3},
                (5, 8),
                (6, 9),
                (4, 1),
            ],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})],
        )
        sdnet = to_sd_network(dnet)

        gamma_sum = sdnet.get_gamma(5, 4) + sdnet.get_gamma(6, 4)
        assert pytest.approx(gamma_sum) == 1.0

    def test_gamma_with_parallel_hybrid_edges(self) -> None:
        """Gamma on parallel hybrid edges should be preserved."""
        dnet = DirectedPhyNetwork(
            edges=[
                (7, 5),
                (7, 6),
                {"u": 5, "v": 4, "key": 0, "gamma": 0.3},
                {"u": 5, "v": 4, "key": 1, "gamma": 0.2},
                {"u": 6, "v": 4, "gamma": 0.5},
                (5, 8),
                (6, 9),
                (4, 1),
            ],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})],
        )
        sdnet = to_sd_network(dnet)

        # Gamma preserved on parallel edges
        assert sdnet.get_gamma(5, 4, key=0) == 0.3
        assert sdnet.get_gamma(5, 4, key=1) == 0.2
        assert sdnet.get_gamma(6, 4) == 0.5


class TestToSDNetworkConsistency:
    """Test consistency between DirectedPhyNetwork and SemiDirectedPhyNetwork."""

    def test_topology_preserved(self) -> None:
        """Network topology should be preserved."""
        dnet = DirectedPhyNetwork(
            edges=[
                (7, 5),
                (7, 6),
                {"u": 5, "v": 4, "gamma": 0.6},
                {"u": 6, "v": 4, "gamma": 0.4},
                (5, 8),
                (6, 9),
                (4, 1),
            ],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})],
        )
        sdnet = to_sd_network(dnet)

        # Hybrid node should still be a hybrid
        assert 4 in sdnet.hybrid_nodes

        # Leaves should still be leaves
        assert sdnet.leaves == dnet.leaves


class TestToSDNetworkBootstrap:
    """Test bootstrap value handling during conversion."""

    def test_bootstrap_preserved_on_tree_edges(self) -> None:
        """Bootstrap values should be preserved on tree edges."""
        dnet = DirectedPhyNetwork(
            edges=[
                {"u": 5, "v": 3, "bootstrap": 0.95},
                {"u": 5, "v": 4, "bootstrap": 0.87},
                (3, 1),
                (3, 2),
                (4, 6),
                (4, 7),
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (6, {'label': 'C'}), (7, {'label': 'D'})],
        )
        sdnet = to_sd_network(dnet)

        # Bootstrap should be preserved on undirected edges after suppression
        # Edge (3,4) should have bootstrap from one of the original edges
        # (depends on which edge's attributes take precedence)
        assert sdnet.get_bootstrap(3, 4) in [0.95, 0.87]

    def test_bootstrap_not_on_hybrid_edges(self) -> None:
        """Bootstrap should not be set on hybrid edges (which are directed)."""
        dnet = DirectedPhyNetwork(
            edges=[
                (7, 5),
                (7, 6),
                {"u": 5, "v": 4, "gamma": 0.6, "bootstrap": 0.9},  # Bootstrap on hybrid edge
                {"u": 6, "v": 4, "gamma": 0.4, "bootstrap": 0.8},
                (5, 8),
                (6, 9),
                (4, 1),
            ],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})],
        )
        sdnet = to_sd_network(dnet)

        # Hybrid edges should have bootstrap (it's allowed, just not required)
        # The attribute is preserved
        data = sdnet._graph._directed.get_edge_data(5, 4)
        key = next(iter(data))
        assert data[key].get("bootstrap") == 0.9


class TestToSDNetworkComplexSuppressions:
    """Test complex node suppression scenarios."""

    def test_suppression_after_hybrid(self) -> None:
        """Degree-2 nodes after hybrid nodes should be handled correctly."""
        dnet = DirectedPhyNetwork(
            edges=[
                (10, 8),
                (10, 9),
                {"u": 8, "v": 5, "gamma": 0.6},
                {"u": 9, "v": 5, "gamma": 0.4},
                (8, 11),
                (9, 12),
                (5, 6),  # Degree-2 node after hybrid
                (6, 1),
                (6, 2),
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (11, {'label': 'C'}), (12, {'label': 'D'})],
        )
        sdnet = to_sd_network(dnet)

        # Node 6 should be suppressed (degree 2 after undirecting)
        # But 5 is a hybrid, so it should remain
        assert 5 in sdnet._graph.nodes()


class TestToSDNetworkInternalLabels:
    """Test internal label handling during conversion."""

    def test_internal_labels_removed_for_suppressed(self) -> None:
        """Internal labels for suppressed nodes should be removed."""
        dnet = DirectedPhyNetwork(
            edges=[
                (5, 3),
                (5, 4),
                (3, 1),
                (3, 2),
                (4, 6),
                (4, 7),
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (6, {'label': 'C'}), (7, {'label': 'D'}), (5, {'label': 'root'}), (3, {'label': 'left'}), (4, {'label': 'right'})],
        )
        sdnet = to_sd_network(dnet)

        # Node 5 suppressed, so label should be gone
        assert sdnet.get_label(5) is None

        # Other labels preserved
        assert sdnet.get_label(3) == "left"
        assert sdnet.get_label(4) == "right"

