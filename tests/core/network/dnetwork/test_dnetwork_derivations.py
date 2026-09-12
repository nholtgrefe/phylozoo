"""
Tests for dnetwork derivations functions.
"""

import numpy as np
import pytest

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.derivations import (
    _switchings,
    displayed_trees,
    displayed_splits,
    displayed_quartets,
    displayed_triplets,
    tree_of_blobs,
    distances,
    induced_splits,
    split_from_cutedge,
    partition_from_blob,
)
from phylozoo.core.split.base import Split
from phylozoo.core.primitives.partition import Partition
from phylozoo.core.network.dnetwork.features import blobs
from phylozoo.core.network.dnetwork.transformations import suppress_2_blobs


class TestTreeOfBlobs:
    """Test tree_of_blobs function for DirectedPhyNetwork."""

    def test_simple_tree_unchanged(self) -> None:
        """Tree with no blobs should remain unchanged."""
        dnet = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})]
        )
        result = tree_of_blobs(dnet)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_nodes() == 3
        assert result.number_of_edges() == 2
        assert result.get_label(1) == "A"
        assert result.get_label(2) == "B"

    def test_basic_functionality(self) -> None:
        """Test that tree_of_blobs runs without errors and returns correct type."""
        from tests.fixtures.directed_networks import LEVEL_1_DNETWORK_SINGLE_HYBRID

        result = tree_of_blobs(LEVEL_1_DNETWORK_SINGLE_HYBRID)
        assert isinstance(result, DirectedPhyNetwork)
        assert result is not LEVEL_1_DNETWORK_SINGLE_HYBRID  # Should return a new network

    def test_single_node_network(self) -> None:
        """Single node network should remain unchanged."""
        dnet = DirectedPhyNetwork(nodes=[(1, {"label": "A"})])
        result = tree_of_blobs(dnet)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_nodes() == 1
        assert result.get_label(1) == "A"

    def test_result_is_copy(self) -> None:
        """Function should return a new network, not modify the original."""
        dnet = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})]
        )
        original_nodes = dnet.number_of_nodes()
        original_edges = dnet.number_of_edges()

        result = tree_of_blobs(dnet)

        # Original should be unchanged
        assert dnet.number_of_nodes() == original_nodes
        assert dnet.number_of_edges() == original_edges
        # Result should be a different object
        assert result is not dnet

    def test_preserves_labels(self) -> None:
        """Leaf labels should be preserved."""
        from tests.fixtures.directed_networks import LEVEL_1_DNETWORK_SINGLE_HYBRID

        result = tree_of_blobs(LEVEL_1_DNETWORK_SINGLE_HYBRID)
        # Check that some label is preserved
        all_labels = [
            result.get_label(node) for node in result._graph.nodes() if result.get_label(node)
        ]
        assert len(all_labels) > 0

    def test_no_2_blobs_preserves_internal_blob_count(self) -> None:
        """When there are no 2-blobs, suppress_2_blobs doesn't change internal blob count."""
        from tests.fixtures.directed_networks import LEVEL_1_DNETWORK_SINGLE_HYBRID

        # Count internal blobs before suppress_2_blobs
        initial_blobs = blobs(LEVEL_1_DNETWORK_SINGLE_HYBRID, trivial=False, leaves=False)
        initial_count = len(initial_blobs)

        # Apply suppress_2_blobs
        suppressed_network = suppress_2_blobs(LEVEL_1_DNETWORK_SINGLE_HYBRID)

        # Count internal blobs after suppress_2_blobs
        after_suppress_blobs = blobs(suppressed_network, trivial=False, leaves=False)
        after_count = len(after_suppress_blobs)

        # Number of internal blobs should remain the same when there are no 2-blobs to suppress
        assert initial_count == after_count

        # But tree_of_blobs should still collapse all internal blobs
        result = tree_of_blobs(LEVEL_1_DNETWORK_SINGLE_HYBRID)
        final_blobs = blobs(result, trivial=False, leaves=False)
        # All internal blobs should be collapsed to single vertices
        assert len(final_blobs) == 0


class TestSwitchings:
    """Test _switchings function for DirectedPhyNetwork."""

    def test_tree_network_single_switching(self) -> None:
        """Tree network (no hybrid nodes) should yield exactly one switching."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})]
        )
        switchings = list(_switchings(net))
        assert len(switchings) == 1
        assert switchings[0].number_of_edges() == 2

        # With probability=True, should have probability 1.0 in all underlying graphs
        switchings_with_prob = list(_switchings(net, probability=True))
        assert len(switchings_with_prob) == 1
        assert switchings_with_prob[0]._graph.graph.get("probability") == 1.0
        assert switchings_with_prob[0]._combined.graph.get("probability") == 1.0

    def test_single_hybrid_two_parents(self) -> None:
        """Network with one hybrid node and two parent edges should yield two switchings."""
        from tests.fixtures.directed_networks import LEVEL_1_DNETWORK_SINGLE_HYBRID

        net = LEVEL_1_DNETWORK_SINGLE_HYBRID
        switchings = list(_switchings(net))
        assert len(switchings) == 2

        # Find the hybrid node
        hybrid_nodes = net.hybrid_nodes
        assert len(hybrid_nodes) == 1
        hybrid = next(iter(hybrid_nodes))

        # Each switching should have exactly one parent edge for the hybrid node
        for sw in switchings:
            parent_edges = list(sw.incident_parent_edges(hybrid, keys=True))
            assert len(parent_edges) == 1

    def test_single_hybrid_three_parents(self) -> None:
        """Network with one hybrid node and three parent edges should yield three switchings."""
        # Create a valid network with one hybrid and three parents
        net = DirectedPhyNetwork(
            edges=[
                (10, 5),
                (10, 6),
                (10, 7),  # Root to tree nodes
                (5, 4),
                (6, 4),
                (7, 4),  # Three parents to hybrid 4
                (4, 1),  # Hybrid to leaf
                (5, 8),
                (6, 9),
                (7, 11),  # Other children for tree nodes
            ],
            nodes=[
                (1, {"label": "A"}),
                (8, {"label": "B"}),
                (9, {"label": "C"}),
                (11, {"label": "D"}),
            ],
        )
        switchings = list(_switchings(net))
        assert len(switchings) == 3

        # Each switching should have exactly one parent edge for the hybrid node
        hybrid = 4
        for sw in switchings:
            parent_edges = list(sw.incident_parent_edges(hybrid, keys=True))
            assert len(parent_edges) == 1

    def test_two_hybrids_independent(self) -> None:
        """Network with two independent hybrid nodes should yield product of switchings."""
        # Create a valid network with two hybrids - ensure all tree nodes have out-degree >= 2
        net = DirectedPhyNetwork(
            edges=[
                (10, 5),
                (10, 6),  # Root to tree nodes
                (5, 4),
                (5, 17),
                (6, 4),
                (6, 18),
                (4, 11),
                (11, 1),
                (11, 12),  # Hybrid 4 with 2 parents, tree node 11
                (10, 8),
                (10, 9),  # Root to more tree nodes (8 and 9 are tree nodes)
                (8, 7),
                (8, 15),  # Tree node 8: one to hybrid 7, one to leaf
                (9, 7),
                (9, 16),  # Tree node 9: one to hybrid 7, one to leaf
                (7, 13),
                (13, 2),
                (13, 14),  # Hybrid 7 with 2 parents, tree node 13
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (12, {"label": "C"}),
                (14, {"label": "D"}),
                (15, {"label": "E"}),
                (16, {"label": "F"}),
                (17, {"label": "G"}),
                (18, {"label": "H"}),
            ],
        )
        switchings = list(_switchings(net))
        assert len(switchings) == 4  # 2 * 2 = 4

        # Each switching should have exactly one parent edge per hybrid
        for sw in switchings:
            parent_edges_4 = list(sw.incident_parent_edges(4, keys=True))
            parent_edges_7 = list(sw.incident_parent_edges(7, keys=True))
            assert len(parent_edges_4) == 1
            assert len(parent_edges_7) == 1

    def test_switchings_preserve_non_hybrid_edges(self) -> None:
        """Switchings should preserve all non-hybrid edges."""
        # Create a valid network with hybrid and non-hybrid edges
        # Ensure all tree nodes have out-degree >= 2
        net = DirectedPhyNetwork(
            edges=[
                (10, 5),
                (10, 6),  # Root to tree nodes
                (5, 4),
                (5, 17),
                (6, 4),
                (6, 18),
                (4, 1),  # Hybrid 4, tree nodes 5 and 6 have out-degree 2
                (10, 3),
                (3, 2),
                (3, 15),  # Non-hybrid tree node 3 with out-degree 2
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (15, {"label": "C"}),
                (17, {"label": "D"}),
                (18, {"label": "E"}),
            ],
        )
        switchings = list(_switchings(net))
        assert len(switchings) == 2

        # All switchings should have the non-hybrid edge (3, 2)
        for sw in switchings:
            assert sw.has_edge(3, 2)

    def test_switchings_are_copies(self) -> None:
        """Switchings should be independent copies, not references."""
        from tests.fixtures.directed_networks import LEVEL_1_DNETWORK_SINGLE_HYBRID

        net = LEVEL_1_DNETWORK_SINGLE_HYBRID
        switchings = list(_switchings(net))

        # Modify one switching and verify others are unchanged
        if len(switchings) > 0:
            original_edge_count = switchings[0].number_of_edges()
            # Remove a non-hybrid edge
            for u, v in switchings[0].edges():
                if switchings[0].indegree(v) < 2:  # Not a hybrid edge
                    switchings[0].remove_edge(u, v)
                    break
            # Other switchings should be unchanged
            for sw in switchings[1:]:
                assert sw.number_of_edges() == original_edge_count

    def test_probability_with_gamma_values(self) -> None:
        """Test probability calculation when gamma values are present."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 5),
                (10, 6),
                {"u": 5, "v": 4, "gamma": 0.6},
                {"u": 6, "v": 4, "gamma": 0.4},
                (4, 8),
                (8, 1),
                (8, 2),
                (5, 3),
                (6, 7),
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (7, {"label": "D"}),
            ],
        )
        switchings = list(_switchings(net, probability=True))
        assert len(switchings) == 2

        # Check probabilities in all underlying graphs
        probs = [sw._graph.graph.get("probability") for sw in switchings]
        assert 0.6 in probs
        assert 0.4 in probs
        # Probabilities should sum to 1.0
        assert abs(sum(probs) - 1.0) < 1e-10

        # Check that probability is set in all underlying graphs
        for sw in switchings:
            prob = sw._graph.graph.get("probability")
            assert prob is not None
            assert sw._combined.graph.get("probability") == prob

    def test_probability_without_gamma_values(self) -> None:
        """Test probability calculation when no gamma values are present."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 5),
                (10, 6),
                (5, 4),
                (6, 4),  # No gamma values
                (4, 8),
                (8, 1),
                (8, 2),
                (5, 3),
                (6, 7),
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (7, {"label": "D"}),
            ],
        )
        switchings = list(_switchings(net, probability=True))
        assert len(switchings) == 2

        # Without gamma, each edge should have probability 1/2 (indegree is 2)
        for sw in switchings:
            prob = sw._graph.graph.get("probability")
            assert prob is not None
            assert abs(prob - 0.5) < 1e-10
            # Check that probability is set in all underlying graphs
            assert sw._combined.graph.get("probability") == prob
        # Probabilities should sum to 1.0
        total_prob = sum(sw._graph.graph.get("probability") for sw in switchings)
        assert abs(total_prob - 1.0) < 1e-10

    def test_probability_with_multiple_hybrids(self) -> None:
        """Test probability calculation with multiple hybrid nodes."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 5),
                (10, 6),
                (10, 8),
                (10, 9),
                (5, 4),
                (5, 17),
                (6, 4),
                (6, 18),  # Hybrid 4 with 2 parents
                (8, 7),
                (8, 15),
                (9, 7),
                (9, 16),  # Hybrid 7 with 2 parents
                (4, 11),
                (11, 1),
                (11, 12),
                (7, 13),
                (13, 2),
                (13, 14),
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (12, {"label": "C"}),
                (14, {"label": "D"}),
                (15, {"label": "E"}),
                (16, {"label": "F"}),
                (17, {"label": "G"}),
                (18, {"label": "H"}),
            ],
        )
        switchings = list(_switchings(net, probability=True))
        assert len(switchings) == 4  # 2 * 2 = 4

        # Each switching should have a probability
        for sw in switchings:
            prob = sw._graph.graph.get("probability")
            assert prob is not None
            assert 0.0 < prob <= 1.0

        # Probabilities should sum to 1.0
        total_prob = sum(sw._graph.graph.get("probability") for sw in switchings)
        assert abs(total_prob - 1.0) < 1e-10

    def test_probability_false_no_attribute(self) -> None:
        """Test that probability=False does not add probability attribute."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 5),
                (10, 6),
                {"u": 5, "v": 4, "gamma": 0.6},
                {"u": 6, "v": 4, "gamma": 0.4},
                (4, 8),
                (8, 1),
                (8, 2),
                (5, 3),
                (6, 7),
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (7, {"label": "D"}),
            ],
        )
        switchings = list(_switchings(net, probability=False))

        # No probability attribute should be present
        for sw in switchings:
            assert "probability" not in sw._graph.graph


class TestDisplayedTrees:
    """Test cases for the displayed_trees function for DirectedPhyNetwork."""

    def test_tree_network_single_displayed_tree(self) -> None:
        """A tree network should yield exactly one displayed tree."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})]
        )
        trees = list(displayed_trees(net))
        assert len(trees) == 1
        assert trees[0].number_of_nodes() == 3
        assert trees[0].number_of_edges() == 2

    def test_single_hybrid_two_displayed_trees(self) -> None:
        """Network with one hybrid node should yield two displayed trees."""
        from tests.fixtures.directed_networks import LEVEL_1_DNETWORK_SINGLE_HYBRID

        net = LEVEL_1_DNETWORK_SINGLE_HYBRID
        trees = list(displayed_trees(net))
        assert len(trees) == 2

        # Each tree should be a valid tree (no hybrid nodes)
        for tree in trees:
            assert len(tree.hybrid_nodes) == 0

    def test_displayed_trees_are_trees(self) -> None:
        """All displayed trees should be trees (no hybrid nodes)."""
        net = DirectedPhyNetwork(
            edges=[(10, 5), (10, 6), (5, 4), (6, 4), (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (7, {"label": "D"}),
            ],
        )
        trees = list(displayed_trees(net))

        for tree in trees:
            assert len(tree.hybrid_nodes) == 0

    def test_probability_with_gamma_values(self) -> None:
        """Test probability calculation when gamma values are present."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 5),
                (10, 6),
                {"u": 5, "v": 4, "gamma": 0.6},
                {"u": 6, "v": 4, "gamma": 0.4},
                (4, 8),
                (8, 1),
                (8, 2),
                (5, 3),
                (6, 7),
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (7, {"label": "D"}),
            ],
        )
        trees = list(displayed_trees(net, probability=True))
        assert len(trees) == 2

        # Check probabilities
        probs = [tree.get_network_attribute("probability") for tree in trees]
        assert 0.6 in probs
        assert 0.4 in probs
        # Probabilities should sum to 1.0
        assert abs(sum(probs) - 1.0) < 1e-10

    def test_probability_without_gamma_values(self) -> None:
        """Test probability calculation when no gamma values are present."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 5),
                (10, 6),
                (5, 4),
                (6, 4),  # No gamma values
                (4, 8),
                (8, 1),
                (8, 2),
                (5, 3),
                (6, 7),
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (7, {"label": "D"}),
            ],
        )
        trees = list(displayed_trees(net, probability=True))
        assert len(trees) == 2

        # Without gamma, each tree should have probability 1/2 (indegree is 2)
        for tree in trees:
            prob = tree.get_network_attribute("probability")
            assert prob is not None
            assert abs(prob - 0.5) < 1e-10
        # Probabilities should sum to 1.0
        total_prob = sum(tree.get_network_attribute("probability") for tree in trees)
        assert abs(total_prob - 1.0) < 1e-10

    def test_probability_tree_network(self) -> None:
        """Tree network (no hybrid nodes) should have probability 1.0."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})]
        )
        trees = list(displayed_trees(net, probability=True))
        assert len(trees) == 1
        assert trees[0].get_network_attribute("probability") == 1.0

    def test_probability_false_no_attribute(self) -> None:
        """Test that probability=False does not add probability attribute."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 5),
                (10, 6),
                {"u": 5, "v": 4, "gamma": 0.6},
                {"u": 6, "v": 4, "gamma": 0.4},
                (4, 8),
                (8, 1),
                (8, 2),
                (5, 3),
                (6, 7),
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (7, {"label": "D"}),
            ],
        )
        trees = list(displayed_trees(net, probability=False))

        # No probability attribute should be present
        for tree in trees:
            assert tree.get_network_attribute("probability") is None

    def test_make_lsa_false_can_produce_unary_root(self) -> None:
        """Without make_lsa, a switching that prunes one side of the root
        yields a tree whose root has out-degree 1 (the degenerate root bug).

        Network layout (root=0, tree node=10, hybrid=3):
          0 → 10, 0 → 3
          10 → 3  (10 is second parent of hybrid 3)
          10 → 20 (B), 10 → 30 (C)
          3  → 40 (D)

        Switching that keeps (10→3) removes (0→3); then 3 is suppressed into
        10, and root 0 is left with a single child 10.
        """
        net = DirectedPhyNetwork(
            edges=[(0, 10), (0, 3), (10, 3), (10, 20), (10, 30), (3, 40)],
            nodes=[(20, {"label": "B"}), (30, {"label": "C"}), (40, {"label": "D"})],
        )
        trees = list(displayed_trees(net))
        assert len(trees) == 2
        root_outdegrees = [t.outdegree(t.root_node) for t in trees]
        assert 1 in root_outdegrees  # at least one tree has a unary root

    def test_make_lsa_true_removes_unary_root(self) -> None:
        """With make_lsa=True every displayed tree has a root of out-degree >= 2."""
        net = DirectedPhyNetwork(
            edges=[(0, 10), (0, 3), (10, 3), (10, 20), (10, 30), (3, 40)],
            nodes=[(20, {"label": "B"}), (30, {"label": "C"}), (40, {"label": "D"})],
        )
        trees = list(displayed_trees(net, make_lsa=True))
        assert len(trees) == 2
        for tree in trees:
            assert tree.outdegree(tree.root_node) >= 2
            assert tree.taxa == {"B", "C", "D"}

    def test_make_lsa_true_on_proper_tree_is_noop(self) -> None:
        """For a tree network make_lsa=True leaves the result unchanged."""
        net = DirectedPhyNetwork(
            edges=[(0, 10), (0, 20), (10, 30), (10, 40)],
            nodes=[(20, {"label": "A"}), (30, {"label": "B"}), (40, {"label": "C"})],
        )
        trees_plain = list(displayed_trees(net))
        trees_lsa = list(displayed_trees(net, make_lsa=True))
        assert len(trees_plain) == len(trees_lsa) == 1
        assert trees_plain[0].taxa == trees_lsa[0].taxa == {"A", "B", "C"}
        assert trees_plain[0].root_node == trees_lsa[0].root_node


class TestDistances:
    """Test distances function for DirectedPhyNetwork."""

    def test_simple_tree_no_branch_lengths(self) -> None:
        """Test distances for a simple tree without branch lengths."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})]
        )
        dm = distances(net, mode="shortest")
        assert len(dm) == 2
        assert dm.get_distance("A", "B") == 2.0  # Two edges, each with default length 1.0
        assert dm.get_distance("A", "A") == 0.0
        assert dm.get_distance("B", "B") == 0.0

    def test_simple_tree_with_branch_lengths(self) -> None:
        """Test distances for a simple tree with branch lengths."""
        net = DirectedPhyNetwork(
            edges=[{"u": 3, "v": 1, "branch_length": 0.5}, {"u": 3, "v": 2, "branch_length": 0.3}],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"})],
        )
        dm = distances(net, mode="shortest")
        assert len(dm) == 2
        assert dm.get_distance("A", "B") == 0.8  # 0.5 + 0.3
        assert dm.get_distance("A", "A") == 0.0
        assert dm.get_distance("B", "B") == 0.0

    def test_network_with_hybrid_shortest(self) -> None:
        """Test shortest distances for a network with hybrid node."""
        net = DirectedPhyNetwork(
            edges=[(10, 5), (10, 6), (5, 4), (6, 4), (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (7, {"label": "D"}),
            ],
        )
        dm = distances(net, mode="shortest")
        assert len(dm) == 4
        # All distances should be finite
        for i, taxon1 in enumerate(dm.labels):
            for j, taxon2 in enumerate(dm.labels):
                if i != j:
                    assert dm.get_distance(taxon1, taxon2) > 0
                    assert dm.get_distance(taxon1, taxon2) < np.inf
                else:
                    assert dm.get_distance(taxon1, taxon2) == 0.0

    def test_network_with_hybrid_longest(self) -> None:
        """Test longest distances for a network with hybrid node."""
        net = DirectedPhyNetwork(
            edges=[(10, 5), (10, 6), (5, 4), (6, 4), (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (7, {"label": "D"}),
            ],
        )
        dm = distances(net, mode="longest")
        assert len(dm) == 4
        # All distances should be finite and >= shortest
        shortest_dm = distances(net, mode="shortest")
        for i, taxon1 in enumerate(dm.labels):
            for j, taxon2 in enumerate(dm.labels):
                if i != j:
                    assert dm.get_distance(taxon1, taxon2) >= shortest_dm.get_distance(
                        taxon1, taxon2
                    )
                else:
                    assert dm.get_distance(taxon1, taxon2) == 0.0

    def test_network_with_hybrid_average(self) -> None:
        """Test average distances for a network with hybrid node."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 5),
                (10, 6),
                {"u": 5, "v": 4, "gamma": 0.6},
                {"u": 6, "v": 4, "gamma": 0.4},
                (4, 8),
                (8, 1),
                (8, 2),
                (5, 3),
                (6, 7),
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (7, {"label": "D"}),
            ],
        )
        dm = distances(net, mode="average")
        assert len(dm) == 4
        # All distances should be finite and between shortest and longest
        shortest_dm = distances(net, mode="shortest")
        longest_dm = distances(net, mode="longest")
        for i, taxon1 in enumerate(dm.labels):
            for j, taxon2 in enumerate(dm.labels):
                if i != j:
                    shortest = shortest_dm.get_distance(taxon1, taxon2)
                    longest = longest_dm.get_distance(taxon1, taxon2)
                    avg = dm.get_distance(taxon1, taxon2)
                    assert shortest <= avg <= longest
                else:
                    assert dm.get_distance(taxon1, taxon2) == 0.0

    def test_single_taxon(self) -> None:
        """Test distances for a network with a single taxon."""
        net = DirectedPhyNetwork(nodes=[(1, {"label": "A"})])
        dm = distances(net, mode="shortest")
        assert len(dm) == 1
        assert dm.get_distance("A", "A") == 0.0

    def test_symmetric_matrix(self) -> None:
        """Test that the distance matrix is symmetric."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})]
        )
        dm = distances(net, mode="shortest")
        assert dm.get_distance("A", "B") == dm.get_distance("B", "A")

    def test_all_modes_produce_valid_matrices(self) -> None:
        """Test that all modes produce valid distance matrices."""
        net = DirectedPhyNetwork(
            edges=[(10, 5), (10, 6), (5, 4), (6, 4), (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (7, {"label": "D"}),
            ],
        )
        for mode in ["shortest", "longest", "average"]:
            dm = distances(net, mode=mode)
            assert len(dm) == 4
            # Check symmetry
            for i, taxon1 in enumerate(dm.labels):
                for j, taxon2 in enumerate(dm.labels):
                    assert dm.get_distance(taxon1, taxon2) == dm.get_distance(taxon2, taxon1)
                    if i == j:
                        assert dm.get_distance(taxon1, taxon2) == 0.0


class TestInducedSplits:
    """Test induced_splits function for DirectedPhyNetwork."""

    def test_simple_tree(self) -> None:
        """Test induced_splits on a simple tree."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})]
        )
        splits = induced_splits(net)
        assert len(splits) >= 1
        # Should have at least one split (the edge between root and leaves)

    def test_tree_with_three_leaves(self) -> None:
        """Test induced_splits on a tree with three leaves."""
        net = DirectedPhyNetwork(
            edges=[(4, 1), (4, 2), (4, 3)],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (3, {"label": "C"})],
        )
        splits = induced_splits(net)
        # Should have splits for each cut-edge
        assert len(splits) >= 1

    def test_network_with_hybrid(self) -> None:
        """Test induced_splits on a network with hybrid node."""
        from tests.fixtures.directed_networks import LEVEL_1_DNETWORK_SINGLE_HYBRID

        splits = induced_splits(LEVEL_1_DNETWORK_SINGLE_HYBRID)
        assert len(splits) >= 1
        # Check that all splits cover all taxa
        all_taxa = set(LEVEL_1_DNETWORK_SINGLE_HYBRID.taxa)
        for split in splits:
            assert split.elements == all_taxa

    def test_empty_network(self) -> None:
        """Test induced_splits on empty network."""
        net = DirectedPhyNetwork()
        splits = induced_splits(net)
        assert len(splits) == 0

    def test_single_taxon_network(self) -> None:
        """Test induced_splits on network with single taxon."""
        net = DirectedPhyNetwork(edges=[(2, 1)], nodes=[(1, {"label": "A"})])
        splits = induced_splits(net)
        # Need at least 2 taxa for splits
        assert len(splits) == 0

    def test_tree_of_blobs_same_splits(self) -> None:
        """Test that tree-of-blobs has the same split system as original network."""
        from tests.fixtures.directed_networks import LEVEL_1_DNETWORK_SINGLE_HYBRID

        original_splits = induced_splits(LEVEL_1_DNETWORK_SINGLE_HYBRID)
        blob_tree = tree_of_blobs(LEVEL_1_DNETWORK_SINGLE_HYBRID)
        blob_tree_splits = induced_splits(blob_tree)

        # Split systems should be equal
        assert original_splits.splits == blob_tree_splits.splits

    def test_tree_of_blobs_same_splits_simple_tree(self) -> None:
        """Test that tree-of-blobs preserves splits for a simple tree."""
        net = DirectedPhyNetwork(
            edges=[(4, 1), (4, 2), (4, 3)],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (3, {"label": "C"})],
        )
        original_splits = induced_splits(net)
        blob_tree = tree_of_blobs(net)
        blob_tree_splits = induced_splits(blob_tree)

        # Split systems should be equal
        assert original_splits.splits == blob_tree_splits.splits

    def test_splits_cover_all_taxa(self) -> None:
        """Test that all splits cover all taxa."""
        from tests.fixtures.directed_networks import LEVEL_1_DNETWORK_SINGLE_HYBRID

        splits = induced_splits(LEVEL_1_DNETWORK_SINGLE_HYBRID)
        all_taxa = set(LEVEL_1_DNETWORK_SINGLE_HYBRID.taxa)

        for split in splits:
            assert split.elements == all_taxa

    def test_splits_are_valid(self) -> None:
        """Test that all returned splits are valid (non-empty sets)."""
        from tests.fixtures.directed_networks import LEVEL_1_DNETWORK_SINGLE_HYBRID

        splits = induced_splits(LEVEL_1_DNETWORK_SINGLE_HYBRID)

        for split in splits:
            assert len(split.set1) > 0
            assert len(split.set2) > 0
            assert split.set1.isdisjoint(split.set2)

    def test_tree_of_blobs_same_splits_multiple_blobs_level2(self) -> None:
        """Test that tree-of-blobs preserves splits for a level 2 network with multiple blobs."""
        from tests.fixtures.directed_networks import LEVEL_2_DNETWORK_MULTIPLE_BLOBS

        original_splits = induced_splits(LEVEL_2_DNETWORK_MULTIPLE_BLOBS)
        blob_tree = tree_of_blobs(LEVEL_2_DNETWORK_MULTIPLE_BLOBS)
        blob_tree_splits = induced_splits(blob_tree)

        # Split systems should be equal
        assert original_splits.splits == blob_tree_splits.splits


class TestDisplayedSplits:
    """Test displayed_splits function for DirectedPhyNetwork."""

    def test_simple_tree(self) -> None:
        """Test displayed_splits on a simple tree."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})]
        )
        splits = displayed_splits(net)
        from phylozoo.core.split import WeightedSplitSystem

        assert isinstance(splits, WeightedSplitSystem)
        # Tree has one displayed tree (itself) with probability 1.0
        assert len(splits) >= 0  # May have splits or be empty for small trees

    def test_single_hybrid_network(self) -> None:
        """Test displayed_splits on a network with a single hybrid."""
        from tests.fixtures.directed_networks import LEVEL_1_DNETWORK_SINGLE_HYBRID

        splits = displayed_splits(LEVEL_1_DNETWORK_SINGLE_HYBRID)
        from phylozoo.core.split import WeightedSplitSystem

        assert isinstance(splits, WeightedSplitSystem)
        # Should have splits from displayed trees
        assert len(splits) > 0

    def test_weights_sum_to_one(self) -> None:
        """Test that weights in displayed_splits sum appropriately."""
        net = DirectedPhyNetwork(
            edges=[(10, 5), (10, 6), (5, 4), (6, 4), (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (7, {"label": "D"}),
            ],
        )
        _ = displayed_splits(net)
        # Get all displayed trees to verify probabilities
        trees = list(displayed_trees(net, probability=True))
        total_prob = sum(tree.get_network_attribute("probability") or 1.0 for tree in trees)
        # Total probability should be 1.0 (or close due to floating point)
        assert abs(total_prob - 1.0) < 1e-10

    def test_splits_from_displayed_trees(self) -> None:
        """Test that displayed_splits contains splits from all displayed trees."""
        net = DirectedPhyNetwork(
            edges=[(10, 5), (10, 6), (5, 4), (6, 4), (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (7, {"label": "D"}),
            ],
        )
        splits = displayed_splits(net)

        # Get all displayed trees and their splits
        all_tree_splits = set()
        for tree in displayed_trees(net, probability=True):
            tree_splits = induced_splits(tree)
            all_tree_splits.update(tree_splits.splits)

        # All splits from displayed trees should be in displayed_splits
        assert splits.splits.issuperset(all_tree_splits)

    def test_weights_accumulate(self) -> None:
        """Test that weights accumulate when splits appear in multiple trees."""
        net = DirectedPhyNetwork(
            edges=[(10, 5), (10, 6), (5, 4), (6, 4), (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (7, {"label": "D"}),
            ],
        )
        splits = displayed_splits(net)

        # Check that if a split appears in multiple trees, its weight is the sum
        # This is verified by checking that weights are positive and reasonable
        for split in splits.splits:
            weight = splits.get_weight(split)
            assert weight > 0
            assert weight <= 1.0  # Should not exceed 1.0

    def test_empty_network(self) -> None:
        """Test displayed_splits on empty network."""
        net = DirectedPhyNetwork()
        splits = displayed_splits(net)
        from phylozoo.core.split import WeightedSplitSystem

        assert isinstance(splits, WeightedSplitSystem)
        assert len(splits) == 0


class TestSplitFromCutedge:
    """Test split_from_cutedge function for DirectedPhyNetwork."""

    def test_simple_tree_cutedge(self) -> None:
        """Test split_from_cutedge on a simple tree."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (4, {"label": "C"})],
        )
        split = split_from_cutedge(net, 3, 1)
        assert isinstance(split, Split)
        # The split should separate A from B and C
        assert "A" in split.set1 or "A" in split.set2
        assert {"B", "C"} == (split.set1 | split.set2) - {"A"}

    def test_cutedge_with_key(self) -> None:
        """Test split_from_cutedge with explicit key."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})]
        )
        split = split_from_cutedge(net, 3, 1, key=0)
        assert isinstance(split, Split)
        assert "A" in split.set1 or "A" in split.set2
        assert "B" in split.set1 or "B" in split.set2

    def test_nonexistent_edge(self) -> None:
        """Test split_from_cutedge raises error for nonexistent edge."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})]
        )
        with pytest.raises(ValueError, match="does not exist"):
            split_from_cutedge(net, 1, 2)

    def test_non_cutedge(self) -> None:
        """Test split_from_cutedge raises error for non-cut-edge."""
        # Use a fixture network with a blob (non-cut-edges exist)
        from tests.fixtures.directed_networks import LEVEL_1_DNETWORK_SINGLE_HYBRID

        net = LEVEL_1_DNETWORK_SINGLE_HYBRID
        # Find a non-cut-edge (edge in a blob)
        from phylozoo.core.network.dnetwork.features import cut_edges

        cut_edges_set = cut_edges(net)
        all_edges = set()
        for u, v, key in net._graph.edges(keys=True):
            all_edges.add((u, v, key))
        # Find an edge that's not a cut-edge
        non_cut_edges = all_edges - cut_edges_set
        if not non_cut_edges:
            pytest.skip("No non-cut-edges found in the test network")
        u, v, key = next(iter(non_cut_edges))
        with pytest.raises(ValueError, match="is not a cut-edge"):
            split_from_cutedge(net, u, v, key=key)

    def test_split_covers_all_taxa(self) -> None:
        """Test that split covers all taxa in the network."""
        net = DirectedPhyNetwork(
            edges=[(5, 1), (5, 2), (5, 3), (5, 4)],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (4, {"label": "D"}),
            ],
        )
        split = split_from_cutedge(net, 5, 1)
        all_taxa = split.set1 | split.set2
        assert all_taxa == {"A", "B", "C", "D"}

    def test_split_set1_on_u_side(self) -> None:
        """Test that set1 contains taxa on the side of u."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})]
        )
        split = split_from_cutedge(net, 3, 1)
        # A should be on the side of node 1
        assert "A" in split.set1 or "A" in split.set2
        # B should be on the side of node 3
        assert "B" in split.set1 or "B" in split.set2
        # They should be in different sets
        assert ("A" in split.set1 and "B" in split.set2) or (
            "A" in split.set2 and "B" in split.set1
        )

    def test_return_node_taxa(self) -> None:
        """Test split_from_cutedge with return_node_taxa parameter."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (4, {"label": "C"})],
        )
        result = split_from_cutedge(net, 3, 1, return_node_taxa=True)
        split, (u_node, u_taxa), (v_node, v_taxa) = result
        assert isinstance(split, Split)
        assert u_node == 3
        assert v_node == 1
        assert isinstance(u_taxa, frozenset)
        assert isinstance(v_taxa, frozenset)
        # Check that taxa match the split
        assert u_taxa == split.set1 or u_taxa == split.set2
        assert v_taxa == split.set1 or v_taxa == split.set2
        # Check that u_taxa and v_taxa are complementary
        assert u_taxa | v_taxa == split.set1 | split.set2
        assert u_taxa & v_taxa == set()

    def test_hybrid_edge_cutedge(self) -> None:
        """Test split_from_cutedge on a hybrid edge."""
        # Use a fixture network that's known to be valid
        from tests.fixtures.directed_networks import LEVEL_1_DNETWORK_SINGLE_HYBRID

        net = LEVEL_1_DNETWORK_SINGLE_HYBRID
        # Find a cut-edge in the network
        from phylozoo.core.network.dnetwork.features import cut_edges

        cut_edges_set = cut_edges(net)
        if cut_edges_set:
            u, v, key = next(iter(cut_edges_set))
            split = split_from_cutedge(net, u, v, key=key)
            assert isinstance(split, Split)
            assert len(split.set1) > 0
            assert len(split.set2) > 0


class TestDisplayedQuartets:
    """Test displayed_quartets function for DirectedPhyNetwork."""

    def test_simple_tree_four_taxa(self) -> None:
        """Test displayed_quartets on a simple tree with exactly 4 taxa."""
        net = DirectedPhyNetwork(
            edges=[(5, 1), (5, 2), (5, 6), (6, 3), (6, 4)],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (4, {"label": "D"}),
            ],
        )
        profileset = displayed_quartets(net)

        # Should have exactly one profile (one 4-taxon set)
        assert len(profileset) == 1
        assert profileset.taxa == frozenset({"A", "B", "C", "D"})

        # Get the profile
        profile, profile_weight = profileset.profiles[frozenset({"A", "B", "C", "D"})]
        # Profile should have default weight 1.0
        assert abs(profile_weight - 1.0) < 1e-10
        # Should have exactly one quartet (tree has one displayed tree)
        assert len(profile) == 1

        # The quartet should be resolved (not a star tree)
        quartet = next(iter(profile.quartets.keys()))
        assert quartet.is_resolved()
        assert not quartet.is_star()
        # Weight should be 1.0 (single displayed tree with probability 1.0)
        assert abs(profile.get_weight(quartet) - 1.0) < 1e-10

    def test_network_with_hybrid_four_taxa(self) -> None:
        """Test displayed_quartets on network with hybrid node and 4 taxa."""
        from tests.fixtures.directed_networks import LEVEL_1_DNETWORK_SINGLE_HYBRID

        net = LEVEL_1_DNETWORK_SINGLE_HYBRID

        # Check if network has at least 4 taxa
        if len(net.taxa) >= 4:
            profileset = displayed_quartets(net)

            # Should have at least one profile
            assert len(profileset) >= 1
            assert profileset.taxa == net.taxa

            # Get the profile
            profile, profile_weight = next(iter(profileset.profiles.values()))
            # Profile should have default weight 1.0
            assert abs(profile_weight - 1.0) < 1e-10

            # Network has multiple displayed trees, weights should sum to 1.0
            total_quartet_weight = sum(profile.quartets.values())
            assert abs(total_quartet_weight - 1.0) < 1e-10

    def test_network_with_hybrid_gamma_values(self) -> None:
        """Test displayed_quartets with explicit gamma values."""
        net = DirectedPhyNetwork(
            edges=[
                {"u": 5, "v": 4, "gamma": 0.6},
                {"u": 6, "v": 4, "gamma": 0.4},
                (4, 8),
                (8, 1),
                (8, 2),
                (5, 3),
                (5, 6),
                (6, 7),
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (7, {"label": "D"}),
            ],
        )
        profileset = displayed_quartets(net)

        # Should have exactly one profile
        assert len(profileset) == 1

        # Get the profile
        profile, _ = next(iter(profileset.profiles.values()))

        # Network has 2 displayed trees with probabilities 0.6 and 0.4
        # Total weight should be 1.0
        total_quartet_weight = sum(profile.quartets.values())
        assert abs(total_quartet_weight - 1.0) < 1e-10

        # Check that individual weights match probabilities
        for quartet, weight in profile.quartets.items():
            assert weight > 0.0
            assert weight <= 1.0

    def test_network_fewer_than_four_taxa(self) -> None:
        """Test displayed_quartets on network with fewer than 4 taxa."""
        from tests.fixtures.directed_networks import DTREE_SMALL_BINARY

        net = DTREE_SMALL_BINARY

        profileset = displayed_quartets(net)

        # Should return empty QuartetProfileSet
        assert len(profileset) == 0
        assert len(profileset.taxa) == 0

    def test_empty_network(self) -> None:
        """Test displayed_quartets on empty network."""
        from tests.fixtures.directed_networks import DTREE_EMPTY

        net = DTREE_EMPTY

        profileset = displayed_quartets(net)

        # Should return empty QuartetProfileSet
        assert len(profileset) == 0
        assert len(profileset.taxa) == 0

    def test_network_more_than_four_taxa(self) -> None:
        """Test displayed_quartets on network with more than 4 taxa."""
        net = DirectedPhyNetwork(
            edges=[(7, 1), (7, 2), (7, 8), (8, 3), (8, 4), (8, 5), (8, 6)],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (4, {"label": "D"}),
                (5, {"label": "E"}),
                (6, {"label": "F"}),
            ],
        )
        profileset = displayed_quartets(net)

        # Should have multiple profiles (one for each 4-taxon combination)
        # For 6 taxa, there are C(6,4) = 15 combinations
        assert len(profileset) == 15
        assert profileset.taxa == net.taxa

        # Each profile should have default weight 1.0
        for profile, profile_weight in profileset.profiles.values():
            assert abs(profile_weight - 1.0) < 1e-10
            # Each profile should have at least one quartet
            assert len(profile) >= 1
            # All quartets in a profile should have the same 4 taxa
            for quartet in profile.quartets.keys():
                assert quartet.taxa in profileset.profiles

    def test_weights_sum_correctly(self) -> None:
        """Test that weights are correctly summed when same quartet appears multiple times."""
        net = DirectedPhyNetwork(
            edges=[
                {"u": 5, "v": 4, "gamma": 0.6},
                {"u": 6, "v": 4, "gamma": 0.4},
                (4, 8),
                (8, 1),
                (8, 2),
                (5, 3),
                (5, 6),
                (6, 7),
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (7, {"label": "D"}),
            ],
        )
        profileset = displayed_quartets(net)

        # Get the profile
        profile, _ = next(iter(profileset.profiles.values()))

        # Total weight of all quartets should sum to 1.0
        # (sum of probabilities of all displayed trees)
        total_weight = sum(profile.quartets.values())
        assert abs(total_weight - 1.0) < 1e-10

    def test_profile_weights_default(self) -> None:
        """Test that profile weights are default (1.0) as specified."""
        net = DirectedPhyNetwork(
            edges=[(5, 1), (5, 2), (5, 6), (6, 3), (6, 4)],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (4, {"label": "D"}),
            ],
        )
        profileset = displayed_quartets(net)

        # All profiles should have default weight 1.0
        for profile, profile_weight in profileset.profiles.values():
            assert abs(profile_weight - 1.0) < 1e-10

    def test_uses_sd_network_conversion(self) -> None:
        """Test that displayed_quartets uses to_sd_network conversion (unrooted quartets)."""
        # Create a directed network with a root
        net = DirectedPhyNetwork(
            edges=[(5, 1), (5, 2), (5, 6), (6, 3), (6, 4)],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (4, {"label": "D"}),
            ],
        )
        profileset = displayed_quartets(net)

        # Should have one profile
        assert len(profileset) == 1

        # Get the quartet
        profile, _ = next(iter(profileset.profiles.values()))
        quartet = next(iter(profile.quartets.keys()))

        # The quartet should be unrooted (not rooted)
        # This is verified by the fact that it's a resolved quartet with a 2|2 split
        assert quartet.is_resolved()
        assert quartet.split is not None
        # The split should be a 2|2 split (unrooted quartet)
        assert len(quartet.split.set1) == 2
        assert len(quartet.split.set2) == 2


class TestDisplayedTriplets:
    """Test displayed_triplets function for DirectedPhyNetwork."""

    def test_simple_tree_three_taxa(self) -> None:
        """Test displayed_triplets on a simple binary rooted tree with exactly 3 taxa."""
        net = DirectedPhyNetwork(
            edges=[(4, 1), (4, 5), (5, 2), (5, 3)],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (3, {"label": "C"})],
        )
        profileset = displayed_triplets(net)

        # Should have exactly one profile (one 3-taxon set)
        assert len(profileset) == 1
        assert profileset.taxa == frozenset({"A", "B", "C"})

        profile, profile_weight = profileset.profiles[frozenset({"A", "B", "C"})]
        assert abs(profile_weight - 1.0) < 1e-10
        # Single displayed tree, so single triplet
        assert len(profile) == 1

        # The triplet should be resolved with A as outgroup
        triplet = next(iter(profile.triplets.keys()))
        assert triplet.is_resolved()
        assert triplet.outgroup == frozenset({"A"})
        assert triplet.cherry == frozenset({"B", "C"})
        # Weight should be 1.0 (single displayed tree with probability 1.0)
        assert abs(profile.get_weight(triplet) - 1.0) < 1e-10

    def test_star_tree_three_taxa(self) -> None:
        """Test displayed_triplets on a rooted star tree with 3 taxa."""
        net = DirectedPhyNetwork(
            edges=[(4, 1), (4, 2), (4, 3)],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (3, {"label": "C"})],
        )
        profileset = displayed_triplets(net)

        assert len(profileset) == 1
        profile, _ = profileset.profiles[frozenset({"A", "B", "C"})]
        assert len(profile) == 1
        triplet = next(iter(profile.triplets.keys()))
        assert triplet.is_star()

    def test_network_with_hybrid_four_taxa(self) -> None:
        """Test displayed_triplets on network with hybrid node and 4 taxa."""
        net = DirectedPhyNetwork(
            edges=[
                {"u": 5, "v": 4, "gamma": 0.6},
                {"u": 6, "v": 4, "gamma": 0.4},
                (4, 8),
                (8, 1),
                (8, 2),
                (5, 3),
                (5, 6),
                (6, 7),
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (7, {"label": "D"}),
            ],
        )
        profileset = displayed_triplets(net)

        # 4 taxa: C(4,3) = 4 profiles
        assert len(profileset) == 4
        assert profileset.taxa == net.taxa

        # Each profile's triplet weights sum to 1.0 (displayed-tree probabilities)
        for profile, profile_weight in profileset.profiles.values():
            assert abs(profile_weight - 1.0) < 1e-10
            total_triplet_weight = sum(profile.triplets.values())
            assert abs(total_triplet_weight - 1.0) < 1e-10

    def test_network_with_hybrid_gamma_weights_split(self) -> None:
        """Test displayed_triplets distributes gamma weights across triplets for an ambiguous 3-taxon set."""
        net = DirectedPhyNetwork(
            edges=[
                {"u": 5, "v": 4, "gamma": 0.6},
                {"u": 6, "v": 4, "gamma": 0.4},
                (4, 8),
                (8, 1),
                (8, 2),
                (5, 3),
                (5, 6),
                (6, 7),
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (7, {"label": "D"}),
            ],
        )
        profileset = displayed_triplets(net)

        # On {A, C, D} the two displayed trees produce different rooted topologies,
        # so the profile should contain two triplets with weights 0.6 and 0.4.
        acd_profile = profileset.get_profile(frozenset({"A", "C", "D"}))
        assert acd_profile is not None
        assert len(acd_profile) == 2
        weights = sorted(acd_profile.triplets.values())
        assert weights == pytest.approx([0.4, 0.6])

    def test_network_fewer_than_three_taxa(self) -> None:
        """Test displayed_triplets on network with fewer than 3 taxa."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)], nodes=[(1, {"label": "A"}), (2, {"label": "B"})]
        )
        profileset = displayed_triplets(net)

        # Should return empty TripletProfileSet
        assert len(profileset) == 0
        assert len(profileset.taxa) == 0

    def test_empty_network(self) -> None:
        """Test displayed_triplets on empty network."""
        from tests.fixtures.directed_networks import DTREE_EMPTY

        net = DTREE_EMPTY

        profileset = displayed_triplets(net)

        # Should return empty TripletProfileSet
        assert len(profileset) == 0
        assert len(profileset.taxa) == 0

    def test_network_more_than_three_taxa(self) -> None:
        """Test displayed_triplets on a tree with more than 3 taxa."""
        net = DirectedPhyNetwork(
            edges=[(7, 1), (7, 8), (8, 2), (8, 9), (9, 3), (9, 4)],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (4, {"label": "D"}),
            ],
        )
        profileset = displayed_triplets(net)

        # For 4 taxa, there are C(4,3) = 4 combinations
        assert len(profileset) == 4
        assert profileset.taxa == net.taxa

        # Each profile should have default weight 1.0 and contain at least one triplet
        for profile, profile_weight in profileset.profiles.values():
            assert abs(profile_weight - 1.0) < 1e-10
            assert len(profile) >= 1
            for triplet in profile.triplets.keys():
                assert triplet.taxa in profileset.profiles

    def test_weights_sum_correctly(self) -> None:
        """Test that weights are correctly summed when same triplet appears multiple times."""
        net = DirectedPhyNetwork(
            edges=[
                {"u": 5, "v": 4, "gamma": 0.6},
                {"u": 6, "v": 4, "gamma": 0.4},
                (4, 8),
                (8, 1),
                (8, 2),
                (5, 3),
                (5, 6),
                (6, 7),
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (7, {"label": "D"}),
            ],
        )
        profileset = displayed_triplets(net)

        # For each profile, total weight should sum to 1.0 across all displayed trees
        for profile, _ in profileset.profiles.values():
            total_weight = sum(profile.triplets.values())
            assert abs(total_weight - 1.0) < 1e-10

    def test_profile_weights_default(self) -> None:
        """Test that profile weights are default (1.0) as specified."""
        net = DirectedPhyNetwork(
            edges=[(4, 1), (4, 5), (5, 2), (5, 3)],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (3, {"label": "C"})],
        )
        profileset = displayed_triplets(net)

        # All profiles should have default weight 1.0
        for profile, profile_weight in profileset.profiles.values():
            assert abs(profile_weight - 1.0) < 1e-10

    def test_rooted_output(self) -> None:
        """Test that displayed_triplets returns ROOTED triplets (with 1|2 trivial splits)."""
        net = DirectedPhyNetwork(
            edges=[(4, 1), (4, 5), (5, 2), (5, 3)],
            nodes=[(1, {"label": "A"}), (2, {"label": "B"}), (3, {"label": "C"})],
        )
        profileset = displayed_triplets(net)

        profile, _ = next(iter(profileset.profiles.values()))
        triplet = next(iter(profile.triplets.keys()))

        # The triplet should be resolved with a trivial 1|2 split
        assert triplet.is_resolved()
        assert triplet.split is not None
        assert triplet.split.is_trivial
        # One side is a singleton (outgroup), other has 2 elements (cherry)
        sizes = sorted((len(triplet.split.set1), len(triplet.split.set2)))
        assert sizes == [1, 2]


class TestPartitionFromBlob:
    """Test partition_from_blob function for DirectedPhyNetwork."""

    def test_basic_partition(self) -> None:
        """Test basic partition_from_blob functionality."""
        # Network with hybrid node creating non-leaf blob
        net = DirectedPhyNetwork(
            edges=[
                (8, 5),
                (8, 6),
                (5, 1),
                (5, 2),
                (6, 3),
                (6, 9),
                (5, 4),
                (6, 4),  # Hybrid node 4
                (4, 7),
                (7, 10),
                (7, 11),
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (9, {"label": "D"}),
                (10, {"label": "E"}),
                (11, {"label": "F"}),
            ],
        )
        # Non-leaf blob is {4, 5, 6, 8}
        partition = partition_from_blob(net, {4, 5, 6, 8})

        assert isinstance(partition, Partition)
        assert len(partition) == 5
        # Leaves E and F are in the same component (both connected via node 7)
        assert {"A"} in partition
        assert {"B"} in partition
        assert {"C"} in partition
        assert {"D"} in partition
        assert {"E", "F"} in partition

    def test_partition_with_edge_taxa(self) -> None:
        """Test partition_from_blob with return_edge_taxa=True."""
        # Network with hybrid node creating non-leaf blob
        net = DirectedPhyNetwork(
            edges=[
                (8, 5),
                (8, 6),
                (5, 1),
                (5, 2),
                (6, 3),
                (6, 9),
                (5, 4),
                (6, 4),  # Hybrid node 4
                (4, 7),
                (7, 10),
                (7, 11),
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (9, {"label": "D"}),
                (10, {"label": "E"}),
                (11, {"label": "F"}),
            ],
        )
        # Non-leaf blob is {4, 5, 6, 8}
        partition, edge_taxa = partition_from_blob(net, {4, 5, 6, 8}, return_edge_taxa=True)

        assert isinstance(partition, Partition)
        assert len(partition) == 5
        assert isinstance(edge_taxa, list)
        assert len(edge_taxa) == 5

        # Check that each tuple has the correct format
        for u, v, taxa_set in edge_taxa:
            assert isinstance(u, (int, str))
            assert isinstance(v, (int, str))
            assert isinstance(taxa_set, frozenset)
            assert v in {4, 5, 6, 8}  # All should connect to blob nodes

        # Check that all taxa are covered
        all_taxa = {taxon for _, _, taxa_set in edge_taxa for taxon in taxa_set}
        assert all_taxa == {"A", "B", "C", "D", "E", "F"}

    def test_empty_blob_raises_error(self) -> None:
        """Test that empty blob raises ValueError."""
        net = DirectedPhyNetwork(
            edges=[
                (8, 5),
                (8, 6),
                (5, 1),
                (5, 2),
                (6, 3),
                (6, 9),
                (5, 4),
                (6, 4),
                (4, 7),
                (7, 10),
                (7, 11),
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (9, {"label": "D"}),
                (10, {"label": "E"}),
                (11, {"label": "F"}),
            ],
        )
        with pytest.raises(ValueError, match="Blob cannot be empty"):
            partition_from_blob(net, set())

    def test_blob_not_in_network_raises_error(self) -> None:
        """Test that blob with nodes not in network raises ValueError."""
        net = DirectedPhyNetwork(
            edges=[
                (8, 5),
                (8, 6),
                (5, 1),
                (5, 2),
                (6, 3),
                (6, 9),
                (5, 4),
                (6, 4),
                (4, 7),
                (7, 10),
                (7, 11),
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (9, {"label": "D"}),
                (10, {"label": "E"}),
                (11, {"label": "F"}),
            ],
        )
        with pytest.raises(ValueError, match="Blob contains nodes not in network"):
            partition_from_blob(net, {99})  # Node 99 not in network

    def test_leaf_blob_raises_error(self) -> None:
        """Test that a leaf blob (single leaf node) raises ValueError."""
        net = DirectedPhyNetwork(
            edges=[
                (8, 5),
                (8, 6),
                (5, 1),
                (5, 2),
                (6, 3),
                (6, 9),
                (5, 4),
                (6, 4),
                (4, 7),
                (7, 10),
                (7, 11),
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (9, {"label": "D"}),
                (10, {"label": "E"}),
                (11, {"label": "F"}),
            ],
        )
        with pytest.raises(ValueError, match="is not a non-leaf blob"):
            partition_from_blob(net, {1})  # Leaf node blob

    def test_multiple_node_blob(self) -> None:
        """Test partition_from_blob with a blob containing multiple nodes."""
        # Network with hybrid node creating non-leaf blob
        net = DirectedPhyNetwork(
            edges=[
                (8, 5),
                (8, 6),
                (5, 1),
                (5, 2),
                (6, 3),
                (6, 9),
                (5, 4),
                (6, 4),  # Hybrid node 4
                (4, 7),
                (7, 10),
                (7, 11),
            ],
            nodes=[
                (1, {"label": "A"}),
                (2, {"label": "B"}),
                (3, {"label": "C"}),
                (9, {"label": "D"}),
                (10, {"label": "E"}),
                (11, {"label": "F"}),
            ],
        )
        # Non-leaf blob is {4, 5, 6, 8}
        partition = partition_from_blob(net, {4, 5, 6, 8})

        assert isinstance(partition, Partition)
        assert len(partition) == 5
        # Leaves E and F are in the same component (both connected via node 7)
        assert {"A"} in partition
        assert {"B"} in partition
        assert {"C"} in partition
        assert {"D"} in partition
        assert {"E", "F"} in partition


class TestPruneDegree1Nodes:
    """Directly exercise the worklist prune used when building displayed trees."""

    def _graph(self, edges):
        from phylozoo.core.primitives.d_multigraph.base import DirectedMultiGraph

        graph = DirectedMultiGraph()
        for u, v in edges:
            graph.add_edge(u, v)
        return graph

    def test_removes_chain_exposed_by_earlier_removal(self):
        """A whole pendant chain goes, not just its tip: removal exposes new degree-1 nodes."""
        from phylozoo.core.network.dnetwork._utils import _prune_degree1_nodes

        graph = self._graph([("r", "a"), ("a", "L1"), ("r", "b"), ("b", "c"), ("c", "d")])
        _prune_degree1_nodes(graph, {"r", "L1"})
        assert set(graph.nodes()) == {"r", "a", "L1"}

    def test_keeps_protected_nodes_however_low_their_degree(self):
        from phylozoo.core.network.dnetwork._utils import _prune_degree1_nodes

        graph = self._graph([("r", "x"), ("x", "L1")])
        _prune_degree1_nodes(graph, {"r", "L1"})
        assert set(graph.nodes()) == {"r", "x", "L1"}

    def test_two_adjacent_degree1_nodes_leave_one_behind(self):
        """Removing one of a connected pair drops the other to degree 0, which is not pruned."""
        from phylozoo.core.network.dnetwork._utils import _prune_degree1_nodes

        graph = self._graph([("a", "b")])
        _prune_degree1_nodes(graph, set())
        assert len(list(graph.nodes())) == 1

    def test_no_degree1_nodes_is_a_noop(self):
        from phylozoo.core.network.dnetwork._utils import _prune_degree1_nodes

        graph = self._graph([("r", "L1"), ("r", "L2")])
        _prune_degree1_nodes(graph, {"r", "L1", "L2"})
        assert set(graph.nodes()) == {"r", "L1", "L2"}


class TestDistancesBlobDecomposition:
    """`distances` enumerates each blob's switchings separately; results must not change.

    The reference here is the definition itself: aggregate over every global switching.
    These networks are built so their reticulations land in *different* blobs, which is
    exactly the case the decomposition changes.
    """

    @staticmethod
    def _switching_matrix_nx(graph, taxa, network):
        """Leaf-to-leaf distances of one switching graph via networkx, independent of the code under test."""
        import networkx as nx
        import numpy as np

        weighted = nx.Graph()
        for u, v, key in graph._combined.edges(keys=True):
            length = network.get_branch_length(u, v, key)
            if length is None:
                length = network.get_branch_length(v, u, key)
            weighted.add_edge(u, v, weight=1.0 if length is None else length)
        leaves = [network._label_to_node[taxon] for taxon in taxa]
        matrix = np.zeros((len(taxa), len(taxa)))
        for i, source in enumerate(leaves):
            lengths = nx.shortest_path_length(weighted, source, weight="weight")
            for j, target in enumerate(leaves):
                matrix[i, j] = lengths[target]
        return matrix

    @staticmethod
    def _brute_force(network, taxa, mode):
        """Aggregate over all global switchings, i.e. the unoptimised definition."""
        import numpy as np

        from phylozoo.core.network.dnetwork.derivations import _switchings

        size = len(taxa)
        result = np.full((size, size), np.inf) if mode == "shortest" else np.zeros((size, size))
        weighted, weight_total = np.zeros((size, size)), 0.0
        for graph in _switchings(network, probability=(mode == "average")):
            matrix = TestDistancesBlobDecomposition._switching_matrix_nx(graph, taxa, network)
            if mode == "shortest":
                result = np.minimum(result, matrix)
            elif mode == "longest":
                result = np.maximum(result, matrix)
            else:
                weight = graph._graph.graph.get("probability", 1.0) or 1.0
                weighted += matrix * weight
                weight_total += weight
        if mode == "average":
            result = weighted / weight_total
        np.fill_diagonal(result, 0.0)
        return result

    @staticmethod
    def _two_blob_network(gamma=False):
        """Two level-1 blobs on separate sides of the root, so they switch independently."""
        edges = [
            ("rho", "a0"),
            ("rho", "b0"),
            # left blob: hybrid ha
            ("a0", "a1"),
            ("a0", "a2"),
            ("a1", "ha"),
            ("a2", "ha"),
            ("a1", "A"),
            ("a2", "B"),
            ("ha", "C"),
            # right blob: hybrid hb
            ("b0", "b1"),
            ("b0", "b2"),
            ("b1", "hb"),
            ("b2", "hb"),
            ("b1", "D"),
            ("b2", "E"),
            ("hb", "F"),
        ]
        if gamma:
            edges = [
                (
                    {"u": u, "v": v, "gamma": 0.25}
                    if v in ("ha", "hb") and u in ("a1", "b1")
                    else ({"u": u, "v": v, "gamma": 0.75} if v in ("ha", "hb") else (u, v))
                )
                for u, v in edges
            ]
        labels = ["A", "B", "C", "D", "E", "F"]
        return DirectedPhyNetwork(edges=edges, nodes=[(name, {"label": name}) for name in labels])

    def test_network_really_has_two_blobs(self):
        """Guard the premise: if this became single-blob the tests below would be vacuous."""
        from phylozoo.core.network.dnetwork.derivations import _hybrid_blob_groups

        network = self._two_blob_network()
        groups = _hybrid_blob_groups(network)
        assert len(groups) == 2
        assert sorted(len(group) for group in groups) == [1, 1]

    @pytest.mark.parametrize("mode", ["shortest", "longest", "average"])
    def test_matches_full_enumeration_across_blobs(self, mode):
        import numpy as np

        network = self._two_blob_network()
        taxa = sorted(network.taxa)
        got = np.asarray(distances(network, mode=mode)._matrix)
        assert np.allclose(got, self._brute_force(network, taxa, mode))

    @pytest.mark.parametrize("mode", ["shortest", "longest", "average"])
    def test_matches_full_enumeration_with_gamma(self, mode):
        """Weighted averaging must normalise per blob exactly as it did globally."""
        import numpy as np

        network = self._two_blob_network(gamma=True)
        taxa = sorted(network.taxa)
        got = np.asarray(distances(network, mode=mode)._matrix)
        assert np.allclose(got, self._brute_force(network, taxa, mode))

    @pytest.mark.parametrize("mode", ["shortest", "longest", "average"])
    def test_single_blob_still_matches(self, mode):
        """Two hybrids in one blob: the decomposition must not split what interacts."""
        import numpy as np

        from phylozoo.core.network.dnetwork.derivations import _hybrid_blob_groups

        network = DirectedPhyNetwork(
            edges=[
                ("rho", "u"),
                ("rho", "w"),
                ("u", "h1"),
                ("w", "h1"),
                ("u", "h2"),
                ("w", "h2"),
                ("u", "A"),
                ("w", "B"),
                ("h1", "C"),
                ("h2", "D"),
            ],
            nodes=[(name, {"label": name}) for name in ("A", "B", "C", "D")],
        )
        assert len(_hybrid_blob_groups(network)) == 1  # they do interact
        taxa = sorted(network.taxa)
        got = np.asarray(distances(network, mode=mode)._matrix)
        assert np.allclose(got, self._brute_force(network, taxa, mode))

    @pytest.mark.parametrize("mode", ["shortest", "longest", "average"])
    def test_tree_has_no_hybrids_to_group(self, mode):
        """A tree has one switching; all three modes must agree on it."""
        import numpy as np

        from phylozoo.core.network.dnetwork.derivations import _hybrid_blob_groups

        network = DirectedPhyNetwork(
            edges=[("r", "x"), ("r", "C"), ("x", "A"), ("x", "B")],
            nodes=[(name, {"label": name}) for name in ("A", "B", "C")],
        )
        assert _hybrid_blob_groups(network) == []
        taxa = sorted(network.taxa)
        got = np.asarray(distances(network, mode=mode)._matrix)
        assert np.allclose(got, self._brute_force(network, taxa, mode))

    def test_shortest_at_most_average_at_most_longest(self):
        """A sanity relation the three modes must satisfy on the same network."""
        import numpy as np

        network = self._two_blob_network()
        shortest = np.asarray(distances(network, mode="shortest")._matrix)
        average = np.asarray(distances(network, mode="average")._matrix)
        longest = np.asarray(distances(network, mode="longest")._matrix)
        assert np.all(shortest <= average + 1e-12)
        assert np.all(average <= longest + 1e-12)

    def test_blobs_are_enumerated_separately_not_jointly(self):
        """With k independent blobs the work is sum(2**r_b), not prod(2**r_b)."""
        from phylozoo.core.network.dnetwork import derivations

        network = self._two_blob_network()
        calls = []
        original = derivations._switching_distance_matrix

        def counting(*args, **kwargs):
            calls.append(1)
            return original(*args, **kwargs)

        derivations._switching_distance_matrix = counting
        try:
            distances(network, mode="average")
        finally:
            derivations._switching_distance_matrix = original
        # 1 reference + 2 local switchings per blob x 2 blobs; never the 4 global ones
        assert len(calls) == 5


class TestDistancesWithParallelEdges:
    """Branch lengths must be looked up by edge key when parallel edges exist.

    An unkeyed lookup raises ``PhyloZooValueError: Multiple parallel edges exist``,
    even though a switching keeps only one of the parallel edges and the network
    constrains them all to share a branch length.
    """

    @staticmethod
    def _network():
        """Hybrid ``h`` reached by two parallel edges from ``u``."""
        return DirectedPhyNetwork(
            edges=[
                ("r", "u"),
                ("r", "B"),
                {"u": "u", "v": "h", "key": 0, "branch_length": 2.0},
                {"u": "u", "v": "h", "key": 1, "branch_length": 2.0},
                ("h", "C"),
                ("u", "A"),
            ],
            nodes=[(name, {"label": name}) for name in ("A", "B", "C")],
        )

    def test_network_really_has_parallel_edges(self):
        """Guard the premise, so this cannot quietly stop testing what it claims to."""
        network = self._network()
        parents = list(network.incident_parent_edges("h", keys=True))
        assert len(parents) == 2
        assert {edge[0] for edge in parents} == {"u"}  # both from the same node

    @pytest.mark.parametrize("mode", ["shortest", "longest", "average"])
    def test_distances_does_not_raise_on_parallel_edges(self, mode):
        distances(self._network(), mode=mode)

    @pytest.mark.parametrize("mode", ["shortest", "longest", "average"])
    def test_parallel_edge_branch_length_is_used(self, mode):
        """A -> u -> h -> C is 1.0 + 2.0 + 1.0; the 2.0 comes from the parallel edge."""
        matrix = distances(self._network(), mode=mode)
        assert matrix.get_distance("A", "C") == pytest.approx(4.0)
        assert matrix.get_distance("A", "B") == pytest.approx(3.0)


class TestDisplayedSplitsBlobDecomposition:
    """`displayed_splits` sums per blob instead of over displayed trees; weights must match.

    The reference is the definition: accumulate each displayed tree's induced splits
    with the tree's probability. The networks put reticulations in *different* blobs,
    which is exactly the case the decomposition changes.
    """

    @staticmethod
    def _by_definition(network):
        weights = {}
        for tree in displayed_trees(network, probability=True):
            probability = tree.get_network_attribute("probability")
            probability = 1.0 if probability is None else probability
            for split in induced_splits(tree).splits:
                weights[split] = weights.get(split, 0.0) + probability
        return weights

    @staticmethod
    def _two_blob_network(gamma=False):
        """Two cherries each carrying one reticulation: two level-1 blobs.

        Switching a reticulation off leaves a path of degree-2 nodes whose edges repeat
        the adjacent cut edge's split, so this also exercises the de-duplication.
        """
        edges = [
            ("rho", "c1"),
            ("rho", "c2"),
            ("c1", "m1"),
            ("m1", "A"),
            ("c1", "m2"),
            ("m2", "B"),
            ("m1", "m2"),
            ("c2", "m3"),
            ("m3", "C"),
            ("c2", "m4"),
            ("m4", "D"),
            ("m3", "m4"),
        ]
        if gamma:
            weighted = {("m1", "m2"): 0.3, ("c1", "m2"): 0.7, ("m3", "m4"): 0.2, ("c2", "m4"): 0.8}
            edges = [
                {"u": u, "v": v, "gamma": weighted[(u, v)]} if (u, v) in weighted else (u, v)
                for u, v in edges
            ]
        return DirectedPhyNetwork(
            edges=edges, nodes=[(name, {"label": name}) for name in ("A", "B", "C", "D")]
        )

    def test_network_really_has_two_blobs(self):
        from phylozoo.core.network.dnetwork.derivations import _hybrid_blob_groups

        assert len(_hybrid_blob_groups(self._two_blob_network())) == 2

    @pytest.mark.parametrize("gamma", [False, True])
    def test_matches_definition_across_blobs(self, gamma):
        network = self._two_blob_network(gamma=gamma)
        got = displayed_splits(network)
        want = self._by_definition(network)
        assert set(got.splits) == set(want)
        for split, weight in want.items():
            assert got.get_weight(split) == pytest.approx(weight)

    def test_single_blob_matches_definition(self):
        network = DirectedPhyNetwork(
            edges=[
                ("rho", "u"),
                ("rho", "w"),
                ("u", "h1"),
                ("w", "h1"),
                ("u", "h2"),
                ("w", "h2"),
                ("u", "A"),
                ("w", "B"),
                ("h1", "C"),
                ("h2", "D"),
            ],
            nodes=[(name, {"label": name}) for name in ("A", "B", "C", "D")],
        )
        got = displayed_splits(network)
        want = self._by_definition(network)
        assert set(got.splits) == set(want)
        for split, weight in want.items():
            assert got.get_weight(split) == pytest.approx(weight)

    def test_tree_gives_its_splits_with_unit_weight(self):
        network = DirectedPhyNetwork(
            edges=[("r", "x"), ("r", "C"), ("x", "A"), ("x", "B")],
            nodes=[(name, {"label": name}) for name in ("A", "B", "C")],
        )
        got = displayed_splits(network)
        assert set(got.splits) == set(induced_splits(network).splits)
        assert all(got.get_weight(split) == pytest.approx(1.0) for split in got.splits)


class TestSwitchingMachinery:
    """The pieces shared by `distances` and `displayed_splits`."""

    @staticmethod
    def _network(gamma=False):
        return TestDistancesBlobDecomposition._two_blob_network(gamma=gamma)

    def test_hybrid_parent_edges_lists_every_parent_edge(self):
        from phylozoo.core.network.dnetwork.derivations import _hybrid_parent_edges

        network = self._network()
        parents = _hybrid_parent_edges(network)
        assert set(parents) == set(network.hybrid_nodes)
        for hybrid, edges in parents.items():
            assert len(edges) == 2
            assert all(v == hybrid for _u, v, _key in edges)

    def test_plan_groups_reference_and_masses(self):
        from phylozoo.core.network.dnetwork.derivations import _BlobSwitchings

        plan = _BlobSwitchings(self._network())
        assert len(plan.groups) == 2
        for hybrid, edges in plan.hybrid_parent_edges.items():
            assert plan.reference[hybrid] == edges[0]
        # two parent edges without gammas: 1/2 each, so every group has mass 1
        assert plan.masses == pytest.approx([1.0, 1.0])
        assert plan.total_mass == pytest.approx(1.0)
        assert plan.other_mass(0) == pytest.approx(plan.masses[1])

    def test_local_switchings_use_gamma_weights(self):
        from phylozoo.core.network.dnetwork.derivations import _BlobSwitchings

        plan = _BlobSwitchings(self._network(gamma=True))
        for index, group in enumerate(plan.groups):
            switchings = list(plan.local_switchings(index))
            assert len(switchings) == 2  # one hybrid with two parent edges
            assert sorted(weight for _choices, weight in switchings) == pytest.approx([0.25, 0.75])
            for choices, _weight in switchings:
                # only this group's hybrid deviates from the reference
                assert all(
                    choices[h] == plan.reference[h] for h in plan.reference if h not in group
                )
        assert plan.total_mass == pytest.approx(1.0)

    def test_removed_edges_are_the_unkept_parent_edges(self):
        from phylozoo.core.network.dnetwork.derivations import _BlobSwitchings

        plan = _BlobSwitchings(self._network())
        removed = plan.removed_edges(plan.reference)
        assert len(removed) == len(plan.hybrid_parent_edges)  # one dropped edge per hybrid
        assert not any(edge in removed for edge in plan.reference.values())

    def test_switching_adjacency_is_a_symmetric_tree_without_removed_edges(self):
        from phylozoo.core.network.dnetwork.derivations import _BlobSwitchings, _switching_adjacency

        network = self._network()
        plan = _BlobSwitchings(network)
        removed = plan.removed_edges(plan.reference)
        adjacency = _switching_adjacency(network, removed)
        assert set(adjacency) == set(network.nodes)
        for node, neighbours in adjacency.items():
            for neighbour, key in neighbours:
                assert (node, key) in adjacency[neighbour]
                assert (node, neighbour, key) not in removed and (
                    neighbour,
                    node,
                    key,
                ) not in removed
        # a switching is a tree: E = V - 1
        assert sum(len(neighbours) for neighbours in adjacency.values()) // 2 == len(adjacency) - 1

    def test_branch_lengths_both_orientations_and_default(self):
        from phylozoo.core.network.dnetwork.derivations import _branch_lengths

        network = DirectedPhyNetwork(
            edges=[{"u": "r", "v": "x", "branch_length": 2.5}, ("r", "C"), ("x", "A"), ("x", "B")],
            nodes=[(name, {"label": name}) for name in ("A", "B", "C")],
        )
        lengths = _branch_lengths(network)
        assert lengths[("r", "x", 0)] == 2.5 and lengths[("x", "r", 0)] == 2.5
        assert lengths[("x", "A", 0)] == 1.0 and lengths[("A", "x", 0)] == 1.0
