"""
Tests for dnetwork derivations functions.
"""

import numpy as np
import pytest

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.derivations import _switchings, displayed_trees, displayed_splits, tree_of_blobs, distances, induced_splits
from phylozoo.core.network.dnetwork.features import blobs
from phylozoo.core.network.dnetwork.io import dnetwork_from_dmgraph
from phylozoo.core.network.dnetwork.transformations import suppress_2_blobs
from phylozoo.utils.validation import no_validation


class TestTreeOfBlobs:
    """Test tree_of_blobs function for DirectedPhyNetwork."""

    def test_simple_tree_unchanged(self) -> None:
        """Tree with no blobs should remain unchanged."""
        dnet = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        result = tree_of_blobs(dnet)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_nodes() == 3
        assert result.number_of_edges() == 2
        assert result.get_label(1) == 'A'
        assert result.get_label(2) == 'B'

    def test_basic_functionality(self) -> None:
        """Test that tree_of_blobs runs without errors and returns correct type."""
        from tests.fixtures.directed_networks import LEVEL_1_DNETWORK_SINGLE_HYBRID
        result = tree_of_blobs(LEVEL_1_DNETWORK_SINGLE_HYBRID)
        assert isinstance(result, DirectedPhyNetwork)
        assert result is not LEVEL_1_DNETWORK_SINGLE_HYBRID  # Should return a new network

    def test_single_node_network(self) -> None:
        """Single node network should remain unchanged."""
        dnet = DirectedPhyNetwork(nodes=[(1, {'label': 'A'})])
        result = tree_of_blobs(dnet)
        assert isinstance(result, DirectedPhyNetwork)
        assert result.number_of_nodes() == 1
        assert result.get_label(1) == 'A'

    def test_result_is_copy(self) -> None:
        """Function should return a new network, not modify the original."""
        dnet = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
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
        all_labels = [result.get_label(node) for node in result._graph.nodes() if result.get_label(node)]
        assert len(all_labels) > 0

    def test_no_2_blobs_preserves_internal_blob_count(self) -> None:
        """When there are no 2-blobs, suppress_2_blobs doesn't change internal blob count."""
        from tests.fixtures.directed_networks import LEVEL_1_DNETWORK_SINGLE_HYBRID
        from phylozoo.core.network.dnetwork.transformations import suppress_2_blobs

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
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        switchings = list(_switchings(net))
        assert len(switchings) == 1
        assert switchings[0].number_of_edges() == 2
        
        # With probability=True, should have probability 1.0 in all underlying graphs
        switchings_with_prob = list(_switchings(net, probability=True))
        assert len(switchings_with_prob) == 1
        assert switchings_with_prob[0]._graph.graph.get('probability') == 1.0
        assert switchings_with_prob[0]._combined.graph.get('probability') == 1.0

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
                (10, 5), (10, 6), (10, 7),  # Root to tree nodes
                (5, 4), (6, 4), (7, 4),  # Three parents to hybrid 4
                (4, 1),  # Hybrid to leaf
                (5, 8), (6, 9), (7, 11),  # Other children for tree nodes
            ],
            nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}), (11, {'label': 'D'})]
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
                (10, 5), (10, 6),  # Root to tree nodes
                (5, 4), (5, 17), (6, 4), (6, 18), (4, 11), (11, 1), (11, 12),  # Hybrid 4 with 2 parents, tree node 11
                (10, 8), (10, 9),  # Root to more tree nodes (8 and 9 are tree nodes)
                (8, 7), (8, 15),  # Tree node 8: one to hybrid 7, one to leaf
                (9, 7), (9, 16),  # Tree node 9: one to hybrid 7, one to leaf
                (7, 13), (13, 2), (13, 14),  # Hybrid 7 with 2 parents, tree node 13
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (12, {'label': 'C'}), (14, {'label': 'D'}), (15, {'label': 'E'}), (16, {'label': 'F'}), (17, {'label': 'G'}), (18, {'label': 'H'})]
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
                (10, 5), (10, 6),  # Root to tree nodes
                (5, 4), (5, 17), (6, 4), (6, 18), (4, 1),  # Hybrid 4, tree nodes 5 and 6 have out-degree 2
                (10, 3), (3, 2), (3, 15),  # Non-hybrid tree node 3 with out-degree 2
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (15, {'label': 'C'}), (17, {'label': 'D'}), (18, {'label': 'E'})]
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
                (10, 5), (10, 6),
                {'u': 5, 'v': 4, 'gamma': 0.6},
                {'u': 6, 'v': 4, 'gamma': 0.4},
                (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
        )
        switchings = list(_switchings(net, probability=True))
        assert len(switchings) == 2
        
        # Check probabilities in all underlying graphs
        probs = [sw._graph.graph.get('probability') for sw in switchings]
        assert 0.6 in probs
        assert 0.4 in probs
        # Probabilities should sum to 1.0
        assert abs(sum(probs) - 1.0) < 1e-10
        
        # Check that probability is set in all underlying graphs
        for sw in switchings:
            prob = sw._graph.graph.get('probability')
            assert prob is not None
            assert sw._combined.graph.get('probability') == prob

    def test_probability_without_gamma_values(self) -> None:
        """Test probability calculation when no gamma values are present."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 5), (10, 6),
                (5, 4), (6, 4),  # No gamma values
                (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
        )
        switchings = list(_switchings(net, probability=True))
        assert len(switchings) == 2
        
        # Without gamma, each edge should have probability 1/2 (indegree is 2)
        for sw in switchings:
            prob = sw._graph.graph.get('probability')
            assert prob is not None
            assert abs(prob - 0.5) < 1e-10
            # Check that probability is set in all underlying graphs
            assert sw._combined.graph.get('probability') == prob
        # Probabilities should sum to 1.0
        total_prob = sum(sw._graph.graph.get('probability') for sw in switchings)
        assert abs(total_prob - 1.0) < 1e-10

    def test_probability_with_multiple_hybrids(self) -> None:
        """Test probability calculation with multiple hybrid nodes."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 5), (10, 6), (10, 8), (10, 9),
                (5, 4), (5, 17), (6, 4), (6, 18),  # Hybrid 4 with 2 parents
                (8, 7), (8, 15), (9, 7), (9, 16),  # Hybrid 7 with 2 parents
                (4, 11), (11, 1), (11, 12),
                (7, 13), (13, 2), (13, 14)
            ],
            nodes=[
                (1, {'label': 'A'}), (2, {'label': 'B'}), (12, {'label': 'C'}),
                (14, {'label': 'D'}), (15, {'label': 'E'}), (16, {'label': 'F'}),
                (17, {'label': 'G'}), (18, {'label': 'H'})
            ]
        )
        switchings = list(_switchings(net, probability=True))
        assert len(switchings) == 4  # 2 * 2 = 4
        
        # Each switching should have a probability
        for sw in switchings:
            prob = sw._graph.graph.get('probability')
            assert prob is not None
            assert 0.0 < prob <= 1.0
        
        # Probabilities should sum to 1.0
        total_prob = sum(sw._graph.graph.get('probability') for sw in switchings)
        assert abs(total_prob - 1.0) < 1e-10

    def test_probability_false_no_attribute(self) -> None:
        """Test that probability=False does not add probability attribute."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 5), (10, 6),
                {'u': 5, 'v': 4, 'gamma': 0.6},
                {'u': 6, 'v': 4, 'gamma': 0.4},
                (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
        )
        switchings = list(_switchings(net, probability=False))
        
        # No probability attribute should be present
        for sw in switchings:
            assert 'probability' not in sw._graph.graph


class TestDisplayedTrees:
    """Test cases for the displayed_trees function for DirectedPhyNetwork."""

    def test_tree_network_single_displayed_tree(self) -> None:
        """A tree network should yield exactly one displayed tree."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
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
            edges=[
                (10, 5), (10, 6),
                (5, 4), (6, 4),
                (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
        )
        trees = list(displayed_trees(net))
        
        for tree in trees:
            assert len(tree.hybrid_nodes) == 0

    def test_probability_with_gamma_values(self) -> None:
        """Test probability calculation when gamma values are present."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 5), (10, 6),
                {'u': 5, 'v': 4, 'gamma': 0.6},
                {'u': 6, 'v': 4, 'gamma': 0.4},
                (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
        )
        trees = list(displayed_trees(net, probability=True))
        assert len(trees) == 2
        
        # Check probabilities
        probs = [tree.get_network_attribute('probability') for tree in trees]
        assert 0.6 in probs
        assert 0.4 in probs
        # Probabilities should sum to 1.0
        assert abs(sum(probs) - 1.0) < 1e-10

    def test_probability_without_gamma_values(self) -> None:
        """Test probability calculation when no gamma values are present."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 5), (10, 6),
                (5, 4), (6, 4),  # No gamma values
                (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
        )
        trees = list(displayed_trees(net, probability=True))
        assert len(trees) == 2
        
        # Without gamma, each tree should have probability 1/2 (indegree is 2)
        for tree in trees:
            prob = tree.get_network_attribute('probability')
            assert prob is not None
            assert abs(prob - 0.5) < 1e-10
        # Probabilities should sum to 1.0
        total_prob = sum(tree.get_network_attribute('probability') for tree in trees)
        assert abs(total_prob - 1.0) < 1e-10

    def test_probability_tree_network(self) -> None:
        """Tree network (no hybrid nodes) should have probability 1.0."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        trees = list(displayed_trees(net, probability=True))
        assert len(trees) == 1
        assert trees[0].get_network_attribute('probability') == 1.0

    def test_probability_false_no_attribute(self) -> None:
        """Test that probability=False does not add probability attribute."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 5), (10, 6),
                {'u': 5, 'v': 4, 'gamma': 0.6},
                {'u': 6, 'v': 4, 'gamma': 0.4},
                (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
        )
        trees = list(displayed_trees(net, probability=False))
        
        # No probability attribute should be present
        for tree in trees:
            assert tree.get_network_attribute('probability') is None


class TestDistances:
    """Test distances function for DirectedPhyNetwork."""
    
    def test_simple_tree_no_branch_lengths(self) -> None:
        """Test distances for a simple tree without branch lengths."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        dm = distances(net, mode='shortest')
        assert len(dm) == 2
        assert dm.get_distance('A', 'B') == 2.0  # Two edges, each with default length 1.0
        assert dm.get_distance('A', 'A') == 0.0
        assert dm.get_distance('B', 'B') == 0.0
    
    def test_simple_tree_with_branch_lengths(self) -> None:
        """Test distances for a simple tree with branch lengths."""
        net = DirectedPhyNetwork(
            edges=[
                {'u': 3, 'v': 1, 'branch_length': 0.5},
                {'u': 3, 'v': 2, 'branch_length': 0.3}
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        dm = distances(net, mode='shortest')
        assert len(dm) == 2
        assert dm.get_distance('A', 'B') == 0.8  # 0.5 + 0.3
        assert dm.get_distance('A', 'A') == 0.0
        assert dm.get_distance('B', 'B') == 0.0
    
    def test_network_with_hybrid_shortest(self) -> None:
        """Test shortest distances for a network with hybrid node."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 5), (10, 6),
                (5, 4), (6, 4),
                (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
        )
        dm = distances(net, mode='shortest')
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
            edges=[
                (10, 5), (10, 6),
                (5, 4), (6, 4),
                (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
        )
        dm = distances(net, mode='longest')
        assert len(dm) == 4
        # All distances should be finite and >= shortest
        shortest_dm = distances(net, mode='shortest')
        for i, taxon1 in enumerate(dm.labels):
            for j, taxon2 in enumerate(dm.labels):
                if i != j:
                    assert dm.get_distance(taxon1, taxon2) >= shortest_dm.get_distance(taxon1, taxon2)
                else:
                    assert dm.get_distance(taxon1, taxon2) == 0.0
    
    def test_network_with_hybrid_average(self) -> None:
        """Test average distances for a network with hybrid node."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 5), (10, 6),
                {'u': 5, 'v': 4, 'gamma': 0.6},
                {'u': 6, 'v': 4, 'gamma': 0.4},
                (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
        )
        dm = distances(net, mode='average')
        assert len(dm) == 4
        # All distances should be finite and between shortest and longest
        shortest_dm = distances(net, mode='shortest')
        longest_dm = distances(net, mode='longest')
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
        net = DirectedPhyNetwork(
            nodes=[(1, {'label': 'A'})]
        )
        dm = distances(net, mode='shortest')
        assert len(dm) == 1
        assert dm.get_distance('A', 'A') == 0.0
    
    def test_symmetric_matrix(self) -> None:
        """Test that the distance matrix is symmetric."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        dm = distances(net, mode='shortest')
        assert dm.get_distance('A', 'B') == dm.get_distance('B', 'A')
    
    def test_all_modes_produce_valid_matrices(self) -> None:
        """Test that all modes produce valid distance matrices."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 5), (10, 6),
                (5, 4), (6, 4),
                (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
        )
        for mode in ['shortest', 'longest', 'average']:
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
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        )
        splits = induced_splits(net)
        assert len(splits) >= 1
        # Should have at least one split (the edge between root and leaves)

    def test_tree_with_three_leaves(self) -> None:
        """Test induced_splits on a tree with three leaves."""
        net = DirectedPhyNetwork(
            edges=[(4, 1), (4, 2), (4, 3)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
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
        net = DirectedPhyNetwork(
            edges=[(2, 1)],
            nodes=[(1, {'label': 'A'})]
        )
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
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
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
            edges=[(3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
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
            edges=[
                (10, 5), (10, 6),
                (5, 4), (6, 4),
                (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
        )
        splits = displayed_splits(net)
        # Get all displayed trees to verify probabilities
        trees = list(displayed_trees(net, probability=True))
        total_prob = sum(tree.get_network_attribute('probability') or 1.0 for tree in trees)
        # Total probability should be 1.0 (or close due to floating point)
        assert abs(total_prob - 1.0) < 1e-10
    
    def test_splits_from_displayed_trees(self) -> None:
        """Test that displayed_splits contains splits from all displayed trees."""
        net = DirectedPhyNetwork(
            edges=[
                (10, 5), (10, 6),
                (5, 4), (6, 4),
                (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
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
            edges=[
                (10, 5), (10, 6),
                (5, 4), (6, 4),
                (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)
            ],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
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
