"""
Tests for sdnetwork derivations functions.
"""

import math

import pytest

from phylozoo.core.network.sdnetwork import MixedPhyNetwork, SemiDirectedPhyNetwork
from phylozoo.core.network.sdnetwork.derivations import _switchings, k_taxon_subnetworks, subnetwork, tree_of_blobs
from phylozoo.core.network.sdnetwork.features import blobs
from phylozoo.core.network.sdnetwork.io import sdnetwork_from_mmgraph
from phylozoo.core.network.sdnetwork.transformations import suppress_2_blobs
from phylozoo.utils.validation import no_validation


class TestTreeOfBlobs:
    """Test tree_of_blobs function for MixedPhyNetwork/SemiDirectedPhyNetwork."""

    def test_basic_functionality(self) -> None:
        """Test that tree_of_blobs runs without errors and returns correct type."""
        from tests.fixtures.sd_networks import LEVEL_1_SDNETWORK_SINGLE_HYBRID
        result = tree_of_blobs(LEVEL_1_SDNETWORK_SINGLE_HYBRID)
        assert isinstance(result, type(LEVEL_1_SDNETWORK_SINGLE_HYBRID))
        assert result is not LEVEL_1_SDNETWORK_SINGLE_HYBRID  # Should return a new network

    def test_empty_network(self) -> None:
        """Empty network should return empty network."""
        mnet = MixedPhyNetwork()
        result = tree_of_blobs(mnet)
        assert isinstance(result, MixedPhyNetwork)
        assert result.number_of_nodes() == 0

    def test_single_node_network(self) -> None:
        """Single node network should remain unchanged."""
        mnet = MixedPhyNetwork(nodes=[(1, {'label': 'A'})])
        result = tree_of_blobs(mnet)
        assert isinstance(result, MixedPhyNetwork)
        assert result.number_of_nodes() == 1
        assert result.get_label(1) == 'A'

    def test_returns_copy(self) -> None:
        """Function should return a new network, not modify the original."""
        from tests.fixtures.sd_networks import LEVEL_1_SDNETWORK_SINGLE_HYBRID
        original_nodes = LEVEL_1_SDNETWORK_SINGLE_HYBRID.number_of_nodes()

        result = tree_of_blobs(LEVEL_1_SDNETWORK_SINGLE_HYBRID)

        # Original should be unchanged
        assert LEVEL_1_SDNETWORK_SINGLE_HYBRID.number_of_nodes() == original_nodes
        # Result should be a different object
        assert result is not LEVEL_1_SDNETWORK_SINGLE_HYBRID

    def test_no_2_blobs_preserves_internal_blob_count(self) -> None:
        """When there are no 2-blobs, suppress_2_blobs doesn't change internal blob count."""
        from tests.fixtures.sd_networks import LEVEL_1_SDNETWORK_SINGLE_HYBRID

        # Count internal blobs before suppress_2_blobs
        initial_blobs = blobs(LEVEL_1_SDNETWORK_SINGLE_HYBRID, trivial=False, leaves=False)
        initial_count = len(initial_blobs)

        # Apply suppress_2_blobs
        suppressed_network = suppress_2_blobs(LEVEL_1_SDNETWORK_SINGLE_HYBRID)

        # Count internal blobs after suppress_2_blobs
        after_suppress_blobs = blobs(suppressed_network, trivial=False, leaves=False)
        after_count = len(after_suppress_blobs)

        # Number of internal blobs should remain the same when there are no 2-blobs to suppress
        assert initial_count == after_count

        # But tree_of_blobs should still collapse all internal blobs
        result = tree_of_blobs(LEVEL_1_SDNETWORK_SINGLE_HYBRID)
        final_blobs = blobs(result, trivial=False, leaves=False)
        # All internal blobs should be collapsed to single vertices
        assert len(final_blobs) == 0


class TestSubnetwork:
    """Test subnetwork function for SemiDirectedPhyNetwork."""

    def test_basic_tree_subnetwork(self) -> None:
        """Test subnetwork extraction from a simple tree."""
        # Use existing fixture
        from tests.fixtures.sd_networks import LEVEL_1_SDNETWORK_SINGLE_HYBRID
        net = LEVEL_1_SDNETWORK_SINGLE_HYBRID
        # Get first two taxa
        taxa_list = sorted(net.taxa)
        if len(taxa_list) >= 2:
            requested_taxa = taxa_list[:2]
            subnet = subnetwork(net, requested_taxa)
            assert isinstance(subnet, SemiDirectedPhyNetwork)
            # Should contain the requested taxa
            for taxon in requested_taxa:
                assert taxon in subnet.taxa

    def test_empty_taxa_list(self) -> None:
        """Empty taxa list should return empty network."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 3), (5, 4), (5, 6), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (6, {'label': 'D'})]
        )
        subnet = subnetwork(net, [])
        assert isinstance(subnet, SemiDirectedPhyNetwork)
        assert subnet.number_of_nodes() == 0
        assert len(subnet.taxa) == 0

    def test_single_taxon(self) -> None:
        """Single taxon should return a network with just that leaf."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 3), (5, 4), (5, 6), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (6, {'label': 'D'})]
        )
        subnet = subnetwork(net, ['A'])
        assert isinstance(subnet, SemiDirectedPhyNetwork)
        assert sorted(subnet.taxa) == ['A']
        assert 1 in subnet._graph.nodes()

    def test_network_with_hybrid(self) -> None:
        """Test subnetwork extraction from a network with hybrid node."""
        # Create a valid network with hybrid: 5 and 6 are roots (connected via undirected edge),
        # 4 is hybrid, 8 connects to leaves
        # Need to connect 5 and 6 to form single source component, and give them degree >= 3
        net = SemiDirectedPhyNetwork(
            directed_edges=[(5, 4), (6, 4)],
            undirected_edges=[(5, 6), (4, 8), (8, 1), (8, 2), (5, 3), (6, 7)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'}), (7, {'label': 'D'})]
        )
        subnet = subnetwork(net, ['A', 'B'])
        assert isinstance(subnet, SemiDirectedPhyNetwork)
        assert sorted(subnet.taxa) == ['A', 'B']
        # Should contain the requested leaves
        assert 1 in subnet._graph.nodes()
        assert 2 in subnet._graph.nodes()
        # Network should be valid
        subnet.validate()

    def test_taxon_not_found(self) -> None:
        """Should raise ValueError if taxon is not found."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 3), (5, 4), (5, 6), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (6, {'label': 'D'})]
        )
        with pytest.raises(ValueError, match="Taxon label 'E' not found"):
            subnetwork(net, ['A', 'E'])

    def test_with_suppress_2_blobs(self) -> None:
        """Test subnetwork with suppress_2_blobs parameter."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 3), (5, 4), (5, 6), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (6, {'label': 'D'})]
        )
        subnet = subnetwork(net, ['A', 'B'], suppress_2_blobs=True)
        assert isinstance(subnet, SemiDirectedPhyNetwork)
        assert sorted(subnet.taxa) == ['A', 'B']
        subnet.validate()  # Should be valid

    def test_with_identify_parallel_edges(self) -> None:
        """Test subnetwork with identify_parallel_edges parameter."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 3), (5, 4), (5, 6), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (6, {'label': 'D'})]
        )
        subnet = subnetwork(net, ['A', 'B'], identify_parallel_edges=True)
        assert isinstance(subnet, SemiDirectedPhyNetwork)
        assert sorted(subnet.taxa) == ['A', 'B']
        subnet.validate()  # Should be valid

    def test_with_both_parameters(self) -> None:
        """Test subnetwork with both suppress_2_blobs and identify_parallel_edges."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 3), (5, 4), (5, 6), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (6, {'label': 'D'})]
        )
        subnet = subnetwork(
            net,
            ['A', 'B'],
            suppress_2_blobs=True,
            identify_parallel_edges=True
        )
        assert isinstance(subnet, SemiDirectedPhyNetwork)
        assert sorted(subnet.taxa) == ['A', 'B']
        subnet.validate()  # Should be valid

    def test_returns_copy(self) -> None:
        """Function should return a new network, not modify the original."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 3), (5, 4), (5, 6), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (6, {'label': 'D'})]
        )
        original_nodes = net.number_of_nodes()
        original_taxa = net.taxa.copy()

        subnet = subnetwork(net, ['A', 'B'])

        # Original should be unchanged
        assert net.number_of_nodes() == original_nodes
        assert net.taxa == original_taxa
        # Result should be a different object
        assert subnet is not net

    def test_all_taxa(self) -> None:
        """Subnetwork with all taxa should include all leaves."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 3), (5, 4), (5, 6), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (6, {'label': 'D'})]
        )
        subnet = subnetwork(net, ['A', 'B', 'C', 'D'])
        assert sorted(subnet.taxa) == ['A', 'B', 'C', 'D']
        # Should contain all nodes
        assert 1 in subnet._graph.nodes()
        assert 2 in subnet._graph.nodes()
        assert 4 in subnet._graph.nodes()
        assert 6 in subnet._graph.nodes()


class TestKTaxonSubnetworks:
    """Test k_taxon_subnetworks function for SemiDirectedPhyNetwork."""

    def test_k_equals_zero_returns_empty_network(self) -> None:
        """Test that k=0 returns a single empty network."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 3), (5, 4), (5, 6), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (6, {'label': 'D'})]
        )
        subnetworks = list(k_taxon_subnetworks(net, k=0))
        
        # Should have exactly 1 combination (empty set)
        assert len(subnetworks) == 1
        assert subnetworks[0].number_of_nodes() == 0
        assert subnetworks[0].number_of_edges() == 0

    def test_k_equals_all_taxa_returns_single_network(self) -> None:
        """Test that k equals number of taxa returns the original network."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 3), (5, 4), (5, 6), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (6, {'label': 'D'})]
        )
        num_taxa = len(net.taxa)
        
        subnetworks = list(k_taxon_subnetworks(net, k=num_taxa))
        
        # Should have exactly 1 combination (all taxa)
        assert len(subnetworks) == 1
        subnet = subnetworks[0]
        
        # Should have the same taxa
        assert subnet.taxa == net.taxa

    def test_k_equals_one_returns_all_single_taxon_networks(self) -> None:
        """Test that k=1 returns all single-taxon subnetworks."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 3), (5, 4), (5, 6), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (6, {'label': 'D'})]
        )
        all_taxa = list(net.taxa)
        num_taxa = len(all_taxa)
        
        subnetworks = list(k_taxon_subnetworks(net, k=1))
        
        # Should have exactly num_taxa combinations
        assert len(subnetworks) == num_taxa
        
        # Each subnetwork should have exactly 1 leaf
        for subnet in subnetworks:
            assert len(subnet.taxa) == 1
        
        # All taxa should be represented
        subnet_taxa = {list(subnet.taxa)[0] for subnet in subnetworks}
        assert subnet_taxa == set(all_taxa)

    def test_k_equals_two_returns_all_pairs(self) -> None:
        """Test that k=2 returns all pairs of taxa."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 3), (5, 4), (5, 6), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (6, {'label': 'D'})]
        )
        num_taxa = len(net.taxa)
        expected_count = math.comb(num_taxa, 2)
        
        subnetworks = list(k_taxon_subnetworks(net, k=2))
        
        # Should have C(num_taxa, 2) combinations
        assert len(subnetworks) == expected_count
        
        # Each subnetwork should have exactly 2 leaves
        for subnet in subnetworks:
            assert len(subnet.taxa) == 2

    def test_negative_k_raises_error(self) -> None:
        """Test that negative k raises ValueError."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 3), (5, 4), (5, 6), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (6, {'label': 'D'})]
        )
        with pytest.raises(ValueError, match="must be non-negative"):
            list(k_taxon_subnetworks(net, k=-1))

    def test_k_greater_than_taxa_raises_error(self) -> None:
        """Test that k > number of taxa raises ValueError."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 3), (5, 4), (5, 6), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (6, {'label': 'D'})]
        )
        num_taxa = len(net.taxa)
        
        with pytest.raises(ValueError, match="cannot exceed"):
            list(k_taxon_subnetworks(net, k=num_taxa + 1))

    def test_all_subnetworks_are_valid(self) -> None:
        """Test that all generated subnetworks are valid."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 3), (5, 4), (5, 6), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (6, {'label': 'D'})]
        )
        num_taxa = len(net.taxa)
        
        # Test for various k values
        for k in [1, 2, min(3, num_taxa // 2)]:
            subnetworks = list(k_taxon_subnetworks(net, k=k))
            
            for subnet in subnetworks:
                # Each subnetwork should be valid
                assert subnet.number_of_nodes() >= 0
                assert subnet.number_of_edges() >= 0
                assert len(subnet.taxa) == k
                
                # All leaves should have labels
                for taxon in subnet.taxa:
                    assert taxon is not None
                
                # Network should be valid
                subnet.validate()

    def test_with_suppress_2_blobs(self) -> None:
        """Test k_taxon_subnetworks with suppress_2_blobs parameter."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 3), (5, 4), (5, 6), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (6, {'label': 'D'})]
        )
        subnetworks = list(k_taxon_subnetworks(net, k=2, suppress_2_blobs=True))
        
        assert len(subnetworks) == math.comb(len(net.taxa), 2)
        for subnet in subnetworks:
            assert len(subnet.taxa) == 2
            subnet.validate()

    def test_with_identify_parallel_edges(self) -> None:
        """Test k_taxon_subnetworks with identify_parallel_edges parameter."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 3), (5, 4), (5, 6), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (6, {'label': 'D'})]
        )
        subnetworks = list(k_taxon_subnetworks(net, k=2, identify_parallel_edges=True))
        
        assert len(subnetworks) == math.comb(len(net.taxa), 2)
        for subnet in subnetworks:
            assert len(subnet.taxa) == 2
            subnet.validate()

    def test_with_both_parameters(self) -> None:
        """Test k_taxon_subnetworks with both parameters."""
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(5, 3), (5, 4), (5, 6), (3, 1), (3, 2)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (6, {'label': 'D'})]
        )
        subnetworks = list(k_taxon_subnetworks(
            net,
            k=2,
            suppress_2_blobs=True,
            identify_parallel_edges=True
        ))
        
        assert len(subnetworks) == math.comb(len(net.taxa), 2)
        for subnet in subnetworks:
            assert len(subnet.taxa) == 2
            subnet.validate()

    def test_with_hybrid_network(self) -> None:
        """Test with a network containing hybrid nodes."""
        from tests.fixtures.sd_networks import LEVEL_1_SDNETWORK_SINGLE_HYBRID
        net = LEVEL_1_SDNETWORK_SINGLE_HYBRID
        num_taxa = len(net.taxa)
        k = min(2, num_taxa)
        
        subnetworks = list(k_taxon_subnetworks(net, k=k))
        
        expected_count = math.comb(num_taxa, k)
        assert len(subnetworks) == expected_count
        
        # All subnetworks should be valid
        for subnet in subnetworks:
            assert len(subnet.taxa) == k
            assert subnet.number_of_nodes() > 0
            subnet.validate()


class TestSwitchings:
    """Test _switchings function for SemiDirectedPhyNetwork."""

    def test_tree_network_single_switching(self) -> None:
        """Tree network (no hybrid nodes) should yield exactly one switching."""
        # Create a valid tree network (internal node needs degree >= 3)
        net = SemiDirectedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
        )
        switchings = list(_switchings(net))
        assert len(switchings) == 1
        assert switchings[0].number_of_edges() == 3
        
        # With probability=True, should have probability 1.0
        switchings_with_prob = list(_switchings(net, probability=True))
        assert len(switchings_with_prob) == 1
        assert switchings_with_prob[0]._directed.graph.get('probability') == 1.0

    def test_single_hybrid_two_parents(self) -> None:
        """Network with one hybrid node and two parent edges should yield two switchings."""
        from tests.fixtures.sd_networks import LEVEL_1_SDNETWORK_SINGLE_HYBRID
        net = LEVEL_1_SDNETWORK_SINGLE_HYBRID
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
        # Use a simpler approach: connect all parents to ensure single source component
        net = SemiDirectedPhyNetwork(
            directed_edges=[(5, 4), (6, 4), (7, 4)],
            undirected_edges=[
                (4, 1),  # Hybrid outgoing edge
                (5, 6), (6, 7),  # Connect parents to ensure single source component
                (5, 8), (5, 9),  # Tree node 5
                (6, 10), (6, 11),  # Tree node 6
                (7, 12), (7, 13),  # Tree node 7
            ],
            nodes=[
                (1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'}),
                (10, {'label': 'D'}), (11, {'label': 'E'}), (12, {'label': 'F'}), (13, {'label': 'G'})
            ]
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
        # Use existing fixture with two hybrids
        from tests.fixtures.sd_networks import LEVEL_1_SDNETWORK_TWO_HYBRIDS_SEPARATE
        net = LEVEL_1_SDNETWORK_TWO_HYBRIDS_SEPARATE
        switchings = list(_switchings(net))
        assert len(switchings) == 4  # 2 * 2 = 4
        
        # Find the hybrid nodes
        hybrid_nodes = net.hybrid_nodes
        assert len(hybrid_nodes) == 2
        hybrids = sorted(hybrid_nodes)
        
        # Each switching should have exactly one parent edge per hybrid
        for sw in switchings:
            for hybrid in hybrids:
                parent_edges = list(sw.incident_parent_edges(hybrid, keys=True))
                assert len(parent_edges) == 1

    def test_switchings_preserve_undirected_edges(self) -> None:
        """Switchings should preserve all undirected edges."""
        # Use existing fixture
        from tests.fixtures.sd_networks import LEVEL_1_SDNETWORK_SINGLE_HYBRID
        net = LEVEL_1_SDNETWORK_SINGLE_HYBRID
        switchings = list(_switchings(net))
        assert len(switchings) == 2
        
        # Get all undirected edges from original network
        original_undirected = set(net._graph._undirected.edges(keys=True))
        
        # All switchings should preserve all undirected edges
        for sw in switchings:
            switching_undirected = set(sw._undirected.edges(keys=True))
            assert original_undirected == switching_undirected

    def test_switchings_are_copies(self) -> None:
        """Switchings should be independent copies, not references."""
        from tests.fixtures.sd_networks import LEVEL_1_SDNETWORK_SINGLE_HYBRID
        net = LEVEL_1_SDNETWORK_SINGLE_HYBRID
        switchings = list(_switchings(net))
        
        # Modify one switching and verify others are unchanged
        if len(switchings) > 0:
            original_edge_count = switchings[0].number_of_edges()
            # Remove an undirected edge (not a hybrid edge)
            for u, v in switchings[0].undirected_edges_iter():
                switchings[0].remove_edge(u, v)
                break
            # Other switchings should be unchanged
            for sw in switchings[1:]:
                assert sw.number_of_edges() == original_edge_count

    def test_probability_with_gamma_values(self) -> None:
        """Test probability calculation when gamma values are present."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {'u': 5, 'v': 4, 'gamma': 0.6},
                {'u': 6, 'v': 4, 'gamma': 0.4}
            ],
            undirected_edges=[
                (5, 3), (5, 6), (6, 7), (4, 8), (8, 1), (8, 2)
            ],
            nodes=[
                (3, {'label': 'C'}), (7, {'label': 'D'}),
                (1, {'label': 'A'}), (2, {'label': 'B'})
            ]
        )
        switchings = list(_switchings(net, probability=True))
        assert len(switchings) == 2
        
        # Check probabilities
        probs = [sw._directed.graph.get('probability') for sw in switchings]
        assert 0.6 in probs
        assert 0.4 in probs
        # Probabilities should sum to 1.0
        assert abs(sum(probs) - 1.0) < 1e-10

    def test_probability_without_gamma_values(self) -> None:
        """Test probability calculation when no gamma values are present."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                (5, 4), (6, 4)  # No gamma values
            ],
            undirected_edges=[
                (5, 3), (5, 6), (6, 7), (4, 8), (8, 1), (8, 2)
            ],
            nodes=[
                (3, {'label': 'C'}), (7, {'label': 'D'}),
                (1, {'label': 'A'}), (2, {'label': 'B'})
            ]
        )
        switchings = list(_switchings(net, probability=True))
        assert len(switchings) == 2
        
        # Without gamma, each edge should have probability 1/2 (indegree is 2)
        for sw in switchings:
            prob = sw._directed.graph.get('probability')
            assert prob is not None
            assert abs(prob - 0.5) < 1e-10
        # Probabilities should sum to 1.0
        total_prob = sum(sw._directed.graph.get('probability') for sw in switchings)
        assert abs(total_prob - 1.0) < 1e-10

    def test_probability_with_multiple_hybrids(self) -> None:
        """Test probability calculation with multiple hybrid nodes."""
        from tests.fixtures.sd_networks import LEVEL_1_SDNETWORK_TWO_HYBRIDS_SEPARATE
        net = LEVEL_1_SDNETWORK_TWO_HYBRIDS_SEPARATE
        switchings = list(_switchings(net, probability=True))
        assert len(switchings) == 4  # 2 * 2 = 4
        
        # Each switching should have a probability
        for sw in switchings:
            prob = sw._directed.graph.get('probability')
            assert prob is not None
            assert 0.0 < prob <= 1.0
        
        # Probabilities should sum to 1.0
        total_prob = sum(sw._directed.graph.get('probability') for sw in switchings)
        assert abs(total_prob - 1.0) < 1e-10

    def test_probability_false_no_attribute(self) -> None:
        """Test that probability=False does not add probability attribute."""
        net = SemiDirectedPhyNetwork(
            directed_edges=[
                {'u': 5, 'v': 4, 'gamma': 0.6},
                {'u': 6, 'v': 4, 'gamma': 0.4}
            ],
            undirected_edges=[
                (5, 3), (5, 6), (6, 7), (4, 8), (8, 1), (8, 2)
            ],
            nodes=[
                (3, {'label': 'C'}), (7, {'label': 'D'}),
                (1, {'label': 'A'}), (2, {'label': 'B'})
            ]
        )
        switchings = list(_switchings(net, probability=False))
        
        # No probability attribute should be present
        for sw in switchings:
            assert 'probability' not in sw._directed.graph
