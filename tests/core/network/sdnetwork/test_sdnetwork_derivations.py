"""
Tests for sdnetwork derivations functions.
"""

import math

import pytest

from phylozoo.core.network.sdnetwork import MixedPhyNetwork, SemiDirectedPhyNetwork
from phylozoo.core.network.sdnetwork.derivations import k_taxon_subnetworks, subnetwork, tree_of_blobs
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
