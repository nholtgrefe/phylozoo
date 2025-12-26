"""
Tests for the k_taxon_subnetworks function in dnetwork derivations.
"""

import math

import pytest

from phylozoo.core.network.dnetwork.derivations import k_taxon_subnetworks
from phylozoo.core.network.dnetwork.features import cut_edges
from tests.fixtures import directed_networks


class TestKTaxonSubnetworksBasic:
    """Basic tests for k_taxon_subnetworks function."""

    def test_k_equals_zero_returns_empty_network(self):
        """Test that k=0 returns a single empty network."""
        network = directed_networks.DTREE_SMALL_BINARY
        subnetworks = list(k_taxon_subnetworks(network, k=0))

        # Should have exactly 1 combination (empty set)
        assert len(subnetworks) == 1
        assert subnetworks[0].number_of_nodes() == 0
        assert subnetworks[0].number_of_edges() == 0

    def test_k_equals_all_taxa_returns_single_network(self):
        """Test that k equals number of taxa returns the original network."""
        network = directed_networks.DTREE_SMALL_BINARY
        num_taxa = len(network.taxa)

        subnetworks = list(k_taxon_subnetworks(network, k=num_taxa))

        # Should have exactly 1 combination (all taxa)
        assert len(subnetworks) == 1
        subnet = subnetworks[0]

        # Should be the same as the original network
        assert subnet.number_of_nodes() == network.number_of_nodes()
        assert subnet.number_of_edges() == network.number_of_edges()
        assert subnet.leaves == network.leaves

    def test_k_equals_one_returns_all_single_taxon_networks(self):
        """Test that k=1 returns all single-taxon subnetworks."""
        network = directed_networks.DTREE_SMALL_BINARY
        all_taxa = list(network.taxa)
        num_taxa = len(all_taxa)

        subnetworks = list(k_taxon_subnetworks(network, k=1))

        # Should have exactly num_taxa combinations
        assert len(subnetworks) == num_taxa

        # Each subnetwork should have exactly 1 leaf
        for subnet in subnetworks:
            assert len(subnet.leaves) == 1

        # All taxa should be represented
        subnet_taxa = {subnet.get_label(next(iter(subnet.leaves))) for subnet in subnetworks}
        assert subnet_taxa == set(all_taxa)

    def test_k_equals_two_returns_all_pairs(self):
        """Test that k=2 returns all pairs of taxa."""
        network = directed_networks.DTREE_SMALL_BINARY
        num_taxa = len(network.taxa)
        expected_count = math.comb(num_taxa, 2)

        subnetworks = list(k_taxon_subnetworks(network, k=2))

        # Should have C(num_taxa, 2) combinations
        assert len(subnetworks) == expected_count

        # Each subnetwork should have exactly 2 leaves
        for subnet in subnetworks:
            assert len(subnet.leaves) == 2

    def test_negative_k_raises_error(self):
        """Test that negative k raises ValueError."""
        network = directed_networks.DTREE_SMALL_BINARY
        with pytest.raises(ValueError, match="must be non-negative"):
            list(k_taxon_subnetworks(network, k=-1))

    def test_k_greater_than_taxa_raises_error(self):
        """Test that k > number of taxa raises ValueError."""
        network = directed_networks.DTREE_SMALL_BINARY
        num_taxa = len(network.taxa)

        with pytest.raises(ValueError, match="cannot exceed"):
            list(k_taxon_subnetworks(network, k=num_taxa + 1))

    def test_all_subnetworks_are_valid(self):
        """Test that all generated subnetworks are valid."""
        network = directed_networks.DTREE_MEDIUM_BINARY
        num_taxa = len(network.taxa)

        # Test for various k values
        for k in [1, 2, min(3, num_taxa // 2)]:
            subnetworks = list(k_taxon_subnetworks(network, k=k))

            for subnet in subnetworks:
                # Each subnetwork should be valid
                assert subnet.number_of_nodes() >= 0
                assert subnet.number_of_edges() >= 0
                assert len(subnet.leaves) == k

                # All leaves should have labels
                for leaf in subnet.leaves:
                    assert subnet.get_label(leaf) is not None


class TestKTaxonSubnetworksWithOptions:
    """Tests for k_taxon_subnetworks with different options."""

    def test_with_suppress_2_blobs(self):
        """Test k_taxon_subnetworks with suppress_2_blobs option."""
        network = directed_networks.LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE
        num_taxa = len(network.taxa)
        k = min(2, num_taxa)

        subnetworks = list(
            k_taxon_subnetworks(network, k=k, suppress_2_blobs=True)
        )

        # All subnetworks should be valid
        assert len(subnetworks) == math.comb(num_taxa, k)
        for subnet in subnetworks:
            assert len(subnet.leaves) == k

    def test_with_identify_parallel_edges(self):
        """Test k_taxon_subnetworks with identify_parallel_edges option."""
        network = directed_networks.LEVEL_1_DNETWORK_PARALLEL_EDGES
        num_taxa = len(network.taxa)
        k = min(2, num_taxa)

        subnetworks = list(
            k_taxon_subnetworks(network, k=k, identify_parallel_edges=True)
        )

        # All subnetworks should be valid and have no parallel edges
        assert len(subnetworks) == math.comb(num_taxa, k)
        for subnet in subnetworks:
            assert len(subnet.leaves) == k
            # Check that parallel edges are identified (if any existed)
            # This is a basic check - the identify_parallel_edges function
            # should have been applied

    def test_with_make_lsa(self):
        """Test k_taxon_subnetworks with make_lsa option."""
        network = directed_networks.LEVEL_1_DNETWORK_SINGLE_HYBRID
        num_taxa = len(network.taxa)
        k = min(2, num_taxa)

        subnetworks = list(k_taxon_subnetworks(network, k=k, make_lsa=True))

        # All subnetworks should be valid and LSA networks
        assert len(subnetworks) == math.comb(num_taxa, k)
        for subnet in subnetworks:
            assert len(subnet.leaves) == k
            # LSA networks have root == LSA_node
            assert subnet.root_node == subnet.LSA_node

    def test_with_all_options(self):
        """Test k_taxon_subnetworks with all options enabled."""
        network = directed_networks.LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE
        num_taxa = len(network.taxa)
        k = min(2, num_taxa)

        subnetworks = list(
            k_taxon_subnetworks(
                network,
                k=k,
                suppress_2_blobs=True,
                identify_parallel_edges=True,
                make_lsa=True,
            )
        )

        # All subnetworks should be valid
        assert len(subnetworks) == math.comb(num_taxa, k)
        for subnet in subnetworks:
            assert len(subnet.leaves) == k
            assert subnet.root_node == subnet.LSA_node


class TestKTaxonSubnetworksWithVariousNetworks:
    """Tests with various fixture networks."""

    def test_with_small_tree(self):
        """Test with a small tree network."""
        network = directed_networks.DTREE_SMALL_BINARY
        num_taxa = len(network.taxa)

        # Test all possible k values
        for k in range(num_taxa + 1):
            subnetworks = list(k_taxon_subnetworks(network, k=k))
            expected_count = math.comb(num_taxa, k)
            assert len(subnetworks) == expected_count

            for subnet in subnetworks:
                assert len(subnet.leaves) == k

    def test_with_medium_tree(self):
        """Test with a medium tree network."""
        network = directed_networks.DTREE_MEDIUM_BINARY
        num_taxa = len(network.taxa)

        # Test a few k values (to avoid too many combinations)
        for k in [1, 2, 3, min(5, num_taxa)]:
            subnetworks = list(k_taxon_subnetworks(network, k=k))
            expected_count = math.comb(num_taxa, k)
            assert len(subnetworks) == expected_count

            for subnet in subnetworks:
                assert len(subnet.leaves) == k

    def test_with_hybrid_network(self):
        """Test with a network containing hybrid nodes."""
        network = directed_networks.LEVEL_1_DNETWORK_SINGLE_HYBRID
        num_taxa = len(network.taxa)
        k = min(2, num_taxa)

        subnetworks = list(k_taxon_subnetworks(network, k=k))

        expected_count = math.comb(num_taxa, k)
        assert len(subnetworks) == expected_count

        # All subnetworks should be valid
        for subnet in subnetworks:
            assert len(subnet.leaves) == k
            assert subnet.number_of_nodes() > 0

    def test_with_multiple_blobs_network(self):
        """Test with a network containing multiple blobs."""
        network = directed_networks.LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE
        num_taxa = len(network.taxa)
        k = min(2, num_taxa)

        subnetworks = list(k_taxon_subnetworks(network, k=k))

        expected_count = math.comb(num_taxa, k)
        assert len(subnetworks) == expected_count

        # All subnetworks should be valid
        for subnet in subnetworks:
            assert len(subnet.leaves) == k


class TestKTaxonSubnetworksGenerator:
    """Tests for generator behavior."""

    def test_is_generator(self):
        """Test that k_taxon_subnetworks is a generator."""
        network = directed_networks.DTREE_SMALL_BINARY
        gen = k_taxon_subnetworks(network, k=1)

        # Should be a generator (has __iter__ and __next__)
        assert hasattr(gen, "__iter__")
        assert hasattr(gen, "__next__")

        # Should yield values lazily
        first = next(gen)
        assert first.number_of_nodes() > 0

    def test_generator_exhaustion(self):
        """Test that generator can be exhausted."""
        network = directed_networks.DTREE_SMALL_BINARY
        num_taxa = len(network.taxa)
        gen = k_taxon_subnetworks(network, k=1)

        # Consume all values
        subnetworks = list(gen)
        assert len(subnetworks) == num_taxa

        # Generator should be exhausted
        with pytest.raises(StopIteration):
            next(gen)

    def test_generator_independence(self):
        """Test that multiple generators are independent."""
        network = directed_networks.DTREE_SMALL_BINARY
        gen1 = k_taxon_subnetworks(network, k=1)
        gen2 = k_taxon_subnetworks(network, k=1)

        # Consuming from one shouldn't affect the other
        first_from_gen1 = next(gen1)
        first_from_gen2 = next(gen2)

        # Both should yield valid networks
        assert first_from_gen1.number_of_nodes() > 0
        assert first_from_gen2.number_of_nodes() > 0

        # They might be the same network (same taxa), but generators are independent
        assert len(first_from_gen1.leaves) == 1
        assert len(first_from_gen2.leaves) == 1


class TestKTaxonSubnetworksEdgeCases:
    """Edge case tests."""

    def test_single_taxon_network(self):
        """Test with a network that has only one taxon."""
        network = directed_networks.DTREE_SINGLE_NODE
        num_taxa = len(network.taxa)

        # k=0 should return empty network
        subnetworks = list(k_taxon_subnetworks(network, k=0))
        assert len(subnetworks) == 1
        assert subnetworks[0].number_of_nodes() == 0

        # k=1 should return the original network
        subnetworks = list(k_taxon_subnetworks(network, k=1))
        assert len(subnetworks) == 1
        assert subnetworks[0].number_of_nodes() == network.number_of_nodes()

    def test_empty_network(self):
        """Test with an empty network."""
        network = directed_networks.DTREE_EMPTY

        # k=0 should return empty network
        subnetworks = list(k_taxon_subnetworks(network, k=0))
        assert len(subnetworks) == 1
        assert subnetworks[0].number_of_nodes() == 0

        # k>0 should raise error
        with pytest.raises(ValueError, match="cannot exceed"):
            list(k_taxon_subnetworks(network, k=1))

    def test_all_combinations_unique(self):
        """Test that all combinations are unique."""
        network = directed_networks.DTREE_MEDIUM_BINARY
        num_taxa = len(network.taxa)
        k = min(2, num_taxa)

        subnetworks = list(k_taxon_subnetworks(network, k=k))

        # Get the set of taxa for each subnetwork
        subnet_taxa_sets = {
            frozenset(subnet.get_label(leaf) for leaf in subnet.leaves)
            for subnet in subnetworks
        }

        # All should be unique
        assert len(subnet_taxa_sets) == len(subnetworks)
        assert len(subnet_taxa_sets) == math.comb(num_taxa, k)

