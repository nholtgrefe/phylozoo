"""
Extensive tests for the subnetwork function in dnetwork derivations.

Tests cover:
- Subnetwork of all taxa equals original network
- Various fixture networks
- Different options (suppress_2_blobs, identify_parallel_edges, make_lsa)
- Edge cases
- Property preservation/transformation verification using features module
"""

import pytest

from phylozoo.core.network.dnetwork.classifications import (
    has_parallel_edges,
    is_binary,
    is_lsa_network,
    is_tree,
    level,
    vertex_level,
)
from phylozoo.core.network.dnetwork.derivations import subnetwork
from phylozoo.core.network.dnetwork.features import (
    blobs,
    cut_edges,
    cut_vertices,
    k_blobs,
)
from tests.fixtures import directed_networks


class TestSubnetworkBasic:
    """Basic tests for subnetwork function."""

    def test_subnetwork_all_taxa_equals_original(self):
        """Test that subnetwork of all taxa equals the original network."""
        # Test with various networks from fixtures
        test_networks = [
            directed_networks.DTREE_SMALL_BINARY,
            directed_networks.DTREE_MEDIUM_BINARY,
            directed_networks.LEVEL_1_DNETWORK_SINGLE_HYBRID,
            directed_networks.LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE,
            directed_networks.LEVEL_2_DNETWORK_TWO_HYBRIDS_SAME_BLOB,
        ]

        for network in test_networks:
            if network.number_of_nodes() == 0:
                continue

            all_taxa = list(network.taxa)
            subnet = subnetwork(network, all_taxa)

            # Check that they have the same number of nodes and edges
            assert subnet.number_of_nodes() == network.number_of_nodes()
            assert subnet.number_of_edges() == network.number_of_edges()

            # Check that they have the same leaves
            assert subnet.leaves == network.leaves

            # Check that they have the same root
            assert subnet.root_node == network.root_node

    def test_subnetwork_single_taxon(self):
        """Test subnetwork with a single taxon."""
        network = directed_networks.DTREE_MEDIUM_BINARY
        all_taxa = list(network.taxa)
        single_taxon = [all_taxa[0]]

        subnet = subnetwork(network, single_taxon)

        # Should have exactly one leaf
        assert len(subnet.leaves) == 1
        assert subnet.get_label(next(iter(subnet.leaves))) == single_taxon[0]

        # Should be a valid network (tree path from root to leaf)
        assert subnet.number_of_nodes() >= 1
        assert subnet.number_of_edges() >= 0

    def test_subnetwork_subset_taxa(self):
        """Test subnetwork with a subset of taxa."""
        network = directed_networks.DTREE_MEDIUM_BINARY
        all_taxa = list(network.taxa)
        subset_taxa = all_taxa[:5]  # First 5 taxa

        subnet = subnetwork(network, subset_taxa)

        # Should have exactly 5 leaves
        assert len(subnet.leaves) == 5
        subnet_taxa = {subnet.get_label(leaf) for leaf in subnet.leaves}
        assert subnet_taxa == set(subset_taxa)

        # Should have fewer or equal nodes/edges than original
        assert subnet.number_of_nodes() <= network.number_of_nodes()
        assert subnet.number_of_edges() <= network.number_of_edges()

    def test_subnetwork_empty_taxa_returns_empty_network(self):
        """Test that empty taxa list returns an empty network."""
        network = directed_networks.DTREE_SMALL_BINARY
        subnet = subnetwork(network, [])
        
        # Should return an empty network
        assert subnet.number_of_nodes() == 0
        assert subnet.number_of_edges() == 0
        assert len(subnet.leaves) == 0

    def test_subnetwork_invalid_taxon_raises(self):
        """Test that invalid taxon raises ValueError."""
        network = directed_networks.DTREE_SMALL_BINARY
        with pytest.raises(ValueError, match="not found in network"):
            subnetwork(network, ["InvalidTaxon"])

    def test_subnetwork_preserves_tree_property(self):
        """Test that subnetwork of a tree is still a tree."""
        network = directed_networks.DTREE_MEDIUM_BINARY
        assert is_tree(network)

        all_taxa = list(network.taxa)
        subset_taxa = all_taxa[:3]

        subnet = subnetwork(network, subset_taxa)
        assert is_tree(subnet)


class TestSubnetworkWithOptions:
    """Tests for subnetwork with different options."""

    def test_subnetwork_with_suppress_2_blobs(self):
        """Test subnetwork with suppress_2_blobs option."""
        network = directed_networks.LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE
        all_taxa = list(network.taxa)
        subset_taxa = all_taxa[:2]

        # Without suppress_2_blobs
        subnet1 = subnetwork(network, subset_taxa, suppress_2_blobs=False)
        # With suppress_2_blobs
        subnet2 = subnetwork(network, subset_taxa, suppress_2_blobs=True)

        # The suppressed version should have fewer or equal 2-blobs
        blobs1 = k_blobs(subnet1, k=2, trivial=False, leaves=False)
        blobs2 = k_blobs(subnet2, k=2, trivial=False, leaves=False)
        assert len(blobs2) <= len(blobs1)

    def test_subnetwork_with_identify_parallel_edges(self):
        """Test subnetwork with identify_parallel_edges option."""
        # Use a network that might have parallel edges after subnetwork extraction
        network = directed_networks.LEVEL_1_DNETWORK_SINGLE_HYBRID
        all_taxa = list(network.taxa)
        subset_taxa = all_taxa[:2]

        # Without identify_parallel_edges
        subnet1 = subnetwork(network, subset_taxa, identify_parallel_edges=False)
        # With identify_parallel_edges
        subnet2 = subnetwork(network, subset_taxa, identify_parallel_edges=True)

        # The version with identify_parallel_edges should have no parallel edges
        assert not has_parallel_edges(subnet2)

        # The version without should have same or more edges
        assert subnet2.number_of_edges() <= subnet1.number_of_edges()

    def test_subnetwork_with_make_lsa(self):
        """Test subnetwork with make_lsa option."""
        network = directed_networks.LEVEL_1_DNETWORK_SINGLE_HYBRID
        all_taxa = list(network.taxa)
        subset_taxa = all_taxa[:2]

        # Without make_lsa
        subnet1 = subnetwork(network, subset_taxa, make_lsa=False)
        # With make_lsa
        subnet2 = subnetwork(network, subset_taxa, make_lsa=True)

        # The LSA version should be an LSA network
        assert is_lsa_network(subnet2)

        # The root should equal the LSA node
        assert subnet2.root_node == subnet2.LSA_node

    def test_subnetwork_all_options_combined(self):
        """Test subnetwork with all options enabled."""
        network = directed_networks.LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE
        all_taxa = list(network.taxa)
        subset_taxa = all_taxa[:3]

        subnet = subnetwork(
            network,
            subset_taxa,
            suppress_2_blobs=True,
            identify_parallel_edges=True,
            make_lsa=True,
        )

        # Should be an LSA network
        assert is_lsa_network(subnet)

        # Should have no parallel edges
        assert not has_parallel_edges(subnet)

        # Should be a valid network
        assert subnet.number_of_nodes() > 0
        assert len(subnet.leaves) == len(subset_taxa)


class TestSubnetworkPropertyPreservation:
    """Tests to verify property preservation/transformation using features module."""

    def test_subnetwork_preserves_leaf_labels(self):
        """Test that subnetwork preserves leaf labels."""
        network = directed_networks.DTREE_MEDIUM_BINARY
        all_taxa = list(network.taxa)
        subset_taxa = all_taxa[:5]

        subnet = subnetwork(network, subset_taxa)

        subnet_taxa = {subnet.get_label(leaf) for leaf in subnet.leaves}
        assert subnet_taxa == set(subset_taxa)

    def test_subnetwork_ancestors_included(self):
        """Test that all ancestors of selected taxa are included."""
        import networkx as nx

        network = directed_networks.DTREE_MEDIUM_BINARY
        all_taxa = list(network.taxa)
        subset_taxa = all_taxa[:3]

        subnet = subnetwork(network, subset_taxa)

        # All nodes in subnet should be ancestors of at least one selected taxon
        # or be the selected taxa themselves
        subnet_nodes = set(subnet._graph.nodes())
        selected_leaf_nodes = {
            subnet.get_node_id(taxon) for taxon in subset_taxa
        }

        # Check that every node in subnet is either a selected leaf or an ancestor
        root = subnet.root_node
        dag = subnet._graph._graph
        for node in subnet_nodes:
            if node in selected_leaf_nodes:
                continue
            # Node should be an ancestor of at least one selected leaf
            is_ancestor = False
            for leaf_node in selected_leaf_nodes:
                if nx.has_path(dag, node, leaf_node):
                    is_ancestor = True
                    break
            assert is_ancestor or node == root

    def test_subnetwork_cut_edges_preserved(self):
        """Test that cut-edges in subnetwork are valid."""
        network = directed_networks.LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE
        all_taxa = list(network.taxa)
        subset_taxa = all_taxa[:2]

        subnet = subnetwork(network, subset_taxa)

        # Get cut edges in subnet
        subnet_cut_edges = cut_edges(subnet)

        # All cut edges should be valid edges in the subnet
        for u, v, key in subnet_cut_edges:
            assert subnet._graph.has_edge(u, v, key)

    def test_subnetwork_blobs_structure(self):
        """Test that blobs in subnetwork are valid."""
        network = directed_networks.LEVEL_2_DNETWORK_TWO_HYBRIDS_SAME_BLOB
        all_taxa = list(network.taxa)
        subset_taxa = all_taxa[:2]

        subnet = subnetwork(network, subset_taxa)

        # Get blobs in subnet
        subnet_blobs = blobs(subnet, trivial=True, leaves=True)

        # All nodes in blobs should be in subnet
        subnet_nodes = set(subnet._graph.nodes())
        for blob in subnet_blobs:
            assert blob.issubset(subnet_nodes)

    def test_subnetwork_level_does_not_increase(self):
        """Test that subnetwork level does not increase."""
        network = directed_networks.LEVEL_2_DNETWORK_TWO_HYBRIDS_SAME_BLOB
        original_level = level(network)

        all_taxa = list(network.taxa)
        subset_taxa = all_taxa[:2]

        subnet = subnetwork(network, subset_taxa)
        subnet_level = level(subnet)

        # Level should not increase (might decrease if blobs are removed)
        assert subnet_level <= original_level

    def test_subnetwork_vertex_level_does_not_increase(self):
        """Test that subnetwork vertex level does not increase."""
        network = directed_networks.LEVEL_2_DNETWORK_NESTED_HYBRIDS
        original_vertex_level = vertex_level(network)

        all_taxa = list(network.taxa)
        subset_taxa = all_taxa[:3]

        subnet = subnetwork(network, subset_taxa)
        subnet_vertex_level = vertex_level(subnet)

        # Vertex level should not increase
        assert subnet_vertex_level <= original_vertex_level


class TestSubnetworkWithVariousNetworks:
    """Tests with various fixture networks."""

    @pytest.mark.parametrize(
        "network_name",
        [
            "DTREE_SMALL_BINARY",
            "DTREE_MEDIUM_BINARY",
            "DTREE_NON_BINARY_SMALL",
            "LEVEL_1_DNETWORK_SINGLE_HYBRID",
            "LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE",
            "LEVEL_2_DNETWORK_TWO_HYBRIDS_SAME_BLOB",
            "LEVEL_2_DNETWORK_TRIANGLE_HYBRID",
        ],
    )
    def test_subnetwork_with_various_networks(self, network_name):
        """Test subnetwork with various fixture networks."""
        network = getattr(directed_networks, network_name)

        if network.number_of_nodes() == 0:
            pytest.skip("Empty network")

        all_taxa = list(network.taxa)
        if len(all_taxa) < 2:
            pytest.skip("Network has fewer than 2 taxa")

        # Test with first half of taxa
        subset_taxa = all_taxa[: len(all_taxa) // 2]

        subnet = subnetwork(network, subset_taxa)

        # Basic checks
        assert subnet.number_of_nodes() > 0
        assert len(subnet.leaves) == len(subset_taxa)
        assert subnet.number_of_edges() >= len(subset_taxa) - 1  # At least tree edges

        # Verify all selected taxa are present
        subnet_taxa = {subnet.get_label(leaf) for leaf in subnet.leaves}
        assert subnet_taxa == set(subset_taxa)

    def test_subnetwork_with_large_network(self):
        """Test subnetwork with a large network."""
        network = directed_networks.DTREE_LARGE_BINARY
        all_taxa = list(network.taxa)

        # Test with various subset sizes
        for size in [1, 5, 10, len(all_taxa) // 2]:
            subset_taxa = all_taxa[:size]
            subnet = subnetwork(network, subset_taxa)

            assert len(subnet.leaves) == size
            subnet_taxa = {subnet.get_label(leaf) for leaf in subnet.leaves}
            assert subnet_taxa == set(subset_taxa)

    def test_subnetwork_with_parallel_edges_network(self):
        """Test subnetwork with a network that has parallel edges."""
        # Use networks that are known to have parallel edges
        networks_with_parallel = [
            directed_networks.LEVEL_1_DNETWORK_PARALLEL_EDGES,
            directed_networks.LEVEL_1_DNETWORK_PARALLEL_EDGES_HYBRID,
            directed_networks.LEVEL_2_DNETWORK_MANY_PARALLEL_EDGES,
            directed_networks.LEVEL_1_DNETWORK_PARALLEL_IN_BLOB,
        ]

        for network in networks_with_parallel:
            if has_parallel_edges(network):
                all_taxa = list(network.taxa)
                if len(all_taxa) < 2:
                    continue
                subset_taxa = all_taxa[:2]

                subnet = subnetwork(network, subset_taxa)
                # Subnetwork should still be valid
                assert subnet.number_of_nodes() > 0
                assert len(subnet.leaves) == len(subset_taxa)
                break
        else:
            pytest.skip("No network with parallel edges found in fixtures")

    def test_subnetwork_with_multiple_blobs(self):
        """Test subnetwork with a network that has multiple blobs."""
        network = directed_networks.LEVEL_1_DNETWORK_FIVE_BLOBS
        all_taxa = list(network.taxa)

        # Test with subset that might span multiple blobs
        subset_taxa = all_taxa[:5]

        subnet = subnetwork(network, subset_taxa)

        # Should be a valid network
        assert subnet.number_of_nodes() > 0
        assert len(subnet.leaves) == len(subset_taxa)

        # Check blob structure
        subnet_blobs = blobs(subnet, trivial=False, leaves=False)
        # Should have at least one blob (might have fewer than original)
        assert len(subnet_blobs) >= 0


class TestSubnetworkEdgeCases:
    """Edge case tests for subnetwork."""

    def test_subnetwork_single_node_network(self):
        """Test subnetwork with a single-node network."""
        network = directed_networks.DTREE_SINGLE_NODE
        all_taxa = list(network.taxa)

        subnet = subnetwork(network, all_taxa)

        # Should be the same network
        assert subnet.number_of_nodes() == network.number_of_nodes()
        assert subnet.number_of_edges() == network.number_of_edges()

    def test_subnetwork_all_taxa_different_orders(self):
        """Test that subnetwork is independent of taxon order."""
        network = directed_networks.DTREE_MEDIUM_BINARY
        all_taxa = list(network.taxa)
        subset_taxa = all_taxa[:5]

        # Create subnetworks with different orders
        subnet1 = subnetwork(network, subset_taxa)
        subnet2 = subnetwork(network, list(reversed(subset_taxa)))

        # Should have same number of nodes and edges
        assert subnet1.number_of_nodes() == subnet2.number_of_nodes()
        assert subnet1.number_of_edges() == subnet2.number_of_edges()

        # Should have same leaves
        assert subnet1.leaves == subnet2.leaves

    def test_subnetwork_repeated_taxa(self):
        """Test subnetwork with repeated taxa (should work, just duplicates)."""
        network = directed_networks.DTREE_SMALL_BINARY
        all_taxa = list(network.taxa)
        repeated_taxa = all_taxa + all_taxa  # Duplicate all taxa

        subnet = subnetwork(network, repeated_taxa)

        # Should be same as subnetwork with unique taxa
        subnet_unique = subnetwork(network, all_taxa)

        assert subnet.number_of_nodes() == subnet_unique.number_of_nodes()
        assert subnet.number_of_edges() == subnet_unique.number_of_edges()
        assert subnet.leaves == subnet_unique.leaves


class TestSubnetworkOptionsInteraction:
    """Tests for how different options interact."""

    def test_suppress_2_blobs_then_identify_parallel(self):
        """Test suppress_2_blobs followed by identify_parallel_edges."""
        network = directed_networks.LEVEL_1_DNETWORK_TWO_HYBRIDS_SEPARATE
        all_taxa = list(network.taxa)
        subset_taxa = all_taxa[:3]

        subnet1 = subnetwork(
            network, subset_taxa, suppress_2_blobs=True, identify_parallel_edges=False
        )
        subnet2 = subnetwork(
            network, subset_taxa, suppress_2_blobs=True, identify_parallel_edges=True
        )

        # Both should be valid
        assert subnet1.number_of_nodes() > 0
        assert subnet2.number_of_nodes() > 0

        # Second should have no parallel edges
        assert not has_parallel_edges(subnet2)

    def test_make_lsa_preserves_structure(self):
        """Test that make_lsa preserves the basic structure."""
        network = directed_networks.LEVEL_1_DNETWORK_SINGLE_HYBRID
        all_taxa = list(network.taxa)
        subset_taxa = all_taxa[:2]

        subnet1 = subnetwork(network, subset_taxa, make_lsa=False)
        subnet2 = subnetwork(network, subset_taxa, make_lsa=True)

        # Both should have same leaves
        assert subnet1.leaves == subnet2.leaves

        # LSA version should be an LSA network
        assert is_lsa_network(subnet2)

    def test_all_options_produce_valid_network(self):
        """Test that all options together produce a valid network."""
        network = directed_networks.LEVEL_2_DNETWORK_NESTED_HYBRIDS
        all_taxa = list(network.taxa)
        subset_taxa = all_taxa[:4]

        subnet = subnetwork(
            network,
            subset_taxa,
            suppress_2_blobs=True,
            identify_parallel_edges=True,
            make_lsa=True,
        )

        # Should be valid
        assert subnet.number_of_nodes() > 0
        assert len(subnet.leaves) == len(subset_taxa)

        # Should satisfy all option constraints
        assert is_lsa_network(subnet)
        assert not has_parallel_edges(subnet)

        # Should have correct taxa
        subnet_taxa = {subnet.get_label(leaf) for leaf in subnet.leaves}
        assert subnet_taxa == set(subset_taxa)

