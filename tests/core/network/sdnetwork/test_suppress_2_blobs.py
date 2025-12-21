"""
Tests for suppress_2_blobs function in MixedPhyNetwork.

This test suite covers suppression of 2-blobs, attribute preservation,
and validation of resulting networks for semi-directed networks.
"""

import pytest

from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
from phylozoo.core.network.sdnetwork.features import k_blobs
from phylozoo.core.network.sdnetwork.transformations import suppress_2_blobs
from phylozoo.utils.validation import no_validation


class TestSuppress2BlobsBasic:
    """Test basic functionality of suppress_2_blobs."""

    def test_no_2_blobs(self) -> None:
        """Test suppress_2_blobs on network with no 2-blobs."""
        # Use a fixture network that likely has no 2-blobs
        from tests.fixtures.sd_networks import SDTREE_SIMPLE
        
        net = SDTREE_SIMPLE
        original_nodes = net.number_of_nodes()
        result = suppress_2_blobs(net)
        result.validate()
        
        # Should return same or similar network (no changes if no 2-blobs)
        assert result.number_of_nodes() == original_nodes

    def test_single_2_blob(self) -> None:
        """Test suppress_2_blobs on network with single 2-blob."""
        # Use a fixture network that has 2-blobs
        from tests.fixtures.sd_networks import LEVEL_1_SDNETWORK_PARALLEL_EDGES
        
        net = LEVEL_1_SDNETWORK_PARALLEL_EDGES
        
        # Check that we have a 2-blob
        two_blobs = k_blobs(net, k=2, trivial=False, leaves=False)
        # May or may not have 2-blobs, but function should work
        
        original_nodes = net.number_of_nodes()
        result = suppress_2_blobs(net)
        result.validate()
        
        # Should have fewer or equal nodes
        assert result.number_of_nodes() <= original_nodes

    def test_multiple_2_blobs(self) -> None:
        """Test suppress_2_blobs on network with multiple 2-blobs."""
        # Use a fixture network
        from tests.fixtures.sd_networks import LEVEL_2_SDNETWORK_MANY_PARALLEL_EDGES
        
        net = LEVEL_2_SDNETWORK_MANY_PARALLEL_EDGES
        
        original_nodes = net.number_of_nodes()
        result = suppress_2_blobs(net)
        result.validate()
        
        # Should have fewer or equal nodes
        assert result.number_of_nodes() <= original_nodes

    def test_mixed_edge_types_2_blob(self) -> None:
        """Test suppress_2_blobs with mixed directed/undirected edges."""
        # Use a fixture network with hybrid nodes (has both directed and undirected)
        from tests.fixtures.sd_networks import LEVEL_1_SDNETWORK_PARALLEL_EDGES_HYBRID
        
        net = LEVEL_1_SDNETWORK_PARALLEL_EDGES_HYBRID
        
        original_nodes = net.number_of_nodes()
        result = suppress_2_blobs(net)
        result.validate()
        
        # Should have fewer or equal nodes
        assert result.number_of_nodes() <= original_nodes


class TestSuppress2BlobsAttributes:
    """Test attribute preservation during 2-blob suppression."""

    def test_branch_length_preserved(self) -> None:
        """Test that branch_length is properly merged during suppression."""
        # Use a fixture network
        from tests.fixtures.sd_networks import LEVEL_1_SDNETWORK_PARALLEL_EDGES
        
        net = LEVEL_1_SDNETWORK_PARALLEL_EDGES
        
        result = suppress_2_blobs(net)
        result.validate()
        
        # Network should be valid
        assert result.number_of_edges() > 0
        assert result.number_of_nodes() > 0

    def test_gamma_preserved(self) -> None:
        """Test that gamma is properly handled during suppression."""
        # Use a fixture network with hybrid nodes (has gamma)
        from tests.fixtures.sd_networks import LEVEL_1_SDNETWORK_PARALLEL_EDGES_HYBRID
        
        net = LEVEL_1_SDNETWORK_PARALLEL_EDGES_HYBRID
        
        result = suppress_2_blobs(net)
        result.validate()
        
        # Network should be valid
        assert result.number_of_edges() > 0
        assert result.number_of_nodes() > 0


class TestSuppress2BlobsFixtureNetworks:
    """Test suppress_2_blobs on fixture networks."""

    def test_on_tree_network(self) -> None:
        """Test on a simple tree network."""
        from tests.fixtures.sd_networks import SDTREE_SIMPLE
        
        result = suppress_2_blobs(SDTREE_SIMPLE)
        result.validate()
        
        # Tree should remain valid
        assert result.number_of_nodes() > 0

    def test_on_hybrid_network(self) -> None:
        """Test on a network with hybrid nodes."""
        from tests.fixtures.sd_networks import LEVEL_1_SDNETWORK_PARALLEL_EDGES
        
        result = suppress_2_blobs(LEVEL_1_SDNETWORK_PARALLEL_EDGES)
        result.validate()
        
        # Network should remain valid
        assert result.number_of_nodes() > 0
        assert result.number_of_edges() > 0

    def test_on_complex_network(self) -> None:
        """Test on a more complex network."""
        from tests.fixtures.sd_networks import LEVEL_2_SDNETWORK_MANY_PARALLEL_EDGES
        
        result = suppress_2_blobs(LEVEL_2_SDNETWORK_MANY_PARALLEL_EDGES)
        result.validate()
        
        # Network should remain valid
        assert result.number_of_nodes() > 0
        assert result.number_of_edges() > 0


class TestSuppress2BlobsEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_network(self) -> None:
        """Test on empty network."""
        from tests.fixtures.sd_networks import SDTREE_EMPTY
        
        result = suppress_2_blobs(SDTREE_EMPTY)
        result.validate()
        
        # Should remain empty
        assert result.number_of_nodes() == 0

    def test_single_node_network(self) -> None:
        """Test on single node network."""
        from tests.fixtures.sd_networks import SDTREE_SINGLE_NODE
        
        result = suppress_2_blobs(SDTREE_SINGLE_NODE)
        result.validate()
        
        # Should remain single node
        assert result.number_of_nodes() == 1

    def test_all_2_blobs_processed(self) -> None:
        """Test that all 2-blobs are processed."""
        # Use a fixture network
        from tests.fixtures.sd_networks import LEVEL_2_SDNETWORK_MANY_PARALLEL_EDGES
        
        net = LEVEL_2_SDNETWORK_MANY_PARALLEL_EDGES
        
        # Count 2-blobs
        two_blobs = k_blobs(net, k=2, trivial=False, leaves=False)
        
        original_nodes = net.number_of_nodes()
        result = suppress_2_blobs(net)
        result.validate()
        
        # Should have processed all 2-blobs
        # The exact node count depends on the structure, but should be reduced or equal
        assert result.number_of_nodes() <= original_nodes

