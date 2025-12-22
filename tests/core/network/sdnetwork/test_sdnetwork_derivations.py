"""
Tests for sdnetwork derivations functions.
"""

import pytest

from phylozoo.core.network.sdnetwork import MixedPhyNetwork, SemiDirectedPhyNetwork
from phylozoo.core.network.sdnetwork.derivations import tree_of_blobs
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
