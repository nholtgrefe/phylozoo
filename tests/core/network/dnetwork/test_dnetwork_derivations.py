"""
Tests for dnetwork derivations functions.
"""

import pytest

from phylozoo.core.network.dnetwork import DirectedPhyNetwork
from phylozoo.core.network.dnetwork.derivations import tree_of_blobs
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
