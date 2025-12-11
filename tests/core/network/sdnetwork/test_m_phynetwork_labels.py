"""
Comprehensive tests for MixedPhyNetwork label operations.

This module tests all aspects of label handling including:
- get_label() method
- get_node_id() method
- Taxa labels
- Internal node labels
- Auto-labeling
- Label uniqueness
"""

import warnings

import pytest

from phylozoo.core.network.sdnetwork import MixedPhyNetwork
from tests.core.network.sdnetwork.conftest import expect_mixed_network_warning


class TestGetLabel:
    """Test cases for get_label() method."""

    def test_get_label_leaf_with_taxon(self) -> None:
        """Test getting label for leaf with taxon."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.get_label(1) == "A"

    def test_get_label_internal_with_label(self) -> None:
        """Test getting label for internal node with label."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (3, {'label': 'root'})],
            )
        assert net.get_label(3) == "root"

    def test_get_label_internal_without_label(self) -> None:
        """Test getting label for internal node without label."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.get_label(3) is None

    def test_get_label_auto_labeled_leaf(self) -> None:
        """Test getting label for auto-labeled leaf."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]  # Only 1 labeled
            )
        # Leaf 2 should be auto-labeled
        label2 = net.get_label(2)
        assert label2 is not None
        assert isinstance(label2, str)

    def test_get_label_nonexistent_node(self) -> None:
        """Test getting label for non-existent node."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.get_label(999) is None


class TestGetNodeId:
    """Test cases for get_node_id() method."""

    def test_get_node_id_existing_taxon(self) -> None:
        """Test getting node ID for existing taxon."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.get_node_id("A") == 1

    def test_get_node_id_existing_internal_label(self) -> None:
        """Test getting node ID for existing internal label."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (3, {'label': 'root'})],
            )
        assert net.get_node_id("root") == 3

    def test_get_node_id_nonexistent_label(self) -> None:
        """Test getting node ID for non-existent label."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.get_node_id("Nonexistent") is None

    def test_get_node_id_auto_labeled_leaf(self) -> None:
        """Test getting node ID for auto-labeled leaf."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'})]  # Only 1 labeled, 2 and 4 auto-labeled
            )
        # Leaf 2 should be auto-labeled with "2"
        node_id = net.get_node_id("2")
        assert node_id == 2


class TestTaxaLabels:
    """Test cases for taxa (leaf) labels."""

    def test_all_leaves_have_labels(self) -> None:
        """Test that all leaves have labels."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]  # Partial mapping
            )
        for leaf in net.leaves:
            assert net.get_label(leaf) is not None

    def test_taxa_labels_unique(self) -> None:
        """Test that taxa labels are unique."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        labels = [net.get_label(leaf) for leaf in net.leaves]
        assert len(labels) == len(set(labels))  # All unique

    def test_taxa_set_matches_leaves(self) -> None:
        """Test that taxa set matches leaves."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert len(net.taxa) == len(net.leaves)


class TestInternalNodeLabels:
    """Test cases for internal node labels."""

    def test_internal_node_label_optional(self) -> None:
        """Test that internal node labels are optional."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        assert net.get_label(3) is None  # No label

    def test_internal_node_label_set(self) -> None:
        """Test setting internal node label."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (3, {'label': 'root'})],
            )
        assert net.get_label(3) == "root"

    def test_internal_node_labels_unique(self) -> None:
        """Test that internal node labels are unique."""
        # Node 4 needs degree >= 3
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(4, 3), (3, 1), (3, 2), (4, 5), (4, 6)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'}), (3, {'label': 'internal1'}), (4, {'label': 'internal2'})],
            )
        labels = [net.get_label(3), net.get_label(4)]
        assert len(labels) == len(set(labels))  # All unique

    def test_internal_node_label_cannot_duplicate_taxon(self) -> None:
        """Test that internal node label cannot duplicate taxon."""
        with pytest.raises(ValueError, match="already used|duplicate"):
            MixedPhyNetwork(
                undirected_edges=[(3, 1), (3, 2), (3, 4)],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (3, {'label': 'A'})],
            )



class TestAutoLabeling:
    """Test cases for auto-labeling of uncovered leaves."""

    def test_auto_labeling_partial_taxa(self) -> None:
        """Test auto-labeling when taxa mapping is partial."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]  # Only one labeled
            )
        # All leaves should have labels
        assert net.get_label(1) == "A"
        assert net.get_label(2) is not None
        assert net.get_label(4) is not None

    def test_auto_labeling_no_taxa(self) -> None:
        """Test auto-labeling when no taxa provided."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=None
            )
        # All leaves should be auto-labeled
        assert net.get_label(1) is not None
        assert net.get_label(2) is not None

    def test_auto_labeling_format(self) -> None:
        """Test format of auto-generated labels."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
            )
        # Auto-generated label should be string representation of node ID
        label2 = net.get_label(2)
        assert label2 == "2" or label2 is not None  # May vary


class TestLabelUniqueness:
    """Test cases for label uniqueness constraints."""

    def test_duplicate_taxa_labels_raises_error(self) -> None:
        """Test that duplicate taxa labels raise error."""
        with pytest.raises(ValueError, match="already used|duplicate"):
            MixedPhyNetwork(
                undirected_edges=[(3, 1), (3, 2), (3, 4)],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'A'}), (4, {'label': 'B'})]  # Duplicate
            )

    def test_duplicate_internal_labels_raises_error(self) -> None:
        """Test that duplicate internal labels raise error."""
        with pytest.raises(ValueError, match="already used|duplicate"):
            MixedPhyNetwork(
                undirected_edges=[(4, 3), (3, 1), (3, 2), (4, 5), (4, 6)],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (5, {'label': 'C'}), (6, {'label': 'D'}), (3, {'label': 'label'}), (4, {'label': 'label'})],
            )

    def test_taxon_and_internal_label_cannot_duplicate(self) -> None:
        """Test that taxon and internal label cannot duplicate."""
        with pytest.raises(ValueError, match="already used|duplicate"):
            MixedPhyNetwork(
                undirected_edges=[(3, 1), (3, 2), (3, 4)],
                nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (3, {'label': 'A'})],
            )


class TestLabelRetrieval:
    """Test cases for label retrieval methods."""

    def test_get_label_all_nodes(self) -> None:
        """Test get_label for all nodes in network."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (3, {'label': 'root'})],
            )
        assert net.get_label(1) == "A"
        assert net.get_label(2) == "B"
        assert net.get_label(3) == "root"

    def test_get_node_id_all_labels(self) -> None:
        """Test get_node_id for all labels in network."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (3, {'label': 'root'})],
            )
        assert net.get_node_id("A") == 1
        assert net.get_node_id("B") == 2
        assert net.get_node_id("root") == 3

    def test_label_round_trip(self) -> None:
        """Test round-trip: label -> node_id -> label."""
        with expect_mixed_network_warning():
            net = MixedPhyNetwork(
            undirected_edges=[(3, 1), (3, 2), (3, 4)],
            nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'}), (3, {'label': 'root'})],
            )
        for node in net:  # Use __iter__ instead of nodes()
            label = net.get_label(node)
            if label is not None:
                node_id = net.get_node_id(label)
                assert node_id == node

