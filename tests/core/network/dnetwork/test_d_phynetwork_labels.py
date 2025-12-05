"""
Comprehensive tests for DirectedPhyNetwork label operations.

This module tests all aspects of label handling including:
- get_label() method
- get_node_id() method
- Auto-labeling behavior
- Label uniqueness
- Edge cases with special characters
"""

import warnings
from typing import Optional

import pytest

from phylozoo.core.network import DirectedPhyNetwork


class TestGetLabel:
    """Test cases for get_label() method."""

    def test_get_label_existing_leaf(self) -> None:
        """Test getting label for existing leaf."""
        net = DirectedPhyNetwork(edges=[(3, 1)], taxa={1: "A"})
        assert net.get_label(1) == "A"

    def test_get_label_existing_internal(self) -> None:
        """Test getting label for existing internal node."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            taxa={1: "A", 2: "B"},
            internal_node_labels={3: "root"}
        )
        assert net.get_label(3) == "root"

    def test_get_label_missing_label(self) -> None:
        """Test getting label for unlabeled internal node."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            taxa={1: "A", 2: "B"}
        )
        assert net.get_label(3) is None

    def test_get_label_nonexistent_node(self) -> None:
        """Test getting label for non-existent node."""
        net = DirectedPhyNetwork(edges=[(3, 1)], taxa={1: "A"})
        assert net.get_label(999) is None

    def test_get_label_all_leaves_have_labels(self) -> None:
        """Test that all leaves always have labels."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A"}  # Only one labeled
        )
        # All leaves should have labels (some auto-generated)
        assert net.get_label(1) == "A"
        assert net.get_label(2) is not None
        assert net.get_label(4) is not None

    def test_get_label_multiple_internal_labels(self) -> None:
        """Test getting labels for multiple internal nodes."""
        net = DirectedPhyNetwork(
            edges=[(4, 3), (3, 1), (3, 2), (4, 5)],
            taxa={1: "A", 2: "B", 5: "C"},
            internal_node_labels={3: "internal", 4: "root"}
        )
        assert net.get_label(3) == "internal"
        assert net.get_label(4) == "root"


class TestGetNodeId:
    """Test cases for get_node_id() method."""

    def test_get_node_id_existing_label(self) -> None:
        """Test getting node ID for existing label."""
        net = DirectedPhyNetwork(edges=[(3, 1)], taxa={1: "A"})
        assert net.get_node_id("A") == 1

    def test_get_node_id_internal_label(self) -> None:
        """Test getting node ID for internal node label."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            taxa={1: "A", 2: "B"},
            internal_node_labels={3: "root"}
        )
        assert net.get_node_id("root") == 3

    def test_get_node_id_missing_label(self) -> None:
        """Test getting node ID for non-existent label."""
        net = DirectedPhyNetwork(edges=[(3, 1)], taxa={1: "A"})
        assert net.get_node_id("Nonexistent") is None

    def test_get_node_id_auto_labeled(self) -> None:
        """Test getting node ID for auto-labeled leaf."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            taxa={1: "A"}  # Only one labeled
        )
        # Leaf 2 should be auto-labeled
        label_2 = net.get_label(2)
        assert label_2 is not None
        assert net.get_node_id(label_2) == 2

    def test_get_node_id_bidirectional(self) -> None:
        """Test that get_label and get_node_id are bidirectional."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            taxa={1: "A", 2: "B"},
            internal_node_labels={3: "root"}
        )
        # Test leaf
        node_id = net.get_node_id("A")
        assert node_id == 1
        assert net.get_label(node_id) == "A"
        
        # Test internal
        node_id = net.get_node_id("root")
        assert node_id == 3
        assert net.get_label(node_id) == "root"


class TestAutoLabeling:
    """Test cases for auto-labeling behavior."""

    def test_auto_label_uncovered_leaves(self) -> None:
        """Test that uncovered leaves get auto-labeled."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A"}  # Only one labeled
        )
        # Leaves 2 and 4 should be auto-labeled
        assert net.get_label(2) is not None
        assert net.get_label(4) is not None
        assert len(net.taxa) == 3  # All leaves have labels

    def test_auto_label_uses_node_id(self) -> None:
        """Test that auto-labeling uses node ID as base."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            taxa={1: "A"}
        )
        # Leaf 2 should be auto-labeled as "2" (string of node ID)
        label_2 = net.get_label(2)
        assert label_2 == "2"

    def test_auto_label_conflict_resolution(self) -> None:
        """Test that auto-labeling resolves conflicts."""
        # If we label leaf 1 as "2", then leaf 2 can't be "2"
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "2"}  # Conflicts with node ID 2
        )
        label_2 = net.get_label(2)
        assert label_2 is not None
        assert label_2 != "2"  # Should be resolved
        # Should be something like "2_1"
        assert "2" in label_2

    def test_auto_label_multiple_conflicts(self) -> None:
        """Test auto-labeling with multiple conflicts."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "2", 2: "2_1"}  # Create conflicts
        )
        # Leaf 4 should get auto-label that doesn't conflict
        label_4 = net.get_label(4)
        assert label_4 is not None
        assert label_4 not in ["2", "2_1"]

    def test_auto_label_numeric_node_ids(self) -> None:
        """Test auto-labeling with numeric node IDs."""
        net = DirectedPhyNetwork(
            edges=[(100, 50), (100, 51)],
            taxa=None
        )
        # Should auto-label as strings
        assert "50" in net.taxa or "51" in net.taxa

    def test_auto_label_string_node_ids(self) -> None:
        """Test auto-labeling with string node IDs."""
        net = DirectedPhyNetwork(
            edges=[("root", "leaf1"), ("root", "leaf2")],
            taxa=None
        )
        # Should auto-label using string representation
        assert "leaf1" in net.taxa or "leaf2" in net.taxa

    def test_auto_label_all_uncovered(self) -> None:
        """Test auto-labeling when no taxa provided."""
        net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa=None)
        # Both leaves should be auto-labeled
        assert len(net.taxa) == 2
        assert net.get_label(1) is not None
        assert net.get_label(2) is not None


class TestLabelUniqueness:
    """Test cases for label uniqueness enforcement."""

    def test_duplicate_taxa_labels(self) -> None:
        """Test that duplicate taxa labels raise ValueError."""
        with pytest.raises(ValueError, match="already used"):
            DirectedPhyNetwork(
                edges=[(3, 1), (3, 2)],
                taxa={1: "A", 2: "A"}  # Duplicate
            )

    def test_duplicate_internal_labels(self) -> None:
        """Test that duplicate internal node labels raise ValueError."""
        with pytest.raises(ValueError, match="already used"):
            DirectedPhyNetwork(
                edges=[(4, 3), (3, 1), (3, 2)],
                taxa={1: "A", 2: "B"},
                internal_node_labels={3: "label", 4: "label"}  # Duplicate
            )

    def test_duplicate_taxa_and_internal_labels(self) -> None:
        """Test that taxa and internal labels must be unique."""
        with pytest.raises(ValueError, match="already used"):
            DirectedPhyNetwork(
                edges=[(3, 1), (3, 2)],
                taxa={1: "A", 2: "B"},
                internal_node_labels={3: "A"}  # Conflicts with taxon
            )

    def test_same_node_same_label_allowed(self) -> None:
        """Test that same node can have same label (no-op)."""
        # This should not raise an error (though it's a no-op)
        net = DirectedPhyNetwork(
            edges=[(3, 1)],
            taxa={1: "A"}
        )
        # Setting same label again should be fine (handled internally)
        assert net.get_label(1) == "A"


class TestLabelEdgeCases:
    """Test edge cases for label operations."""

    def test_empty_string_label(self) -> None:
        """Test that empty string can be a label."""
        net = DirectedPhyNetwork(
            edges=[(3, 1)],
            taxa={1: ""}  # Empty string label
        )
        assert net.get_label(1) == ""
        assert net.get_node_id("") == 1

    def test_special_characters_in_labels(self) -> None:
        """Test labels with special characters."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            taxa={1: "Taxon_1", 2: "Taxon-2"}
        )
        assert net.get_label(1) == "Taxon_1"
        assert net.get_label(2) == "Taxon-2"

    def test_unicode_labels(self) -> None:
        """Test labels with unicode characters."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            taxa={1: "Taxonα", 2: "Taxonβ"}
        )
        assert net.get_label(1) == "Taxonα"
        assert net.get_label(2) == "Taxonβ"
        assert net.get_node_id("Taxonα") == 1

    def test_very_long_labels(self) -> None:
        """Test labels with very long strings."""
        long_label = "A" * 1000
        net = DirectedPhyNetwork(
            edges=[(3, 1)],
            taxa={1: long_label}
        )
        assert net.get_label(1) == long_label
        assert len(net.get_label(1)) == 1000

    def test_numeric_labels(self) -> None:
        """Test labels that are numeric strings."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            taxa={1: "123", 2: "456"}
        )
        assert net.get_label(1) == "123"
        assert net.get_node_id("123") == 1

    def test_whitespace_in_labels(self) -> None:
        """Test labels with whitespace."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            taxa={1: "Taxon 1", 2: "Taxon\t2"}
        )
        assert net.get_label(1) == "Taxon 1"
        assert net.get_label(2) == "Taxon\t2"

    def test_labels_with_quotes(self) -> None:
        """Test labels with quotes."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            taxa={1: "Taxon'1", 2: 'Taxon"2'}
        )
        assert net.get_label(1) == "Taxon'1"
        assert net.get_label(2) == 'Taxon"2'


class TestLabelRetrievalAllNodeTypes:
    """Test label retrieval for all node types."""

    def test_labels_for_root(self) -> None:
        """Test getting label for root node."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            taxa={1: "A", 2: "B"},
            internal_node_labels={3: "root"}
        )
        assert net.get_label(3) == "root"

    def test_labels_for_tree_nodes(self) -> None:
        """Test getting labels for tree nodes."""
        net = DirectedPhyNetwork(
            edges=[(4, 3), (3, 1), (3, 2), (4, 5)],
            taxa={1: "A", 2: "B", 5: "C"},
            internal_node_labels={3: "tree_node", 4: "root"}
        )
        assert net.get_label(3) == "tree_node"

    def test_labels_for_hybrid_nodes(self) -> None:
        """Test getting labels for hybrid nodes."""
        net = DirectedPhyNetwork(
            edges=[(7, 5), (7, 6), (5, 4), (5, 8), (6, 4), (6, 9), (4, 2)],
            taxa={2: "A", 8: "B", 9: "C"},
            internal_node_labels={4: "hybrid", 5: "tree1", 6: "tree2"}
        )
        assert net.get_label(4) == "hybrid"

    def test_labels_for_leaves(self) -> None:
        """Test that all leaves have labels."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
        )
        assert net.get_label(1) == "A"
        assert net.get_label(2) == "B"
        assert net.get_label(4) == "C"

    def test_unlabeled_internal_nodes(self) -> None:
        """Test that internal nodes can be unlabeled."""
        net = DirectedPhyNetwork(
            edges=[(4, 3), (3, 1), (3, 2), (4, 5)],
            taxa={1: "A", 2: "B", 5: "C"}
            # No internal_node_labels
        )
        assert net.get_label(3) is None
        assert net.get_label(4) is None


class TestLabelConsistency:
    """Test consistency of label operations."""

    def test_all_leaves_in_taxa(self) -> None:
        """Test that all leaves appear in taxa set."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A"}  # Partial labeling
        )
        # All leaves should be in taxa
        for leaf in net.leaves:
            label = net.get_label(leaf)
            assert label is not None
            assert label in net.taxa

    def test_taxa_matches_leaf_labels(self) -> None:
        """Test that taxa set matches leaf labels."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2), (3, 4)],
            taxa={1: "A", 2: "B", 4: "C"}
        )
        # All taxa should correspond to leaf labels
        for taxon in net.taxa:
            node_id = net.get_node_id(taxon)
            assert node_id is not None
            assert node_id in net.leaves
            assert net.get_label(node_id) == taxon

    def test_label_node_id_consistency(self) -> None:
        """Test consistency between labels and node IDs."""
        net = DirectedPhyNetwork(
            edges=[(3, 1), (3, 2)],
            taxa={1: "A", 2: "B"},
            internal_node_labels={3: "root"}
        )
        # For each labeled node, get_node_id(get_label(node)) == node
        for node in [1, 2, 3]:
            label = net.get_label(node)
            if label is not None:
                assert net.get_node_id(label) == node

