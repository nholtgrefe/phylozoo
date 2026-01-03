"""
Tests for the split module.
"""

import pytest
from phylozoo.core.split import Split, SplitSystem, is_compatible, is_subsplit


class TestSplit:
    """Test cases for the Split class."""

    def test_split_creation(self) -> None:
        """Test creating a valid split."""
        split = Split({1, 2}, {3, 4})
        assert split.set1 == {1, 2}
        assert split.set2 == {3, 4}
        assert split.elements == {1, 2, 3, 4}

    def test_split_empty_set1_raises_error(self) -> None:
        """Test that creating a split with empty set1 raises ValueError."""
        with pytest.raises(ValueError, match="Split sets cannot be empty"):
            Split(set(), {3, 4})

    def test_split_empty_set2_raises_error(self) -> None:
        """Test that creating a split with empty set2 raises ValueError."""
        with pytest.raises(ValueError, match="Split sets cannot be empty"):
            Split({1, 2}, set())

    def test_split_canonical_ordering(self) -> None:
        """Test that splits use canonical ordering (set1 is always the same)."""
        # Create splits with different input orders
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({3, 4}, {1, 2})
        
        # Both should have the same canonical ordering
        assert split1.set1 == split2.set1
        assert split1.set2 == split2.set2
        assert split1 == split2

    def test_split_canonical_ordering_size_based(self) -> None:
        """Test that canonical ordering prioritizes smaller sets first."""
        # Smaller set should be set1
        split = Split({1, 2, 3}, {4})
        # Since {4} is smaller, it should be set1 in canonical form
        assert len(split.set1) <= len(split.set2)

    def test_is_trivial(self) -> None:
        """Test the is_trivial method."""
        trivial_split = Split({1}, {2, 3, 4})
        assert trivial_split.is_trivial() is True

        non_trivial_split = Split({1, 2}, {3, 4})
        assert non_trivial_split.is_trivial() is False

    def test_is_compatible_true_same_elements(self) -> None:
        """Test is_compatible returns True for compatible splits with same elements."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1}, {2, 3, 4})
        assert is_compatible(split1, split2) is True

    def test_is_compatible_true_other_way(self) -> None:
        """Test is_compatible returns True when other split has subset."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 2, 3}, {4})
        assert is_compatible(split1, split2) is True

    def test_is_compatible_false_different_elements(self) -> None:
        """Test is_compatible returns False for splits with different elements."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 2}, {3, 5})
        assert is_compatible(split1, split2) is False

    def test_is_compatible_false_same_elements_incompatible(self) -> None:
        """Test is_compatible returns False for incompatible splits with same elements."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        assert is_compatible(split1, split2) is False

    def test_is_compatible_raises_error_first_arg(self) -> None:
        """Test that is_compatible raises error for non-Split as first argument."""
        split = Split({1, 2}, {3, 4})
        with pytest.raises(ValueError, match="First argument must be a Split instance"):
            is_compatible("not a split", split)  # type: ignore

    def test_is_compatible_raises_error_second_arg(self) -> None:
        """Test that is_compatible raises error for non-Split as second argument."""
        split = Split({1, 2}, {3, 4})
        with pytest.raises(ValueError, match="Second argument must be a Split instance"):
            is_compatible(split, "not a split")  # type: ignore

    def test_is_subsplit_true(self) -> None:
        """Test is_subsplit returns True when one split is subsplit of another."""
        split1 = Split({1, 2, 6}, {3, 4, 5})
        split2 = Split({1, 2}, {3, 4})
        assert is_subsplit(split2, split1) is True

    def test_is_subsplit_false(self) -> None:
        """Test is_subsplit returns False when not a subsplit."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        assert is_subsplit(split1, split2) is False

    def test_is_subsplit_raises_error_first_arg(self) -> None:
        """Test that is_subsplit raises error for non-Split as first argument."""
        split = Split({1, 2}, {3, 4})
        with pytest.raises(ValueError, match="First argument must be a Split instance"):
            is_subsplit("not a split", split)  # type: ignore

    def test_is_subsplit_raises_error_second_arg(self) -> None:
        """Test that is_subsplit raises error for non-Split as second argument."""
        split = Split({1, 2}, {3, 4})
        with pytest.raises(ValueError, match="Second argument must be a Split instance"):
            is_subsplit(split, "not a split")  # type: ignore

    def test_induced_quartetsplits(self) -> None:
        """Test induced_quartetsplits generates correct quartet splits."""
        split = Split({1, 2, 3}, {4, 5, 6})
        quartets = split.induced_quartetsplits()
        
        # Should generate C(3,2) * C(3,2) = 3 * 3 = 9 quartet splits
        assert len(quartets) == 9
        
        # All should be 2|2 splits
        for q in quartets:
            assert len(q.set1) == 2
            assert len(q.set2) == 2

    def test_induced_quartetsplits_include_trivial(self) -> None:
        """Test induced_quartetsplits with include_trivial=True."""
        split = Split({1, 2}, {3, 4, 5})
        quartets = split.induced_quartetsplits(include_trivial=True)
        
        # Should have 2|2 splits: C(2,2) * C(3,2) = 1 * 3 = 3
        # Plus trivial 1|3 splits: C(2,1) * C(3,3) + C(2,3) * C(3,1) = 2 * 1 + 0 * 3 = 2
        # Total: 3 + 2 = 5
        assert len(quartets) >= 3  # At least the 2|2 splits
        
        # Check that trivial splits are included
        trivial_count = sum(1 for q in quartets if q.is_trivial())
        assert trivial_count > 0

    def test_split_repr(self) -> None:
        """Test string representation of a split."""
        split = Split({1, 2}, {3, 4})
        repr_str = repr(split)
        assert "Split" in repr_str

    def test_split_equality(self) -> None:
        """Test that splits with same sets are equal regardless of input order."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({3, 4}, {1, 2})
        assert split1 == split2

    def test_split_hash(self) -> None:
        """Test that equal splits have the same hash."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({3, 4}, {1, 2})
        assert hash(split1) == hash(split2)

    def test_split_overlapping_sets_raises_error(self) -> None:
        """Test that overlapping sets raise ValueError."""
        with pytest.raises(ValueError, match="Invalid partition: sets overlap"):
            Split({1, 2}, {2, 3})


class TestSplitSystem:
    """Test cases for the SplitSystem class."""

    def test_split_system_creation(self) -> None:
        """Test creating a split system."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        system = SplitSystem([split1, split2])
        assert len(system) == 2

    def test_split_system_with_single_split(self) -> None:
        """Test creating a split system with a single split."""
        split = Split({1, 2}, {3, 4})
        system = SplitSystem([split])
        assert len(system) == 1

    def test_split_system_repr(self) -> None:
        """Test string representation of a split system."""
        system = SplitSystem()
        repr_str = repr(system)
        assert "SplitSystem" in repr_str


class TestSplitSystemImmutability:
    """Test cases for SplitSystem immutability."""

    def test_split_system_immutability(self) -> None:
        """Test that SplitSystem is immutable."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        system = SplitSystem([split1, split2])
        
        # Try to modify splits attribute
        with pytest.raises(AttributeError, match="Cannot modify attribute"):
            system.splits = set()
        
        # Try to modify elements attribute
        with pytest.raises(AttributeError, match="Cannot modify attribute"):
            system.elements = set()
    
    def test_split_system_contains(self) -> None:
        """Test that SplitSystem supports 'in' operator."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        system = SplitSystem([split1, split2])
        
        assert split1 in system
        assert split2 in system
        
        split3 = Split({1, 4}, {2, 3})
        assert split3 not in system


class TestSplitClassifications:
    """Test cases for split system classification functions."""

    def test_is_pairwise_compatible_true(self) -> None:
        """Test is_pairwise_compatible returns True for compatible splits."""
        from phylozoo.core.split.classifications import is_pairwise_compatible
        
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1}, {2, 3, 4})
        split3 = Split({1, 2, 3}, {4})
        system = SplitSystem([split1, split2, split3])
        assert is_pairwise_compatible(system) is True

    def test_is_pairwise_compatible_false(self) -> None:
        """Test is_pairwise_compatible returns False for incompatible splits."""
        from phylozoo.core.split.classifications import is_pairwise_compatible
        
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})  # Incompatible with split1
        system = SplitSystem([split1, split2])
        assert is_pairwise_compatible(system) is False

    def test_is_pairwise_compatible_empty(self) -> None:
        """Test is_pairwise_compatible returns True for empty system."""
        from phylozoo.core.split.classifications import is_pairwise_compatible
        
        system = SplitSystem()
        assert is_pairwise_compatible(system) is True

    def test_is_pairwise_compatible_single_split(self) -> None:
        """Test is_pairwise_compatible returns True for single split."""
        from phylozoo.core.split.classifications import is_pairwise_compatible
        
        split = Split({1, 2}, {3, 4})
        system = SplitSystem([split])
        assert is_pairwise_compatible(system) is True

    def test_has_all_trivial_splits_true(self) -> None:
        """Test has_all_trivial_splits returns True when all trivial splits present."""
        from phylozoo.core.split.classifications import has_all_trivial_splits
        
        split1 = Split({1}, {2, 3})
        split2 = Split({2}, {1, 3})
        split3 = Split({3}, {1, 2})
        system = SplitSystem([split1, split2, split3])
        assert has_all_trivial_splits(system) is True

    def test_has_all_trivial_splits_false(self) -> None:
        """Test has_all_trivial_splits returns False when missing trivial splits."""
        from phylozoo.core.split.classifications import has_all_trivial_splits
        
        split1 = Split({1}, {2, 3})
        split2 = Split({2}, {1, 3})
        # Missing split for element 3
        system = SplitSystem([split1, split2])
        assert has_all_trivial_splits(system) is False

    def test_has_all_trivial_splits_empty(self) -> None:
        """Test has_all_trivial_splits returns True for empty system."""
        from phylozoo.core.split.classifications import has_all_trivial_splits
        
        system = SplitSystem()
        assert has_all_trivial_splits(system) is True

    def test_has_all_trivial_splits_with_non_trivial(self) -> None:
        """Test has_all_trivial_splits works with non-trivial splits present."""
        from phylozoo.core.split.classifications import has_all_trivial_splits
        
        split1 = Split({1}, {2, 3, 4})
        split2 = Split({2}, {1, 3, 4})
        split3 = Split({3}, {1, 2, 4})
        split4 = Split({4}, {1, 2, 3})
        split5 = Split({1, 2}, {3, 4})  # Non-trivial
        system = SplitSystem([split1, split2, split3, split4, split5])
        assert has_all_trivial_splits(system) is True

    def test_is_tree_compatible_true(self) -> None:
        """Test is_tree_compatible returns True for tree-compatible system."""
        from phylozoo.core.split.classifications import is_tree_compatible
        
        split1 = Split({1}, {2, 3, 4})
        split2 = Split({2}, {1, 3, 4})
        split3 = Split({3}, {1, 2, 4})
        split4 = Split({4}, {1, 2, 3})
        split5 = Split({1, 2}, {3, 4})
        system = SplitSystem([split1, split2, split3, split4, split5])
        assert is_tree_compatible(system) is True

    def test_is_tree_compatible_false_missing_trivial(self) -> None:
        """Test is_tree_compatible returns False when trivial splits missing."""
        from phylozoo.core.split.classifications import is_tree_compatible
        
        split1 = Split({1}, {2, 3, 4})
        split2 = Split({2}, {1, 3, 4})
        # Missing trivial splits for 3 and 4
        split5 = Split({1, 2}, {3, 4})
        system = SplitSystem([split1, split2, split5])
        assert is_tree_compatible(system) is False

    def test_is_tree_compatible_false_incompatible(self) -> None:
        """Test is_tree_compatible returns False when splits are incompatible."""
        from phylozoo.core.split.classifications import is_tree_compatible
        
        split1 = Split({1}, {2, 3, 4})
        split2 = Split({2}, {1, 3, 4})
        split3 = Split({3}, {1, 2, 4})
        split4 = Split({4}, {1, 2, 3})
        split5 = Split({1, 2}, {3, 4})
        split6 = Split({1, 3}, {2, 4})  # Incompatible with split5
        system = SplitSystem([split1, split2, split3, split4, split5, split6])
        assert is_tree_compatible(system) is False

    def test_induced_splits_of_tree_are_tree_compatible(self) -> None:
        """Test that induced splits from a tree are tree-compatible."""
        from phylozoo.core.split.classifications import is_tree_compatible
        from phylozoo.core.network.sdnetwork.derivations import induced_splits
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        
        # Create a simple tree
        tree = SemiDirectedPhyNetwork(
            undirected_edges=[(0, 1), (0, 2), (0, 3), (0, 4)],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (4, {'label': 'D'})
            ]
        )
        
        splits = induced_splits(tree)
        assert is_tree_compatible(splits) is True

    def test_induced_splits_of_larger_tree_are_tree_compatible(self) -> None:
        """Test that induced splits from a larger tree are tree-compatible."""
        from phylozoo.core.split.classifications import is_tree_compatible
        from phylozoo.core.network.sdnetwork.derivations import induced_splits
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        
        # Create a larger tree with more structure (all internal nodes have degree >= 3)
        tree = SemiDirectedPhyNetwork(
            undirected_edges=[
                (5, 0), (5, 6), (5, 7),
                (0, 1), (0, 2),
                (6, 3), (6, 4)
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (4, {'label': 'D'}),
                (7, {'label': 'E'})
            ]
        )
        
        splits = induced_splits(tree)
        assert is_tree_compatible(splits) is True

    def test_induced_splits_of_directed_tree_are_tree_compatible(self) -> None:
        """Test that induced splits from a directed tree are tree-compatible."""
        from phylozoo.core.split.classifications import is_tree_compatible
        from phylozoo.core.network.dnetwork.derivations import induced_splits
        from phylozoo.core.network.dnetwork import DirectedPhyNetwork
        
        # Create a directed tree
        tree = DirectedPhyNetwork(
            edges=[
                (5, 0), (5, 6), (5, 7),
                (0, 1), (0, 2),
                (6, 3), (6, 4)
            ],
            nodes=[
                (1, {'label': 'A'}),
                (2, {'label': 'B'}),
                (3, {'label': 'C'}),
                (4, {'label': 'D'}),
                (7, {'label': 'E'})
            ]
        )
        
        splits = induced_splits(tree)
        assert is_tree_compatible(splits) is True
    
    def test_split_system_iteration(self) -> None:
        """Test that SplitSystem is iterable."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        system = SplitSystem([split1, split2])
        
        splits_list = list(system)
        assert len(splits_list) == 2
        assert split1 in splits_list
        assert split2 in splits_list


class TestSplitsystemToTree:
    """Test cases for the splitsystem_to_tree function."""
    
    def test_empty_system(self) -> None:
        """Test splitsystem_to_tree with empty system."""
        from phylozoo.core.split.algorithms import splitsystem_to_tree
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        from phylozoo.core.network.sdnetwork.classifications import is_tree
        
        system = SplitSystem([])
        tree = splitsystem_to_tree(system)
        
        assert isinstance(tree, SemiDirectedPhyNetwork)
        assert is_tree(tree)
        assert tree.number_of_nodes() == 0
        assert tree.number_of_edges() == 0
    
    def test_three_element_system(self) -> None:
        """Test splitsystem_to_tree with three elements (minimal non-trivial case)."""
        from phylozoo.core.split.algorithms import splitsystem_to_tree
        from phylozoo.core.network.sdnetwork.classifications import is_tree
        
        # Three elements: need trivial splits for all three
        split1 = Split({1}, {2, 3})
        split2 = Split({2}, {1, 3})
        split3 = Split({3}, {1, 2})
        system = SplitSystem([split1, split2, split3])
        tree = splitsystem_to_tree(system)
        
        assert is_tree(tree)
        # Should have at least 3 leaves
        assert tree.number_of_nodes() >= 3
        assert tree.number_of_edges() >= 2
    
    def test_simple_tree_compatible_system(self) -> None:
        """Test splitsystem_to_tree with simple compatible system."""
        from phylozoo.core.split.algorithms import splitsystem_to_tree
        from phylozoo.core.network.sdnetwork.classifications import is_tree
        from phylozoo.core.network.sdnetwork.derivations import induced_splits
        
        # Create a simple system with all trivial splits + one non-trivial
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1}, {2, 3, 4})
        split3 = Split({2}, {1, 3, 4})
        split4 = Split({3}, {1, 2, 4})
        split5 = Split({4}, {1, 2, 3})
        system = SplitSystem([split1, split2, split3, split4, split5])
        
        tree = splitsystem_to_tree(system)
        
        assert is_tree(tree)
        assert tree.number_of_nodes() >= 4  # At least 4 leaves
        assert tree.number_of_edges() >= 3  # At least 3 edges for 4 leaves
        
        # Check that all original splits are in the tree
        tree_splits = induced_splits(tree)
        original_str_splits = {
            Split({str(x) for x in s.set1}, {str(x) for x in s.set2})
            for s in system.splits
        }
        assert original_str_splits.issubset(tree_splits.splits)
    
    def test_incompatible_system_raises_error(self) -> None:
        """Test splitsystem_to_tree raises error for incompatible system."""
        from phylozoo.core.split.algorithms import splitsystem_to_tree
        
        # Create incompatible splits
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})  # Incompatible with split1
        split3 = Split({1}, {2, 3, 4})
        split4 = Split({2}, {1, 3, 4})
        split5 = Split({3}, {1, 2, 4})
        split6 = Split({4}, {1, 2, 3})
        system = SplitSystem([split1, split2, split3, split4, split5, split6])
        
        with pytest.raises(ValueError, match="not compatible with a tree"):
            splitsystem_to_tree(system)
    
    def test_missing_trivial_splits_raises_error(self) -> None:
        """Test splitsystem_to_tree raises error when trivial splits missing."""
        from phylozoo.core.split.algorithms import splitsystem_to_tree
        
        # Missing some trivial splits
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1}, {2, 3, 4})
        # Missing splits for 2, 3, 4
        system = SplitSystem([split1, split2])
        
        with pytest.raises(ValueError, match="not compatible with a tree"):
            splitsystem_to_tree(system)
    
    def test_check_compatibility_false(self) -> None:
        """Test splitsystem_to_tree with check_compatibility=False."""
        from phylozoo.core.split.algorithms import splitsystem_to_tree
        from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
        
        # System that would fail compatibility check
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1}, {2, 3, 4})
        system = SplitSystem([split1, split2])
        
        # Should not raise error but may produce invalid tree
        tree = splitsystem_to_tree(system, check_compatibility=False)
        assert isinstance(tree, SemiDirectedPhyNetwork)
    
    def test_small_binary_tree_roundtrip(self) -> None:
        """Test roundtrip: tree -> splits -> tree -> splits."""
        from phylozoo.core.split.algorithms import splitsystem_to_tree
        from phylozoo.core.network.sdnetwork.classifications import is_tree
        from phylozoo.core.network.sdnetwork.derivations import induced_splits
        from phylozoo.core.split import SplitSystem
        from tests.fixtures.sd_networks import SDTREE_SMALL_BINARY
        
        original_tree = SDTREE_SMALL_BINARY
        original_splits = induced_splits(original_tree)
        original_taxa = original_tree.taxa
        
        # Use all splits from the original tree (not filtered)
        split_system = SplitSystem(original_splits.splits)
        
        # Convert splits back to tree
        reconstructed_tree = splitsystem_to_tree(split_system)
        
        # Check that reconstructed tree is valid
        assert is_tree(reconstructed_tree)
        assert reconstructed_tree.number_of_nodes() > 0
        assert reconstructed_tree.number_of_edges() > 0
        
        # Check that taxa match (most important check)
        assert original_taxa == reconstructed_tree.taxa, \
            f"Taxa mismatch: original={original_taxa}, reconstructed={reconstructed_tree.taxa}"
        
        # Check that splits match
        reconstructed_splits = induced_splits(reconstructed_tree)
        assert original_splits.splits == reconstructed_splits.splits
    
    def test_medium_binary_tree_roundtrip(self) -> None:
        """Test roundtrip with medium binary tree."""
        from phylozoo.core.split.algorithms import splitsystem_to_tree
        from phylozoo.core.network.sdnetwork.classifications import is_tree
        from phylozoo.core.network.sdnetwork.derivations import induced_splits
        from phylozoo.core.split import SplitSystem
        from tests.fixtures.sd_networks import SDTREE_MEDIUM_BINARY
        
        original_tree = SDTREE_MEDIUM_BINARY
        original_splits = induced_splits(original_tree)
        original_taxa = original_tree.taxa
        
        # Use all splits from the original tree (not filtered)
        split_system = SplitSystem(original_splits.splits)
        
        # Convert splits back to tree
        reconstructed_tree = splitsystem_to_tree(split_system)
        
        # Check that reconstructed tree is valid
        assert is_tree(reconstructed_tree)
        assert reconstructed_tree.number_of_nodes() > 0
        assert reconstructed_tree.number_of_edges() > 0
        
        # Check that taxa match (most important check)
        assert original_taxa == reconstructed_tree.taxa, \
            f"Taxa mismatch: original={original_taxa}, reconstructed={reconstructed_tree.taxa}"
        
        # Check that splits match
        reconstructed_splits = induced_splits(reconstructed_tree)
        assert original_splits.splits == reconstructed_splits.splits
    
    def test_large_binary_tree_roundtrip(self) -> None:
        """Test roundtrip with large binary tree."""
        from phylozoo.core.split.algorithms import splitsystem_to_tree
        from phylozoo.core.network.sdnetwork.classifications import is_tree
        from phylozoo.core.network.sdnetwork.derivations import induced_splits
        from phylozoo.core.split import SplitSystem
        from tests.fixtures.sd_networks import SDTREE_LARGE_BINARY
        
        original_tree = SDTREE_LARGE_BINARY
        original_splits = induced_splits(original_tree)
        original_taxa = original_tree.taxa
        
        # Use all splits from the original tree (not filtered)
        split_system = SplitSystem(original_splits.splits)
        
        # Convert splits back to tree
        reconstructed_tree = splitsystem_to_tree(split_system)
        
        # Check that reconstructed tree is valid
        assert is_tree(reconstructed_tree)
        assert reconstructed_tree.number_of_nodes() > 0
        assert reconstructed_tree.number_of_edges() > 0
        
        # Check that taxa match (most important check)
        assert original_taxa == reconstructed_tree.taxa, \
            f"Taxa mismatch: original={original_taxa}, reconstructed={reconstructed_tree.taxa}"
        
        # Check that splits match
        reconstructed_splits = induced_splits(reconstructed_tree)
        assert original_splits.splits == reconstructed_splits.splits
    
    def test_non_binary_tree_roundtrip(self) -> None:
        """Test roundtrip with non-binary tree."""
        from phylozoo.core.split.algorithms import splitsystem_to_tree
        from phylozoo.core.network.sdnetwork.classifications import is_tree
        from phylozoo.core.network.sdnetwork.derivations import induced_splits
        from phylozoo.core.split import SplitSystem
        from tests.fixtures.sd_networks import SDTREE_NON_BINARY_SMALL
        
        original_tree = SDTREE_NON_BINARY_SMALL
        original_splits = induced_splits(original_tree)
        original_taxa = original_tree.taxa
        
        # Use all splits from the original tree (not filtered)
        split_system = SplitSystem(original_splits.splits)
        
        # Convert splits back to tree
        reconstructed_tree = splitsystem_to_tree(split_system)
        
        # Check that reconstructed tree is valid
        assert is_tree(reconstructed_tree)
        assert reconstructed_tree.number_of_nodes() > 0
        assert reconstructed_tree.number_of_edges() > 0
        
        # Check that taxa match (most important check)
        assert original_taxa == reconstructed_tree.taxa, \
            f"Taxa mismatch: original={original_taxa}, reconstructed={reconstructed_tree.taxa}"
        
        # Check that splits match
        reconstructed_splits = induced_splits(reconstructed_tree)
        assert original_splits.splits == reconstructed_splits.splits
    
    def test_directed_tree_roundtrip(self) -> None:
        """Test roundtrip with directed tree (converted to semi-directed)."""
        from phylozoo.core.split.algorithms import splitsystem_to_tree
        from phylozoo.core.network.sdnetwork.classifications import is_tree
        from phylozoo.core.network.sdnetwork.derivations import induced_splits
        from phylozoo.core.network.dnetwork.derivations import to_sd_network
        from phylozoo.core.split import SplitSystem
        from tests.fixtures.directed_networks import DTREE_MEDIUM_BINARY
        
        # Convert directed tree to semi-directed
        original_dtree = DTREE_MEDIUM_BINARY
        original_sdtree = to_sd_network(original_dtree)
        original_splits = induced_splits(original_sdtree)
        original_taxa = original_sdtree.taxa
        
        # Use all splits from the original tree (not filtered)
        split_system = SplitSystem(original_splits.splits)
        
        # Convert splits back to tree
        reconstructed_tree = splitsystem_to_tree(split_system)
        
        # Check that reconstructed tree is valid
        assert is_tree(reconstructed_tree)
        assert reconstructed_tree.number_of_nodes() > 0
        assert reconstructed_tree.number_of_edges() > 0
        
        # Check that taxa match (most important check)
        assert original_taxa == reconstructed_tree.taxa, \
            f"Taxa mismatch: original={original_taxa}, reconstructed={reconstructed_tree.taxa}"
        
        # Check that splits match
        reconstructed_splits = induced_splits(reconstructed_tree)
        assert original_splits.splits == reconstructed_splits.splits
    
    def test_balanced_tree_roundtrip(self) -> None:
        """Test roundtrip with balanced tree."""
        from phylozoo.core.split.algorithms import splitsystem_to_tree
        from phylozoo.core.network.sdnetwork.classifications import is_tree
        from phylozoo.core.network.sdnetwork.derivations import induced_splits
        from phylozoo.core.split import SplitSystem
        from tests.fixtures.sd_networks import SDTREE_BALANCED
        
        original_tree = SDTREE_BALANCED
        original_splits = induced_splits(original_tree)
        original_taxa = original_tree.taxa
        
        # Use all splits from the original tree (not filtered)
        split_system = SplitSystem(original_splits.splits)
        
        # Convert splits back to tree
        reconstructed_tree = splitsystem_to_tree(split_system)
        
        # Check that reconstructed tree is valid
        assert is_tree(reconstructed_tree)
        assert reconstructed_tree.number_of_nodes() > 0
        assert reconstructed_tree.number_of_edges() > 0
        
        # Check that taxa match (most important check)
        assert original_taxa == reconstructed_tree.taxa, \
            f"Taxa mismatch: original={original_taxa}, reconstructed={reconstructed_tree.taxa}"
        
        # Check that splits match
        reconstructed_splits = induced_splits(reconstructed_tree)
        assert original_splits.splits == reconstructed_splits.splits
    
    def test_star_tree(self) -> None:
        """Test splitsystem_to_tree with star tree (only trivial splits)."""
        from phylozoo.core.split.algorithms import splitsystem_to_tree
        from phylozoo.core.network.sdnetwork.classifications import is_tree
        from phylozoo.core.network.sdnetwork.derivations import induced_splits
        
        # Create star tree splits (only trivial splits)
        split1 = Split({1}, {2, 3, 4})
        split2 = Split({2}, {1, 3, 4})
        split3 = Split({3}, {1, 2, 4})
        split4 = Split({4}, {1, 2, 3})
        system = SplitSystem([split1, split2, split3, split4])
        
        tree = splitsystem_to_tree(system)
        
        assert is_tree(tree)
        # Star tree should have 5 nodes (center + 4 leaves)
        assert tree.number_of_nodes() == 5
        # Star tree should have 4 edges
        assert tree.number_of_edges() == 4
        
        # Check that all splits are present
        tree_splits = induced_splits(tree)
        original_str_splits = {
            Split({str(x) for x in s.set1}, {str(x) for x in s.set2})
            for s in system.splits
        }
        assert original_str_splits.issubset(tree_splits.splits)


