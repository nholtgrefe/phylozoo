"""
Tests for the splits module.
"""

import pytest
from phylozoo.core.structure import Split, SplitSystem


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
        assert split1.is_compatible(split2) is True

    def test_is_compatible_true_other_way(self) -> None:
        """Test is_compatible returns True when other split has subset."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 2, 3}, {4})
        assert split1.is_compatible(split2) is True

    def test_is_compatible_false_different_elements(self) -> None:
        """Test is_compatible returns False for splits with different elements."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 2}, {3, 5})
        assert split1.is_compatible(split2) is False

    def test_is_compatible_false_same_elements_incompatible(self) -> None:
        """Test is_compatible returns False for incompatible splits with same elements."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        assert split1.is_compatible(split2) is False

    def test_is_compatible_raises_error(self) -> None:
        """Test that is_compatible raises error for non-Split."""
        split = Split({1, 2}, {3, 4})
        with pytest.raises(ValueError, match="Not a Split instance"):
            split.is_compatible("not a split")

    def test_is_subsplit_true(self) -> None:
        """Test is_subsplit returns True when one split is subsplit of another."""
        split1 = Split({1, 2, 6}, {3, 4, 5})
        split2 = Split({1, 2}, {3, 4})
        assert split2.is_subsplit(split1) is True

    def test_is_subsplit_false(self) -> None:
        """Test is_subsplit returns False when not a subsplit."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        assert split1.is_subsplit(split2) is False

    def test_is_subsplit_raises_error(self) -> None:
        """Test that is_subsplit raises error for non-Split."""
        split = Split({1, 2}, {3, 4})
        with pytest.raises(ValueError, match="Not a Split instance"):
            split.is_subsplit("not a split")

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
    
    def test_split_system_iteration(self) -> None:
        """Test that SplitSystem is iterable."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        system = SplitSystem([split1, split2])
        
        splits_list = list(system)
        assert len(splits_list) == 2
        assert split1 in splits_list
        assert split2 in splits_list


