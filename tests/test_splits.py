"""
Tests for the splits module.
"""

import pytest
from phylozoo.core import Split, SplitSet, SplitSystem, QuartetSplitSet


class TestSplit:
    """Test cases for the Split class."""

    def test_split_creation(self) -> None:
        """Test creating a valid split."""
        split = Split({1, 2}, {3, 4})
        assert split.set1 == {1, 2}
        assert split.set2 == {3, 4}
        assert split.elements == {1, 2, 3, 4}

    def test_is_trivial(self) -> None:
        """Test the is_trivial method."""
        trivial_split = Split({1}, {2, 3, 4})
        assert trivial_split.is_trivial() is True

        non_trivial_split = Split({1, 2}, {3, 4})
        assert non_trivial_split.is_trivial() is False

    def test_split_repr(self) -> None:
        """Test string representation of a split."""
        split = Split({1, 2}, {3, 4})
        repr_str = repr(split)
        assert "Split" in repr_str


class TestSplitSystem:
    """Test cases for the SplitSystem class."""

    def test_split_system_creation(self) -> None:
        """Test creating a split system."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        system = SplitSystem([split1, split2])
        assert len(system) == 2

    def test_add_split(self) -> None:
        """Test adding a split to the system."""
        system = SplitSystem()
        split = Split({1, 2}, {3, 4})
        system.add(split)
        assert len(system) == 1

    def test_split_system_repr(self) -> None:
        """Test string representation of a split system."""
        system = SplitSystem()
        repr_str = repr(system)
        assert "SplitSystem" in repr_str


class TestSplitSet:
    """Test cases for the SplitSet class."""

    def test_split_set_creation(self) -> None:
        """Test creating a split set."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        split_set = SplitSet([split1, split2])
        assert len(split_set) == 2


class TestQuartetSplitSet:
    """Test cases for the QuartetSplitSet class."""

    def test_quartet_split_set_creation(self, sample_quartet_split_set: QuartetSplitSet) -> None:
        """Test creating a quartet split set."""
        assert len(sample_quartet_split_set) == 2
        assert sample_quartet_split_set.elements == {1, 2, 3, 4}
