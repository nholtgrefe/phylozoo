"""
Tests for the weighted split system module.
"""

import pytest
from phylozoo.core.split import (
    Split,
    SplitSystem,
    WeightedSplitSystem,
    to_weightedsplitsystem,
)


class TestWeightedSplitSystem:
    """Test cases for the WeightedSplitSystem class."""

    def test_weighted_split_system_creation(self) -> None:
        """Test creating a weighted split system."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        weights = {split1: 2.5, split2: 1.0}
        system = WeightedSplitSystem([split1, split2], weights=weights)
        assert len(system) == 2
        assert split1 in system
        assert split2 in system

    def test_weighted_split_system_get_weight(self) -> None:
        """Test getting weights for splits in the system."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        weights = {split1: 2.5, split2: 1.0}
        system = WeightedSplitSystem([split1, split2], weights=weights)
        assert system.get_weight(split1) == 2.5
        assert system.get_weight(split2) == 1.0

    def test_weighted_split_system_get_weight_not_in_system(self) -> None:
        """Test getting weight for split not in system returns 0.0."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        split3 = Split({1, 4}, {2, 3})
        weights = {split1: 2.5, split2: 1.0}
        system = WeightedSplitSystem([split1, split2], weights=weights)
        assert system.get_weight(split3) == 0.0

    def test_weighted_split_system_get_weight_different_elements_raises_error(self) -> None:
        """Test that getting weight for split with different elements raises ValueError."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        weights = {split1: 2.5, split2: 1.0}
        system = WeightedSplitSystem([split1, split2], weights=weights)
        # Split with different elements
        split_different = Split({1, 2}, {5, 6})
        with pytest.raises(ValueError, match="does not cover the same elements"):
            system.get_weight(split_different)

    def test_weighted_split_system_total_weight(self) -> None:
        """Test total_weight property."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        weights = {split1: 2.5, split2: 1.0}
        system = WeightedSplitSystem([split1, split2], weights=weights)
        assert system.total_weight == 3.5

    def test_weighted_split_system_weights_property(self) -> None:
        """Test weights property returns a copy."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        weights = {split1: 2.5, split2: 1.0}
        system = WeightedSplitSystem([split1, split2], weights=weights)
        weights_dict = system.weights
        assert weights_dict == {split1: 2.5, split2: 1.0}
        # Verify it's a copy (modifying it shouldn't affect the system)
        weights_dict[split1] = 999.0
        assert system.get_weight(split1) == 2.5

    def test_weighted_split_system_zero_weight_raises_error(self) -> None:
        """Test that zero weights raise ValueError."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        split3 = Split({1, 4}, {2, 3})
        # Zero weight should raise an error
        weights = {split1: 2.5, split2: 1.0, split3: 0.0}
        with pytest.raises(ValueError, match="Weight must be positive"):
            WeightedSplitSystem([split1, split2, split3], weights=weights)

    def test_weighted_split_system_no_weights_provided(self) -> None:
        """Test creating weighted split system with no weights."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        # No weights provided - should result in empty system
        system = WeightedSplitSystem([split1, split2])
        assert len(system) == 0
        assert system.total_weight == 0.0

    def test_weighted_split_system_weights_from_keys(self) -> None:
        """Test creating weighted split system from weights dict only."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        weights = {split1: 2.5, split2: 1.0}
        # Don't provide splits parameter - should derive from weights.keys()
        system = WeightedSplitSystem(weights=weights)
        assert len(system) == 2
        assert split1 in system
        assert split2 in system

    def test_weighted_split_system_negative_weight_raises_error(self) -> None:
        """Test that negative weights raise ValueError."""
        split1 = Split({1, 2}, {3, 4})
        weights = {split1: -1.0}
        with pytest.raises(ValueError, match="Weight must be positive"):
            WeightedSplitSystem([split1], weights=weights)

    def test_weighted_split_system_weight_not_in_splits_raises_error(self) -> None:
        """Test that weight for split not in splits raises ValueError."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        weights = {split1: 2.5, split2: 1.0}
        # split2 not in splits list
        with pytest.raises(ValueError, match="Weight provided for split.*not in splits"):
            WeightedSplitSystem([split1], weights=weights)

    def test_weighted_split_system_inherits_from_splitsystem(self) -> None:
        """Test that WeightedSplitSystem inherits from SplitSystem."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        weights = {split1: 2.5, split2: 1.0}
        system = WeightedSplitSystem([split1, split2], weights=weights)
        assert isinstance(system, SplitSystem)

    def test_weighted_split_system_elements(self) -> None:
        """Test that elements property works correctly."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        weights = {split1: 2.5, split2: 1.0}
        system = WeightedSplitSystem([split1, split2], weights=weights)
        assert system.elements == {1, 2, 3, 4}

    def test_weighted_split_system_iteration(self) -> None:
        """Test that WeightedSplitSystem is iterable."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        weights = {split1: 2.5, split2: 1.0}
        system = WeightedSplitSystem([split1, split2], weights=weights)
        splits_list = list(system)
        assert len(splits_list) == 2
        assert split1 in splits_list
        assert split2 in splits_list

    def test_weighted_split_system_repr(self) -> None:
        """Test string representation of weighted split system."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        weights = {split1: 2.5, split2: 1.0}
        system = WeightedSplitSystem([split1, split2], weights=weights)
        repr_str = repr(system)
        assert "WeightedSplitSystem" in repr_str
        assert "2.5" in repr_str or "1.0" in repr_str

    def test_weighted_split_system_immutability(self) -> None:
        """Test that WeightedSplitSystem is immutable."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        weights = {split1: 2.5, split2: 1.0}
        system = WeightedSplitSystem([split1, split2], weights=weights)
        
        # Try to modify attributes
        with pytest.raises(AttributeError, match="Cannot modify attribute"):
            system.splits = set()
        
        with pytest.raises(AttributeError, match="Cannot modify attribute"):
            system.elements = set()
        
        with pytest.raises(AttributeError, match="Cannot modify attribute"):
            system.weights = {}

    def test_weighted_split_system_empty_system(self) -> None:
        """Test creating an empty weighted split system."""
        system = WeightedSplitSystem()
        assert len(system) == 0
        assert system.total_weight == 0.0
        assert system.weights == {}

    def test_weighted_split_system_total_weight_cached(self) -> None:
        """Test that total_weight is cached (computed once)."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        weights = {split1: 2.5, split2: 1.0}
        system = WeightedSplitSystem([split1, split2], weights=weights)
        # Access multiple times - should return same value
        assert system.total_weight == 3.5
        assert system.total_weight == 3.5
        assert system.total_weight == 3.5


class TestToWeightedSplitSystem:
    """Test cases for the to_weightedsplitsystem function."""

    def test_to_weightedsplitsystem_basic(self) -> None:
        """Test converting SplitSystem to WeightedSplitSystem."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        system = SplitSystem([split1, split2])
        weighted = to_weightedsplitsystem(system)
        assert isinstance(weighted, WeightedSplitSystem)
        assert len(weighted) == 2
        assert weighted.get_weight(split1) == 1.0
        assert weighted.get_weight(split2) == 1.0

    def test_to_weightedsplitsystem_custom_weight(self) -> None:
        """Test converting with custom default weight."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        system = SplitSystem([split1, split2])
        weighted = to_weightedsplitsystem(system, default_weight=2.5)
        assert weighted.get_weight(split1) == 2.5
        assert weighted.get_weight(split2) == 2.5
        assert weighted.total_weight == 5.0

    def test_to_weightedsplitsystem_negative_weight_raises_error(self) -> None:
        """Test that negative default_weight raises ValueError."""
        split1 = Split({1, 2}, {3, 4})
        system = SplitSystem([split1])
        with pytest.raises(ValueError, match="default_weight must be positive"):
            to_weightedsplitsystem(system, default_weight=-1.0)

    def test_to_weightedsplitsystem_zero_weight_raises_error(self) -> None:
        """Test that zero default_weight raises ValueError."""
        split1 = Split({1, 2}, {3, 4})
        system = SplitSystem([split1])
        with pytest.raises(ValueError, match="default_weight must be positive"):
            to_weightedsplitsystem(system, default_weight=0.0)

    def test_to_weightedsplitsystem_empty_system(self) -> None:
        """Test converting empty SplitSystem."""
        system = SplitSystem()
        weighted = to_weightedsplitsystem(system)
        assert isinstance(weighted, WeightedSplitSystem)
        assert len(weighted) == 0
        assert weighted.total_weight == 0.0

    def test_to_weightedsplitsystem_preserves_elements(self) -> None:
        """Test that conversion preserves elements."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        system = SplitSystem([split1, split2])
        weighted = to_weightedsplitsystem(system)
        assert weighted.elements == {1, 2, 3, 4}

