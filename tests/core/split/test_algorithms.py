"""
Tests for split system algorithms.
"""

import numpy as np
import pytest

from phylozoo.core.distance import DistanceMatrix
from phylozoo.core.split.algorithms import distances_from_splitsystem, splitsystem_to_tree
from phylozoo.core.split.base import Split
from phylozoo.core.split.weighted_splitsystem import WeightedSplitSystem


class TestDistancesFromSplitsystem:
    """Test distances_from_splitsystem function."""
    
    def test_simple_two_splits(self) -> None:
        """Test distances_from_splitsystem with two simple splits."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        weights = {split1: 2.0, split2: 1.5}
        system = WeightedSplitSystem([split1, split2], weights=weights)
        
        dm = distances_from_splitsystem(system)
        
        # Check that it's a DistanceMatrix
        assert isinstance(dm, DistanceMatrix)
        
        # Check elements
        assert set(dm.labels) == {1, 2, 3, 4}
        
        # Check distances
        # 1-2: separated by split2 only (1 in set1, 2 in set2)
        assert abs(dm.get_distance(1, 2) - 1.5) < 1e-10
        
        # 1-3: separated by split1 only (1 in set1, 3 in set2)
        assert abs(dm.get_distance(1, 3) - 2.0) < 1e-10
        
        # 1-4: separated by both splits (1 in set1, 4 in set2 for both)
        assert abs(dm.get_distance(1, 4) - 3.5) < 1e-10
        
        # 2-3: separated by both splits (2 in set1, 3 in set2 for split1; 2 in set2, 3 in set1 for split2)
        assert abs(dm.get_distance(2, 3) - 3.5) < 1e-10
        
        # 2-4: separated by split1 only (2 in set1, 4 in set2)
        assert abs(dm.get_distance(2, 4) - 2.0) < 1e-10
        
        # 3-4: separated by split2 only (3 in set1, 4 in set2)
        assert abs(dm.get_distance(3, 4) - 1.5) < 1e-10
    
    def test_symmetric_matrix(self) -> None:
        """Test that the distance matrix is symmetric."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        weights = {split1: 2.0, split2: 1.5}
        system = WeightedSplitSystem([split1, split2], weights=weights)
        
        dm = distances_from_splitsystem(system)
        
        # Check symmetry for all pairs
        for i, x in enumerate(dm.labels):
            for j, y in enumerate(dm.labels):
                assert abs(dm.get_distance(x, y) - dm.get_distance(y, x)) < 1e-10
    
    def test_zero_distance_diagonal(self) -> None:
        """Test that diagonal elements (distance to self) are zero."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        weights = {split1: 2.0, split2: 1.5}
        system = WeightedSplitSystem([split1, split2], weights=weights)
        
        dm = distances_from_splitsystem(system)
        
        # Check that distance to self is zero
        for label in dm.labels:
            assert abs(dm.get_distance(label, label)) < 1e-10
    
    def test_single_split(self) -> None:
        """Test with a single split."""
        split = Split({1, 2}, {3, 4})
        weights = {split: 5.0}
        system = WeightedSplitSystem([split], weights=weights)
        
        dm = distances_from_splitsystem(system)
        
        # Elements on same side should have distance 0
        assert abs(dm.get_distance(1, 2)) < 1e-10
        assert abs(dm.get_distance(3, 4)) < 1e-10
        
        # Elements on different sides should have distance equal to weight
        assert abs(dm.get_distance(1, 3) - 5.0) < 1e-10
        assert abs(dm.get_distance(1, 4) - 5.0) < 1e-10
        assert abs(dm.get_distance(2, 3) - 5.0) < 1e-10
        assert abs(dm.get_distance(2, 4) - 5.0) < 1e-10
    
    def test_empty_system(self) -> None:
        """Test with empty split system."""
        system = WeightedSplitSystem()
        dm = distances_from_splitsystem(system)
        
        assert isinstance(dm, DistanceMatrix)
        assert len(dm) == 0
        assert len(dm.labels) == 0
    
    def test_single_element_impossible(self) -> None:
        """Test that a split system must have at least 2 elements."""
        # A split system with one element is not possible (splits require at least 2 elements)
        # So we test that empty system works correctly
        system = WeightedSplitSystem()
        dm = distances_from_splitsystem(system)
        
        assert isinstance(dm, DistanceMatrix)
        assert len(dm) == 0
        assert len(dm.labels) == 0
    
    def test_three_elements(self) -> None:
        """Test with three elements."""
        split1 = Split({1}, {2, 3})
        split2 = Split({2}, {1, 3})
        weights = {split1: 1.0, split2: 2.0}
        system = WeightedSplitSystem([split1, split2], weights=weights)
        
        dm = distances_from_splitsystem(system)
        
        # 1-2: separated by split1 (1 in set1, 2 in set2) and split2 (1 in set2, 2 in set1)
        assert abs(dm.get_distance(1, 2) - 3.0) < 1e-10
        
        # 1-3: separated by split1 only (1 in set1, 3 in set2)
        assert abs(dm.get_distance(1, 3) - 1.0) < 1e-10
        
        # 2-3: separated by split2 only (2 in set1, 3 in set2)
        assert abs(dm.get_distance(2, 3) - 2.0) < 1e-10
    
    def test_multiple_splits_same_separation(self) -> None:
        """Test that multiple splits separating the same pair sum correctly."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 2}, {3, 4})  # Same split (but different object)
        # Actually, splits are compared by content, so we need different splits
        split2 = Split({1}, {2, 3, 4})  # This also separates 1 from 3 and 4
        split3 = Split({1, 3}, {2, 4})  # This separates 1 from 2 and 4
        
        weights = {split1: 1.0, split2: 2.0, split3: 0.5}
        system = WeightedSplitSystem([split1, split2, split3], weights=weights)
        
        dm = distances_from_splitsystem(system)
        
        # 1-3: separated by split1 (1 in set1, 3 in set2) and split2 (1 in set1, 3 in set2)
        #      but NOT by split3 (1 in set1, 3 in set1)
        assert abs(dm.get_distance(1, 3) - 3.0) < 1e-10
        
        # 1-4: separated by split1 (1 in set1, 4 in set2), split2 (1 in set1, 4 in set2),
        #      and split3 (1 in set1, 4 in set2)
        assert abs(dm.get_distance(1, 4) - 3.5) < 1e-10
    
    def test_larger_system(self) -> None:
        """Test with a larger system (5 elements)."""
        split1 = Split({1, 2}, {3, 4, 5})
        split2 = Split({1, 3}, {2, 4, 5})
        split3 = Split({1, 2, 3}, {4, 5})
        weights = {split1: 1.0, split2: 2.0, split3: 1.5}
        system = WeightedSplitSystem([split1, split2, split3], weights=weights)
        
        dm = distances_from_splitsystem(system)
        
        # Check that all elements are present
        assert set(dm.labels) == {1, 2, 3, 4, 5}
        
        # Check some distances
        # 1-4: separated by all three splits (1 in set1, 4 in set2 for all)
        assert abs(dm.get_distance(1, 4) - 4.5) < 1e-10
        
        # 1-2: separated by split2 (1 in set1, 2 in set2) and split3 (1 in set1, 2 in set1) - wait, no
        # Let's recalculate: split1 = {1,2}|{3,4,5}, split2 = {1,3}|{2,4,5}, split3 = {1,2,3}|{4,5}
        # 1-2: split1: 1 in set1, 2 in set1 -> NOT separated
        #      split2: 1 in set1, 2 in set2 -> separated (add 2.0)
        #      split3: 1 in set1, 2 in set1 -> NOT separated
        #      Total: 2.0
        assert abs(dm.get_distance(1, 2) - 2.0) < 1e-10
        
        # 4-5: split1: 4 in set2, 5 in set2 -> NOT separated
        #      split2: 4 in set2, 5 in set2 -> NOT separated
        #      split3: 4 in set2, 5 in set2 -> NOT separated
        #      Total: 0.0
        assert abs(dm.get_distance(4, 5)) < 1e-10
    
    def test_weights_accumulate(self) -> None:
        """Test that weights accumulate correctly when multiple splits separate the same pair."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        split3 = Split({1, 4}, {2, 3})
        weights = {split1: 0.5, split2: 1.0, split3: 1.5}
        system = WeightedSplitSystem([split1, split2, split3], weights=weights)
        
        dm = distances_from_splitsystem(system)
        
        # 1-2: separated by split2 (1 in set1, 2 in set2) and split3 (1 in set1, 2 in set2)
        #      but NOT by split1 (1 in set1, 2 in set1)
        assert abs(dm.get_distance(1, 2) - 2.5) < 1e-10
        
        # 1-3: separated by split1 (1 in set1, 3 in set2) and split3 (1 in set1, 3 in set2)
        #      but NOT by split2 (1 in set1, 3 in set1)
        assert abs(dm.get_distance(1, 3) - 2.0) < 1e-10
        
        # 1-4: separated by split1 (1 in set1, 4 in set2) and split2 (1 in set1, 4 in set2)
        #      but NOT by split3 (1 in set1, 4 in set1)
        assert abs(dm.get_distance(1, 4) - 1.5) < 1e-10

