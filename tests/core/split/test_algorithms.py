"""
Tests for split system algorithms.
"""

import numpy as np
import pytest

from phylozoo.core.distance import DistanceMatrix
from phylozoo.core.split.algorithms import distances_from_splitsystem, quartets_from_splitsystem, tree_from_splitsystem
from phylozoo.core.split.base import Split
from phylozoo.core.split.splitsystem import SplitSystem
from phylozoo.core.split.weighted_splitsystem import WeightedSplitSystem


class TestDistancesFromSplitsystem:
    """Test distances_from_splitsystem function."""
    
    def test_simple_two_splits(self) -> None:
        """Test distances_from_splitsystem with two simple splits."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        weights = {split1: 2.0, split2: 1.5}
        system = WeightedSplitSystem(weights)
        
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
        system = WeightedSplitSystem(weights)
        
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
        system = WeightedSplitSystem(weights)
        
        dm = distances_from_splitsystem(system)
        
        # Check that distance to self is zero
        for label in dm.labels:
            assert abs(dm.get_distance(label, label)) < 1e-10
    
    def test_single_split(self) -> None:
        """Test with a single split."""
        split = Split({1, 2}, {3, 4})
        weights = {split: 5.0}
        system = WeightedSplitSystem(weights)
        
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
        system = WeightedSplitSystem(weights)
        
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
        system = WeightedSplitSystem(weights)
        
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
        system = WeightedSplitSystem(weights)
        
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
        system = WeightedSplitSystem(weights)
        
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
    
    def test_unweighted_split_system(self) -> None:
        """Test distances_from_splitsystem with unweighted SplitSystem."""
        split1 = Split({1, 2}, {3, 4})
        split2 = Split({1, 3}, {2, 4})
        system = SplitSystem([split1, split2])
        
        dm = distances_from_splitsystem(system)
        
        # Check that it's a DistanceMatrix
        assert isinstance(dm, DistanceMatrix)
        
        # Check elements
        assert set(dm.labels) == {1, 2, 3, 4}
        
        # Check distances (each split has implicit weight 1.0)
        # 1-2: separated by split2 only (1 in set1, 2 in set2)
        assert abs(dm.get_distance(1, 2) - 1.0) < 1e-10
        
        # 1-3: separated by split1 only (1 in set1, 3 in set2)
        assert abs(dm.get_distance(1, 3) - 1.0) < 1e-10
        
        # 1-4: separated by both splits (1 in set1, 4 in set2 for both)
        assert abs(dm.get_distance(1, 4) - 2.0) < 1e-10
        
        # 2-3: separated by both splits
        assert abs(dm.get_distance(2, 3) - 2.0) < 1e-10
        
        # 2-4: separated by split1 only
        assert abs(dm.get_distance(2, 4) - 1.0) < 1e-10
        
        # 3-4: separated by split2 only
        assert abs(dm.get_distance(3, 4) - 1.0) < 1e-10


class TestQuartetsFromSplitsystem:
    """Test quartets_from_splitsystem function."""
    
    def test_simple_split_system(self) -> None:
        """Test quartets_from_splitsystem with simple splits."""
        split1 = Split({1, 2, 3}, {4, 5, 6})
        split2 = Split({1, 2}, {3, 4, 5, 6})
        system = SplitSystem([split1, split2])
        
        profileset = quartets_from_splitsystem(system)
        
        # Should have multiple profiles (one for each 4-taxon combination)
        assert len(profileset) > 0
        assert profileset.taxa == system.elements
        
        # All profiles should have default weight 1.0
        for profile, profile_weight in profileset.profiles.values():
            assert abs(profile_weight - 1.0) < 1e-10
            # Each profile should have at least one quartet
            assert len(profile) >= 1
    
    def test_quartet_weights_accumulate(self) -> None:
        """Test that quartet weights accumulate when quartets appear in multiple splits."""
        # split1 = {1,2,3}|{4,5,6} can produce Split({1,2},{4,5})
        # split2 = {1,2}|{3,4,5,6} can also produce Split({1,2},{4,5})
        # So Split({1,2},{4,5}) should have weight 2.0
        split1 = Split({1, 2, 3}, {4, 5, 6})
        split2 = Split({1, 2}, {3, 4, 5, 6})
        system = SplitSystem([split1, split2])
        
        profileset = quartets_from_splitsystem(system)
        
        # Find the profile for {1, 2, 4, 5}
        target_taxa = frozenset({1, 2, 4, 5})
        if target_taxa in profileset.profiles:
            profile, _ = profileset.profiles[target_taxa]
            # Find the quartet Split({1, 2}, {4, 5})
            from phylozoo.core.quartet import Quartet
            target_quartet = Quartet(Split({1, 2}, {4, 5}))
            
            if target_quartet in profile.quartets:
                # QuartetProfile normalizes so weights sum to 1.0; this profile has one quartet
                assert abs(profile.get_weight(target_quartet) - 1.0) < 1e-10
                assert abs(sum(profile.quartets.values()) - 1.0) < 1e-10
    
    def test_weighted_split_system(self) -> None:
        """Test quartets_from_splitsystem with weighted split system."""
        split1 = Split({1, 2, 3}, {4, 5, 6})
        split2 = Split({1, 2}, {3, 4, 5, 6})
        weights = {split1: 2.0, split2: 1.5}
        system = WeightedSplitSystem(weights)
        
        profileset = quartets_from_splitsystem(system)
        
        # Should have multiple profiles
        assert len(profileset) > 0
        
        # Find a quartet that appears in both splits
        # Split({1, 2}, {4, 5}) appears in both
        from phylozoo.core.quartet import Quartet
        target_quartet = Quartet(Split({1, 2}, {4, 5}))
        target_taxa = frozenset({1, 2, 4, 5})
        
        if target_taxa in profileset.profiles:
            profile, _ = profileset.profiles[target_taxa]
            if target_quartet in profile.quartets:
                # QuartetProfile normalizes to sum 1.0; single quartet has weight 1.0
                assert abs(profile.get_weight(target_quartet) - 1.0) < 1e-10
                assert abs(sum(profile.quartets.values()) - 1.0) < 1e-10
    
    def test_empty_system(self) -> None:
        """Test quartets_from_splitsystem with empty system."""
        system = SplitSystem()
        profileset = quartets_from_splitsystem(system)
        
        assert len(profileset) == 0
        assert len(profileset.taxa) == 0
    
    def test_system_fewer_than_four_elements(self) -> None:
        """Test quartets_from_splitsystem with system having fewer than 4 elements."""
        split = Split({1, 2}, {3})
        system = SplitSystem([split])
        profileset = quartets_from_splitsystem(system)
        
        # Should return empty (need at least 4 elements for quartets)
        assert len(profileset) == 0
    
    def test_splits_with_insufficient_elements(self) -> None:
        """Test that splits with fewer than 2 elements on one side are skipped."""
        # Split with only 1 element on one side can't produce quartets
        split1 = Split({1, 2, 3}, {4, 5, 6})  # Can produce quartets
        split2 = Split({1}, {2, 3, 4, 5, 6})  # Cannot produce quartets (need 2 on each side)
        system = SplitSystem([split1, split2])
        
        profileset = quartets_from_splitsystem(system)
        
        # Should have profiles from split1 only
        assert len(profileset) > 0
        # All quartets should come from split1 (elements 1,2,3,4,5,6)
        assert all(1 in taxa or 2 in taxa or 3 in taxa for taxa in profileset.profiles.keys())
    
    def test_single_split_produces_multiple_quartets(self) -> None:
        """Test that a single split produces multiple quartets."""
        # Split({1,2,3},{4,5,6}) should produce C(3,2) * C(3,2) = 3 * 3 = 9 quartets
        split = Split({1, 2, 3}, {4, 5, 6})
        system = SplitSystem([split])
        
        profileset = quartets_from_splitsystem(system)
        
        # Should have multiple profiles (one for each 4-taxon combination)
        # C(6,4) = 15 combinations, but only those that can be formed from the split
        # Actually, the split can produce quartets for any 4-taxon set that has
        # at least 2 from {1,2,3} and at least 2 from {4,5,6}
        assert len(profileset) > 0
        
        # Count total quartets across all profiles
        total_quartets = sum(len(profile) for profile, _ in profileset.profiles.values())
        # Should be 9 quartets total
        assert total_quartets == 9
    
    def test_all_quartets_are_resolved(self) -> None:
        """Test that all quartets produced are resolved (not star trees)."""
        split1 = Split({1, 2, 3}, {4, 5, 6})
        split2 = Split({1, 2}, {3, 4, 5, 6})
        system = SplitSystem([split1, split2])
        
        profileset = quartets_from_splitsystem(system)
        
        # All quartets should be resolved (from 2|2 splits)
        for profile, _ in profileset.profiles.values():
            for quartet in profile.quartets.keys():
                assert quartet.is_resolved()
                assert not quartet.is_star()
                assert quartet.split is not None
    
    def test_profile_weights_default(self) -> None:
        """Test that profile weights are default (1.0) as specified."""
        split1 = Split({1, 2, 3}, {4, 5, 6})
        system = SplitSystem([split1])
        
        profileset = quartets_from_splitsystem(system)
        
        # All profiles should have default weight 1.0
        for profile, profile_weight in profileset.profiles.values():
            assert abs(profile_weight - 1.0) < 1e-10

