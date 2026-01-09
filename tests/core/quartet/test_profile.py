"""
Tests for QuartetProfile class.
"""

from typing import Mapping

import pytest

from phylozoo.core.quartet import Quartet, QuartetProfile
from phylozoo.core.split.base import Split


class TestQuartetProfileInit:
    """Tests for QuartetProfile initialization."""
    
    def test_init_from_dict(self) -> None:
        """Test creating a profile from a dictionary."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        profile = QuartetProfile({q1: 0.8, q2: 0.2})
        
        assert profile.taxa == frozenset({1, 2, 3, 4})
        assert len(profile) == 2
        assert profile.get_weight(q1) == 0.8
        assert profile.get_weight(q2) == 0.2
        assert profile.total_weight == 1.0
    
    def test_init_from_list_quartets(self) -> None:
        """Test creating a profile from a list of quartets."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        profile = QuartetProfile([q1, q2])
        
        assert profile.taxa == frozenset({1, 2, 3, 4})
        assert len(profile) == 2
        assert profile.get_weight(q1) == 1.0
        assert profile.get_weight(q2) == 1.0
        assert profile.total_weight == 2.0
    
    def test_init_from_list_tuples(self) -> None:
        """Test creating a profile from a list of (quartet, weight) tuples."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        profile = QuartetProfile([(q1, 0.7), (q2, 0.3)])
        
        assert profile.taxa == frozenset({1, 2, 3, 4})
        assert len(profile) == 2
        assert profile.get_weight(q1) == 0.7
        assert profile.get_weight(q2) == 0.3
        assert profile.total_weight == 1.0
    
    def test_init_empty_error(self) -> None:
        """Test that empty profile raises ValueError."""
        with pytest.raises(ValueError, match="QuartetProfile must have at least one quartet"):
            QuartetProfile([])
        
        with pytest.raises(ValueError, match="QuartetProfile must have at least one quartet"):
            QuartetProfile({})
    
    def test_init_different_taxa_error(self) -> None:
        """Test that quartets with different taxa raise ValueError."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))
        
        with pytest.raises(ValueError, match="All quartets must have the same taxa"):
            QuartetProfile({q1: 0.5, q2: 0.5})
        
        with pytest.raises(ValueError, match="All quartets must have the same taxa"):
            QuartetProfile([q1, q2])
    
    def test_init_non_positive_weight_error(self) -> None:
        """Test that non-positive weights raise ValueError."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        
        with pytest.raises(ValueError, match="Weight must be positive"):
            QuartetProfile({q1: 0.0, q2: 1.0})
        
        with pytest.raises(ValueError, match="Weight must be positive"):
            QuartetProfile({q1: -0.5, q2: 1.0})
        
        with pytest.raises(ValueError, match="Weight must be positive"):
            QuartetProfile([(q1, 0.0), (q2, 1.0)])
    
    def test_init_duplicate_quartet_error(self) -> None:
        """Test that duplicate quartets raise ValueError."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        
        # Duplicate in list of quartets
        with pytest.raises(ValueError, match="appears multiple times in the input"):
            QuartetProfile([q1, q2, q1])
        
        # Duplicate in list of tuples
        with pytest.raises(ValueError, match="appears multiple times in the input"):
            QuartetProfile([(q1, 0.5), (q2, 0.3), (q1, 0.2)])
        
        # Duplicate in dict (can't actually happen with dict syntax, but test explicit dict)
        # Note: dict syntax automatically overwrites, so we can't test this directly
        # But we can test that the validation would catch it if we construct it manually
        quartets_dict = {q1: 0.5, q2: 0.3}
        quartets_dict[q1] = 0.2  # This overwrites, so no error
        # The dict itself doesn't have duplicates, so this is fine
        profile = QuartetProfile(quartets_dict)
        assert profile.get_weight(q1) == 0.2  # Last value wins
    
    def test_init_single_quartet(self) -> None:
        """Test creating a profile with a single quartet."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profile = QuartetProfile([q1])
        
        assert profile.taxa == frozenset({1, 2, 3, 4})
        assert len(profile) == 1
        assert profile.get_weight(q1) == 1.0
        assert profile.total_weight == 1.0
    
    def test_init_mixed_star_and_resolved(self) -> None:
        """Test creating a profile with both star and resolved quartets on same taxa."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet({1, 2, 3, 4})  # Star tree
        profile = QuartetProfile({q1: 0.6, q2: 0.4})
        
        assert profile.taxa == frozenset({1, 2, 3, 4})
        assert len(profile) == 2
        assert profile.get_weight(q1) == 0.6
        assert profile.get_weight(q2) == 0.4


class TestQuartetProfileProperties:
    """Tests for QuartetProfile properties."""
    
    def test_taxa_property(self) -> None:
        """Test taxa property."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profile = QuartetProfile([q1])
        
        assert profile.taxa == frozenset({1, 2, 3, 4})
        assert isinstance(profile.taxa, frozenset)
    
    def test_quartets_property(self) -> None:
        """Test quartets property."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        profile = QuartetProfile({q1: 0.8, q2: 0.2})
        
        quartets = profile.quartets
        assert isinstance(quartets, Mapping)
        assert q1 in quartets
        assert q2 in quartets
        assert quartets[q1] == 0.8
        assert quartets[q2] == 0.2
    
    def test_total_weight_property(self) -> None:
        """Test total_weight property."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        profile = QuartetProfile({q1: 0.3, q2: 0.7})
        
        assert profile.total_weight == 1.0
    
    def test_immutability(self) -> None:
        """Test that profile is immutable."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profile = QuartetProfile([q1])
        
        with pytest.raises(AttributeError, match="Cannot modify attribute"):
            profile._taxa = frozenset({5, 6, 7, 8})
        
        with pytest.raises(AttributeError, match="Cannot modify attribute"):
            profile._quartets = {}
    
    def test_split_single_resolved(self) -> None:
        """Test split property with single resolved quartet."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profile = QuartetProfile([q1])
        
        split = profile.split
        assert split is not None
        assert split.set1 == {1, 2} or split.set1 == {3, 4}
        assert split.set2 == {3, 4} or split.set2 == {1, 2}
    
    def test_split_single_star(self) -> None:
        """Test split property with single star tree."""
        star = Quartet({1, 2, 3, 4})
        profile = QuartetProfile([star])
        
        assert profile.split is None
    
    def test_split_multiple_quartets(self) -> None:
        """Test split property with multiple quartets."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        profile = QuartetProfile([q1, q2])
        
        assert profile.split is None
    
    def test_split_multiple_with_star(self) -> None:
        """Test split property with multiple quartets including star."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        star = Quartet({1, 2, 3, 4})
        profile = QuartetProfile([q1, star])
        
        assert profile.split is None
    
    def test_split_cached(self) -> None:
        """Test that split property is cached."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profile = QuartetProfile([q1])
        
        split1 = profile.split
        split2 = profile.split
        
        # Should return the same object (cached)
        assert split1 is split2
    
    def test_circular_orderings_single_resolved(self) -> None:
        """Test circular_orderings property with single resolved quartet."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profile = QuartetProfile([q1])
        
        orderings = profile.circular_orderings
        assert orderings is not None
        assert len(orderings) == 2
        # Should match the quartet's circular orderings
        assert orderings == q1.circular_orderings
    
    def test_circular_orderings_single_star(self) -> None:
        """Test circular_orderings property with single star tree."""
        star = Quartet({1, 2, 3, 4})
        profile = QuartetProfile([star])
        
        orderings = profile.circular_orderings
        assert orderings is not None
        assert len(orderings) == 3
        # Should match the star's circular orderings
        assert orderings == star.circular_orderings
    
    def test_circular_orderings_two_resolved(self) -> None:
        """Test circular_orderings property with two resolved quartets."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        profile = QuartetProfile([q1, q2])
        
        orderings = profile.circular_orderings
        assert orderings is not None
        # Intersection of two quartets' orderings
        expected = q1.circular_orderings & q2.circular_orderings
        assert orderings == expected
        # Should have at least 1 ordering (the one that satisfies both)
        assert len(orderings) >= 1
    
    def test_circular_orderings_three_resolved(self) -> None:
        """Test circular_orderings property with three resolved quartets."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        q3 = Quartet(Split({1, 4}, {2, 3}))
        profile = QuartetProfile([q1, q2, q3])
        
        orderings = profile.circular_orderings
        # All three resolved quartets: no ordering satisfies all
        assert orderings is None
    
    def test_circular_orderings_resolved_and_star(self) -> None:
        """Test circular_orderings property with resolved quartet and star."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        star = Quartet({1, 2, 3, 4})
        profile = QuartetProfile([q1, star])
        
        orderings = profile.circular_orderings
        assert orderings is not None
        # Star doesn't constrain, so should return resolved quartet's orderings
        assert orderings == q1.circular_orderings
        assert len(orderings) == 2
    
    def test_circular_orderings_two_stars(self) -> None:
        """Test circular_orderings property with star tree."""
        # Note: Can't have duplicate stars (same quartet), so test with one star
        star1 = Quartet({1, 2, 3, 4})
        profile = QuartetProfile([star1])
        
        orderings = profile.circular_orderings
        assert orderings is not None
        # Stars don't constrain, should return all 3 orderings
        assert len(orderings) == 3
    
    def test_circular_orderings_cached(self) -> None:
        """Test that circular_orderings property is cached."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profile = QuartetProfile([q1])
        
        orderings1 = profile.circular_orderings
        orderings2 = profile.circular_orderings
        
        # Should return the same object (cached)
        assert orderings1 is orderings2
    
    def test_circular_orderings_intersection_empty(self) -> None:
        """Test circular_orderings when intersection of two resolved quartets is empty."""
        # Create two resolved quartets with no common ordering
        # This is tricky - let's use quartets that have different splits
        # Actually, any two different resolved quartets should have at least one common ordering
        # But three resolved quartets should have none
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        q3 = Quartet(Split({1, 4}, {2, 3}))
        
        # Two quartets should have intersection
        profile12 = QuartetProfile([q1, q2])
        orderings12 = profile12.circular_orderings
        assert orderings12 is not None
        assert len(orderings12) >= 1
        
        # Three quartets should have no intersection
        profile123 = QuartetProfile([q1, q2, q3])
        orderings123 = profile123.circular_orderings
        assert orderings123 is None


class TestQuartetProfileMethods:
    """Tests for QuartetProfile methods."""
    
    def test_get_weight(self) -> None:
        """Test get_weight method."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        profile = QuartetProfile({q1: 0.8, q2: 0.2})
        
        assert profile.get_weight(q1) == 0.8
        assert profile.get_weight(q2) == 0.2
    
    def test_get_weight_not_found(self) -> None:
        """Test get_weight for quartet not in profile."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        q3 = Quartet(Split({1, 4}, {2, 3}))
        profile = QuartetProfile({q1: 0.8, q2: 0.2})
        
        assert profile.get_weight(q3) == 0.0
    
    def test_len(self) -> None:
        """Test __len__ method."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        q3 = Quartet(Split({1, 4}, {2, 3}))
        profile = QuartetProfile([q1, q2, q3])
        
        assert len(profile) == 3
    
    def test_iter(self) -> None:
        """Test __iter__ method."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        profile = QuartetProfile([q1, q2])
        
        quartets_list = list(profile)
        assert len(quartets_list) == 2
        assert q1 in quartets_list
        assert q2 in quartets_list
    
    def test_contains(self) -> None:
        """Test __contains__ method."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        q3 = Quartet(Split({1, 4}, {2, 3}))
        profile = QuartetProfile([q1, q2])
        
        assert q1 in profile
        assert q2 in profile
        assert q3 not in profile
    
    def test_repr(self) -> None:
        """Test __repr__ method."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profile = QuartetProfile([q1])
        
        repr_str = repr(profile)
        assert "QuartetProfile" in repr_str
        assert "taxa" in repr_str or "1" in repr_str


class TestQuartetProfileEdgeCases:
    """Tests for QuartetProfile edge cases."""
    
    def test_duplicate_quartets_in_list(self) -> None:
        """Test that duplicate quartets in list raise an error."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        # Same quartet appears twice - should raise error
        with pytest.raises(ValueError, match="appears multiple times in the input"):
            QuartetProfile([q1, q1, q2])
    
    def test_duplicate_quartets_in_dict(self) -> None:
        """Test that duplicate quartets in dict use last weight."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        # Same quartet with different weights (dict behavior)
        profile = QuartetProfile({q1: 0.5, q1: 0.8})
        
        # Dict will only have one entry (last one wins)
        assert len(profile) == 1
        assert profile.get_weight(q1) == 0.8
    
    def test_large_weights(self) -> None:
        """Test profile with large weights."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        profile = QuartetProfile({q1: 100.0, q2: 200.0})
        
        assert profile.total_weight == 300.0
        assert profile.get_weight(q1) == 100.0
        assert profile.get_weight(q2) == 200.0
    
    def test_small_weights(self) -> None:
        """Test profile with very small weights."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        profile = QuartetProfile({q1: 0.0001, q2: 0.0002})
        
        assert profile.total_weight == pytest.approx(0.0003)
        assert profile.get_weight(q1) == 0.0001
        assert profile.get_weight(q2) == 0.0002
    
    def test_many_quartets(self) -> None:
        """Test profile with many quartets."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        q3 = Quartet(Split({1, 4}, {2, 3}))
        star = Quartet({1, 2, 3, 4})
        
        profile = QuartetProfile([q1, q2, q3, star])
        
        assert len(profile) == 4
        assert profile.total_weight == 4.0
        assert all(q in profile for q in [q1, q2, q3, star])


class TestQuartetProfileIsResolved:
    """Tests for QuartetProfile.is_resolved method."""
    
    def test_all_resolved_quartets(self) -> None:
        """Test is_resolved with all resolved quartets."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        profile = QuartetProfile([q1, q2])
        
        assert profile.is_resolved() is True
    
    def test_single_resolved_quartet(self) -> None:
        """Test is_resolved with single resolved quartet."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profile = QuartetProfile([q1])
        
        assert profile.is_resolved() is True
    
    def test_mixed_resolved_and_star(self) -> None:
        """Test is_resolved with mix of resolved and star quartets."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        star = Quartet({1, 2, 3, 4})
        profile = QuartetProfile([q1, star])
        
        assert profile.is_resolved() is False
    
    def test_all_star_quartets(self) -> None:
        """Test is_resolved with all star quartets."""
        star1 = Quartet({1, 2, 3, 4})
        # Can't have duplicate stars, so just use one
        profile = QuartetProfile([star1])
        
        assert profile.is_resolved() is False
    
    def test_single_star_quartet(self) -> None:
        """Test is_resolved with single star quartet."""
        star = Quartet({1, 2, 3, 4})
        profile = QuartetProfile([star])
        
        assert profile.is_resolved() is False
    
    def test_three_resolved_quartets(self) -> None:
        """Test is_resolved with three resolved quartets."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        q3 = Quartet(Split({1, 4}, {2, 3}))
        profile = QuartetProfile([q1, q2, q3])
        
        assert profile.is_resolved() is True
    
    def test_multiple_resolved_one_star(self) -> None:
        """Test is_resolved with multiple resolved and one star."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        q3 = Quartet(Split({1, 4}, {2, 3}))
        star = Quartet({1, 2, 3, 4})
        profile = QuartetProfile([q1, q2, q3, star])
        
        assert profile.is_resolved() is False
    
    def test_is_resolved_with_weights(self) -> None:
        """Test is_resolved works correctly with weighted quartets."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        profile = QuartetProfile({q1: 0.7, q2: 0.3})
        
        assert profile.is_resolved() is True
        
        star = Quartet({1, 2, 3, 4})
        profile2 = QuartetProfile({q1: 0.7, star: 0.3})
        
        assert profile2.is_resolved() is False

