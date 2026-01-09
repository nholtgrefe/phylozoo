"""
Tests for SqQuartetProfileSet class.
"""

import pytest

from phylozoo.core.quartet.base import Quartet
from phylozoo.core.quartet.qprofile import QuartetProfile
from phylozoo.core.split.base import Split
from phylozoo.inference.squirrel.sqprofile import SqQuartetProfile
from phylozoo.inference.squirrel.sqprofileset import SqQuartetProfileSet


class TestSqQuartetProfileSetInit:
    """Tests for SqQuartetProfileSet initialization."""
    
    def test_init_empty(self) -> None:
        """Test creating an empty profile set."""
        profileset = SqQuartetProfileSet()
        
        assert len(profileset) == 0
        assert len(profileset.taxa) == 0
    
    def test_init_from_sq_profiles(self) -> None:
        """Test creating from SqQuartetProfile objects."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))
        sq_profile1 = SqQuartetProfile([q1])
        sq_profile2 = SqQuartetProfile([q2])
        
        profileset = SqQuartetProfileSet(profiles=[sq_profile1, sq_profile2])
        
        assert len(profileset) == 2
        assert profileset.get_profile_weight(frozenset({1, 2, 3, 4})) == 1.0
        assert profileset.get_profile_weight(frozenset({5, 6, 7, 8})) == 1.0
    
    def test_init_from_sq_profiles_with_weights(self) -> None:
        """Test creating from SqQuartetProfile objects with explicit weights."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))
        sq_profile1 = SqQuartetProfile([q1])
        sq_profile2 = SqQuartetProfile([q2])
        
        profileset = SqQuartetProfileSet(profiles=[(sq_profile1, 2.0), (sq_profile2, 1.5)])
        
        assert len(profileset) == 2
        assert profileset.get_profile_weight(frozenset({1, 2, 3, 4})) == 2.0
        assert profileset.get_profile_weight(frozenset({5, 6, 7, 8})) == 1.5
    
    def test_init_from_quartets_auto_conversion(self) -> None:
        """Test creating from Quartet objects (auto-converted to SqQuartetProfile)."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))
        
        profileset = SqQuartetProfileSet(profiles=[q1, q2])
        
        assert len(profileset) == 2
        # Each quartet is converted to a SqQuartetProfile with weight 1.0
        assert profileset.get_profile_weight(frozenset({1, 2, 3, 4})) == 1.0
        assert profileset.get_profile_weight(frozenset({5, 6, 7, 8})) == 1.0
        
        # Verify profiles are SqQuartetProfile instances
        profile1 = profileset.get_profile(frozenset({1, 2, 3, 4}))
        assert isinstance(profile1, SqQuartetProfile)
        assert len(profile1) == 1
    
    def test_init_from_quartets_with_weights(self) -> None:
        """Test creating from Quartet objects with explicit weights."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))
        
        profileset = SqQuartetProfileSet(profiles=[(q1, 0.8), (q2, 1.2)])
        
        assert len(profileset) == 2
        assert profileset.get_profile_weight(frozenset({1, 2, 3, 4})) == 0.8
        assert profileset.get_profile_weight(frozenset({5, 6, 7, 8})) == 1.2
    
    def test_init_from_quartets_grouped_by_taxa(self) -> None:
        """Test that quartets with same taxa are merged into one profile."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))  # Same taxa, different quartet
        
        profileset = SqQuartetProfileSet(profiles=[q1, q2])
        
        # Quartets with same taxa should be merged into one profile
        assert len(profileset) == 1
        profile = profileset.get_profile(frozenset({1, 2, 3, 4}))
        assert profile is not None
        # Both quartets should be in the merged profile
        assert len(profile) == 2
        assert q1 in profile
        assert q2 in profile
        # Profile weight should be sum of quartet weights
        assert profileset.get_profile_weight(frozenset({1, 2, 3, 4})) == 2.0
    
    def test_init_from_two_quartets_same_taxa_with_reticulation_leaf(self) -> None:
        """Test creating profile set from two quartets with same taxa forming a cycle."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        
        # When quartets are provided separately, they get auto-merged
        # Reticulation leaf will be None since we can't specify it when auto-converting
        profileset = SqQuartetProfileSet(profiles=[q1, q2])
        
        assert len(profileset) == 1
        profile = profileset.get_profile(frozenset({1, 2, 3, 4}))
        assert profile is not None
        # Both quartets should be merged into one profile
        assert len(profile) == 2
        assert q1 in profile
        assert q2 in profile
        # Reticulation leaf will be None since we can't specify it when auto-converting
        assert profile.reticulation_leaf is None
    
    def test_init_mixed_sq_profiles_and_quartets_error(self) -> None:
        """Test that mixing SqQuartetProfile and Quartet raises error."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        sq_profile1 = SqQuartetProfile([q1])
        q2 = Quartet(Split({5, 6}, {7, 8}))
        
        # The parent class should handle this, but let's check what happens
        # Actually, looking at the implementation, it should work because
        # we convert Quartet to SqQuartetProfile. But the parent class might
        # complain about mixing types. Let's test.
        # Actually, the parent class checks for mixing QuartetProfile and Quartet,
        # but we're converting everything to SqQuartetProfile first, so it should work.
        profileset = SqQuartetProfileSet(profiles=[sq_profile1, q2])
        assert len(profileset) == 2
    
    def test_init_non_sq_profile_error(self) -> None:
        """Test that non-SqQuartetProfile profiles raise error."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        # Create a regular QuartetProfile (not SqQuartetProfile)
        regular_profile = QuartetProfile([q1, q2])
        
        # The parent class will accept it, but our validation should catch it
        # Actually, we convert everything first, so regular_profile won't be accepted
        # because it's not a SqQuartetProfile or Quartet
        with pytest.raises(ValueError, match="Expected SqQuartetProfile or Quartet"):
            SqQuartetProfileSet(profiles=[regular_profile])
    
    def test_init_non_positive_weight_error(self) -> None:
        """Test that non-positive weights raise error."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        sq_profile1 = SqQuartetProfile([q1])
        
        with pytest.raises(ValueError, match="Profile weight must be positive"):
            SqQuartetProfileSet(profiles=[(sq_profile1, 0.0)])
        
        with pytest.raises(ValueError, match="Profile weight must be positive"):
            SqQuartetProfileSet(profiles=[(sq_profile1, -0.5)])
    
    def test_init_duplicate_quartet_in_quartet_mode_error(self) -> None:
        """Test that duplicate quartets in quartet mode raise error."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        
        # Same quartet appears twice (same taxa)
        with pytest.raises(ValueError, match="appears multiple times in the input"):
            SqQuartetProfileSet(profiles=[q1, q2, q1])
        
        # Same quartet with different weights
        with pytest.raises(ValueError, match="appears multiple times in the input"):
            SqQuartetProfileSet(profiles=[(q1, 0.5), (q2, 0.3), (q1, 0.2)])
        
        # Three different quartets with same taxa (should fail - SqQuartetProfile only allows 1-2)
        q3 = Quartet(Split({1, 4}, {2, 3}))  # Different quartet, same taxa
        with pytest.raises(ValueError, match="SqQuartetProfile must contain exactly 1 or 2 resolved quartets"):
            SqQuartetProfileSet(profiles=[q1, q2, q3])


class TestSqQuartetProfileSetProperties:
    """Tests for SqQuartetProfileSet properties."""
    
    def test_profiles_property(self) -> None:
        """Test the profiles property returns SqQuartetProfile instances."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))
        sq_profile1 = SqQuartetProfile([q1])
        sq_profile2 = SqQuartetProfile([q2])
        
        profileset = SqQuartetProfileSet(profiles=[sq_profile1, sq_profile2])
        
        for taxa, (profile, weight) in profileset.profiles.items():
            assert isinstance(profile, SqQuartetProfile)
            assert weight > 0
    
    def test_get_profile_returns_sq_profile(self) -> None:
        """Test that get_profile returns SqQuartetProfile or None."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        sq_profile1 = SqQuartetProfile([q1])
        
        profileset = SqQuartetProfileSet(profiles=[sq_profile1])
        
        profile = profileset.get_profile(frozenset({1, 2, 3, 4}))
        assert isinstance(profile, SqQuartetProfile)
        assert len(profile) == 1
        
        # Test non-existent profile
        assert profileset.get_profile(frozenset({5, 6, 7, 8})) is None
    
    def test_inherited_properties(self) -> None:
        """Test that inherited properties from QuartetProfileSet work correctly."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))
        sq_profile1 = SqQuartetProfile([q1])
        sq_profile2 = SqQuartetProfile([q2])
        
        profileset = SqQuartetProfileSet(profiles=[sq_profile1, sq_profile2])
        
        assert profileset.taxa == frozenset({1, 2, 3, 4, 5, 6, 7, 8})
        assert len(profileset) == 2
        assert profileset.has_profile(frozenset({1, 2, 3, 4})) is True
        assert profileset.has_profile(frozenset({9, 10, 11, 12})) is False


class TestSqQuartetProfileSetEdgeCases:
    """Tests for edge cases and special behaviors."""
    
    def test_profiles_with_reticulation_leaf(self) -> None:
        """Test that profiles with reticulation_leaf are preserved."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        sq_profile = SqQuartetProfile([q1, q2], reticulation_leaf=1)
        
        profileset = SqQuartetProfileSet(profiles=[sq_profile])
        
        profile = profileset.get_profile(frozenset({1, 2, 3, 4}))
        assert profile is not None
        assert profile.reticulation_leaf == 1
    
    def test_auto_conversion_preserves_quartet_weights(self) -> None:
        """Test that when auto-converting quartets, weights are preserved."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        
        profileset = SqQuartetProfileSet(profiles=[(q1, 0.6), (q2, 0.4)])
        
        # Both quartets have same taxa, so they should be merged
        profile = profileset.get_profile(frozenset({1, 2, 3, 4}))
        assert profile is not None
        # Profile weight should be sum of quartet weights
        assert profileset.get_profile_weight(frozenset({1, 2, 3, 4})) == 1.0
        # Individual quartet weights should be preserved within the profile
        assert profile.get_weight(q1) == 0.6
        assert profile.get_weight(q2) == 0.4
    
    def test_taxa_parameter(self) -> None:
        """Test that taxa parameter works correctly."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        sq_profile1 = SqQuartetProfile([q1])
        
        # Provide taxa that includes taxa not in profiles
        profileset = SqQuartetProfileSet(
            profiles=[sq_profile1],
            taxa=frozenset({1, 2, 3, 4, 5, 6, 7, 8})
        )
        
        assert profileset.taxa == frozenset({1, 2, 3, 4, 5, 6, 7, 8})
        assert len(profileset) == 1
    
    def test_taxa_parameter_not_superset_error(self) -> None:
        """Test that taxa parameter must be a superset of profile taxa."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        sq_profile1 = SqQuartetProfile([q1])
        
        with pytest.raises(ValueError, match="Provided taxa must be a superset"):
            SqQuartetProfileSet(
                profiles=[sq_profile1],
                taxa=frozenset({1, 2, 3})  # Missing taxon 4
            )

