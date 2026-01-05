"""
Tests for QuartetProfileSet class.
"""

from typing import Mapping

import pytest

from phylozoo.core.quartet import Quartet, QuartetProfile, QuartetProfileSet
from phylozoo.core.split.base import Split


class TestQuartetProfileSetInit:
    """Tests for QuartetProfileSet initialization."""
    
    def test_init_empty(self) -> None:
        """Test creating an empty profile set."""
        profileset = QuartetProfileSet()
        
        assert len(profileset) == 0
        assert len(profileset.taxa) == 0
        assert profileset.is_dense is True
    
    def test_init_from_quartet_profiles(self) -> None:
        """Test creating from QuartetProfile objects."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        profile1 = QuartetProfile({q1: 0.8, q2: 0.2})
        profile2 = QuartetProfile([Quartet(Split({5, 6}, {7, 8}))])
        
        profileset = QuartetProfileSet(profiles=[profile1, profile2])
        
        assert len(profileset) == 2
        assert profileset.get_profile_weight(frozenset({1, 2, 3, 4})) == 1.0
        assert profileset.get_profile_weight(frozenset({5, 6, 7, 8})) == 1.0
    
    def test_init_from_quartet_profiles_with_weights(self) -> None:
        """Test creating from QuartetProfile objects with explicit weights."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        profile1 = QuartetProfile({q1: 0.8, q2: 0.2})
        profile2 = QuartetProfile([Quartet(Split({5, 6}, {7, 8}))])
        
        profileset = QuartetProfileSet(profiles=[(profile1, 2.0), (profile2, 1.5)])
        
        assert len(profileset) == 2
        assert profileset.get_profile_weight(frozenset({1, 2, 3, 4})) == 2.0
        assert profileset.get_profile_weight(frozenset({5, 6, 7, 8})) == 1.5
    
    def test_init_from_quartets(self) -> None:
        """Test creating from Quartet objects (grouped into profiles)."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        q3 = Quartet(Split({5, 6}, {7, 8}))
        
        profileset = QuartetProfileSet(profiles=[q1, q2, q3])
        
        assert len(profileset) == 2
        # Profile weight should be sum of quartet weights (both default to 1.0)
        assert profileset.get_profile_weight(frozenset({1, 2, 3, 4})) == 2.0
        assert profileset.get_profile_weight(frozenset({5, 6, 7, 8})) == 1.0
    
    def test_init_from_quartets_with_weights(self) -> None:
        """Test creating from Quartet objects with explicit weights."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        q3 = Quartet(Split({5, 6}, {7, 8}))
        
        profileset = QuartetProfileSet(profiles=[(q1, 0.8), (q2, 0.2), (q3, 1.0)])
        
        assert len(profileset) == 2
        # Profile weight should be sum of quartet weights
        assert profileset.get_profile_weight(frozenset({1, 2, 3, 4})) == 1.0
        assert profileset.get_profile_weight(frozenset({5, 6, 7, 8})) == 1.0
        
        # Check that quartets within profile keep their weights
        profile = profileset.get_profile(frozenset({1, 2, 3, 4}))
        assert profile is not None
        assert profile.get_weight(q1) == 0.8
        assert profile.get_weight(q2) == 0.2
    
    def test_init_mixed_quartet_profiles_error(self) -> None:
        """Test that mixing QuartetProfile and Quartet raises error."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profile1 = QuartetProfile([q1])
        
        with pytest.raises(ValueError, match="Cannot mix QuartetProfile and Quartet"):
            QuartetProfileSet(profiles=[profile1, q1])
    
    def test_init_non_positive_profile_weight_error(self) -> None:
        """Test that non-positive profile weights raise error."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profile1 = QuartetProfile([q1])
        
        with pytest.raises(ValueError, match="Profile weight must be positive"):
            QuartetProfileSet(profiles=[(profile1, 0.0)])
        
        with pytest.raises(ValueError, match="Profile weight must be positive"):
            QuartetProfileSet(profiles=[(profile1, -0.5)])
    
    def test_init_non_positive_quartet_weight_error(self) -> None:
        """Test that non-positive quartet weights raise error."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        
        with pytest.raises(ValueError, match="Quartet weight must be positive"):
            QuartetProfileSet(profiles=[(q1, 0.0)])
        
        with pytest.raises(ValueError, match="Quartet weight must be positive"):
            QuartetProfileSet(profiles=[(q1, -0.5)])
    
    def test_init_with_taxa_parameter(self) -> None:
        """Test initialization with explicit taxa parameter."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))
        
        profileset = QuartetProfileSet(
            profiles=[q1, q2],
            taxa=frozenset({1, 2, 3, 4, 5, 6, 7, 8, 9, 10})
        )
        
        assert len(profileset.taxa) == 10
        assert 9 in profileset.taxa
        assert 10 in profileset.taxa
    
    def test_init_with_taxa_not_superset_error(self) -> None:
        """Test that taxa parameter must be a superset."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        
        with pytest.raises(ValueError, match="Provided taxa must be a superset"):
            QuartetProfileSet(profiles=[q1], taxa=frozenset({1, 2, 3}))
    
    def test_init_duplicate_taxa_profiles(self) -> None:
        """Test that duplicate taxa in profiles overwrite previous profile."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        profile1 = QuartetProfile([q1])
        profile2 = QuartetProfile([q2])
        
        profileset = QuartetProfileSet(profiles=[(profile1, 1.0), (profile2, 2.0)])
        
        # Should only have one profile (last one wins)
        assert len(profileset) == 1
        assert profileset.get_profile_weight(frozenset({1, 2, 3, 4})) == 2.0
        # Should have q2, not q1
        profile = profileset.get_profile(frozenset({1, 2, 3, 4}))
        assert profile is not None
        assert q2 in profile
        assert q1 not in profile


class TestQuartetProfileSetProperties:
    """Tests for QuartetProfileSet properties."""
    
    def test_profiles_property(self) -> None:
        """Test profiles property."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profile1 = QuartetProfile([q1])
        
        profileset = QuartetProfileSet(profiles=[profile1])
        
        profiles = profileset.profiles
        assert isinstance(profiles, Mapping)
        assert frozenset({1, 2, 3, 4}) in profiles
        profile, weight = profiles[frozenset({1, 2, 3, 4})]
        assert isinstance(profile, QuartetProfile)
        assert weight == 1.0
    
    def test_taxa_property(self) -> None:
        """Test taxa property."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))
        
        profileset = QuartetProfileSet(profiles=[q1, q2])
        
        assert profileset.taxa == frozenset({1, 2, 3, 4, 5, 6, 7, 8})
        assert isinstance(profileset.taxa, frozenset)
    
    def test_taxa_property_with_explicit_taxa(self) -> None:
        """Test taxa property when explicit taxa is provided."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        
        profileset = QuartetProfileSet(
            profiles=[q1],
            taxa=frozenset({1, 2, 3, 4, 5, 6, 7, 8})
        )
        
        assert len(profileset.taxa) == 8
        assert 5 in profileset.taxa
        assert 6 in profileset.taxa
    
    def test_is_dense_empty(self) -> None:
        """Test is_dense for empty profile set."""
        profileset = QuartetProfileSet()
        
        assert profileset.is_dense is True
    
    def test_is_dense_single_profile(self) -> None:
        """Test is_dense for single profile."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profileset = QuartetProfileSet(profiles=[q1])
        
        # 4 taxa, should have C(4,4) = 1 profile for dense
        assert profileset.is_dense is True
    
    def test_is_dense_not_dense(self) -> None:
        """Test is_dense for non-dense profile set."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))
        profileset = QuartetProfileSet(profiles=[q1, q2])
        
        # 8 taxa, should have C(8,4) = 70 profiles for dense, but only 2
        assert profileset.is_dense is False
    
    def test_is_dense_with_explicit_taxa(self) -> None:
        """Test is_dense when explicit taxa is provided."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profileset = QuartetProfileSet(
            profiles=[q1],
            taxa=frozenset({1, 2, 3, 4, 5, 6, 7, 8})
        )
        
        # 8 taxa, should have C(8,4) = 70 profiles for dense, but only 1
        assert profileset.is_dense is False
    
    def test_is_all_resolved_all_resolved(self) -> None:
        """Test is_all_resolved when all profiles are resolved."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        q3 = Quartet(Split({5, 6}, {7, 8}))
        profileset = QuartetProfileSet(profiles=[q1, q2, q3])
        
        assert profileset.is_all_resolved is True
    
    def test_is_all_resolved_mixed(self) -> None:
        """Test is_all_resolved when some profiles are not resolved."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        star = Quartet({1, 2, 3, 4})
        q2 = Quartet(Split({5, 6}, {7, 8}))
        profileset = QuartetProfileSet(profiles=[q1, star, q2])
        
        assert profileset.is_all_resolved is False
    
    def test_is_all_resolved_all_star(self) -> None:
        """Test is_all_resolved when all profiles are star trees."""
        star1 = Quartet({1, 2, 3, 4})
        star2 = Quartet({5, 6, 7, 8})
        profileset = QuartetProfileSet(profiles=[star1, star2])
        
        assert profileset.is_all_resolved is False
    
    def test_is_all_resolved_empty(self) -> None:
        """Test is_all_resolved for empty profile set."""
        profileset = QuartetProfileSet()
        
        assert profileset.is_all_resolved is True
    
    def test_is_all_resolved_profile_with_multiple_quartets(self) -> None:
        """Test is_all_resolved with profiles containing multiple resolved quartets."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        q3 = Quartet(Split({1, 4}, {2, 3}))
        profile = QuartetProfile([q1, q2, q3])
        profileset = QuartetProfileSet(profiles=[profile])
        
        assert profile.is_resolved() is True
        assert profileset.is_all_resolved is True
    
    def test_max_profile_len_single_quartet(self) -> None:
        """Test max_profile_len with profiles containing single quartets."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))
        profileset = QuartetProfileSet(profiles=[q1, q2])
        
        assert profileset.max_profile_len == 1
    
    def test_max_profile_len_multiple_quartets(self) -> None:
        """Test max_profile_len with profiles containing multiple quartets."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        q3 = Quartet(Split({1, 4}, {2, 3}))
        profile1 = QuartetProfile([q1, q2, q3])
        profile2 = QuartetProfile([Quartet(Split({5, 6}, {7, 8}))])
        profileset = QuartetProfileSet(profiles=[profile1, profile2])
        
        assert profileset.max_profile_len == 3
    
    def test_max_profile_len_empty(self) -> None:
        """Test max_profile_len for empty profile set."""
        profileset = QuartetProfileSet()
        
        assert profileset.max_profile_len == 0
    
    def test_max_profile_len_mixed_lengths(self) -> None:
        """Test max_profile_len with profiles of different lengths."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        profile1 = QuartetProfile([q1])  # length 1
        profile2 = QuartetProfile([q1, q2])  # length 2
        q3 = Quartet(Split({5, 6}, {7, 8}))
        profile3 = QuartetProfile([q3])  # length 1
        profileset = QuartetProfileSet(profiles=[profile1, profile2, profile3])
        
        assert profileset.max_profile_len == 2
    
    def test_max_profile_len_with_star_trees(self) -> None:
        """Test max_profile_len with profiles containing star trees."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        star = Quartet({1, 2, 3, 4})
        profile = QuartetProfile([q1, star])
        profileset = QuartetProfileSet(profiles=[profile])
        
        assert profileset.max_profile_len == 2


class TestQuartetProfileSetMethods:
    """Tests for QuartetProfileSet methods."""
    
    def test_get_profile_existing(self) -> None:
        """Test get_profile for existing profile."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profile1 = QuartetProfile([q1])
        
        profileset = QuartetProfileSet(profiles=[profile1])
        
        profile = profileset.get_profile(frozenset({1, 2, 3, 4}))
        assert profile is not None
        assert isinstance(profile, QuartetProfile)
        assert q1 in profile
    
    def test_get_profile_not_existing(self) -> None:
        """Test get_profile for non-existing profile."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profileset = QuartetProfileSet(profiles=[q1])
        
        profile = profileset.get_profile(frozenset({5, 6, 7, 8}))
        assert profile is None
    
    def test_get_profile_weight_existing(self) -> None:
        """Test get_profile_weight for existing profile."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profile1 = QuartetProfile([q1])
        
        profileset = QuartetProfileSet(profiles=[(profile1, 2.5)])
        
        weight = profileset.get_profile_weight(frozenset({1, 2, 3, 4}))
        assert weight == 2.5
    
    def test_get_profile_weight_not_existing(self) -> None:
        """Test get_profile_weight for non-existing profile."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profileset = QuartetProfileSet(profiles=[q1])
        
        weight = profileset.get_profile_weight(frozenset({5, 6, 7, 8}))
        assert weight is None
    
    def test_has_profile_existing(self) -> None:
        """Test has_profile for existing profile."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profileset = QuartetProfileSet(profiles=[q1])
        
        assert profileset.has_profile(frozenset({1, 2, 3, 4})) is True
    
    def test_has_profile_not_existing(self) -> None:
        """Test has_profile for non-existing profile."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profileset = QuartetProfileSet(profiles=[q1])
        
        assert profileset.has_profile(frozenset({5, 6, 7, 8})) is False
    
    def test_all_profile_taxon_sets(self) -> None:
        """Test all_profile_taxon_sets generator."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))
        profileset = QuartetProfileSet(profiles=[q1, q2])
        
        taxa_sets = list(profileset.all_profile_taxon_sets())
        assert len(taxa_sets) == 2
        assert frozenset({1, 2, 3, 4}) in taxa_sets
        assert frozenset({5, 6, 7, 8}) in taxa_sets
    
    def test_all_profile_taxon_sets_empty(self) -> None:
        """Test all_profile_taxon_sets for empty profile set."""
        profileset = QuartetProfileSet()
        
        taxa_sets = list(profileset.all_profile_taxon_sets())
        assert len(taxa_sets) == 0


class TestQuartetProfileSetMagicMethods:
    """Tests for QuartetProfileSet magic methods."""
    
    def test_len(self) -> None:
        """Test __len__ method."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))
        profileset = QuartetProfileSet(profiles=[q1, q2])
        
        assert len(profileset) == 2
    
    def test_len_empty(self) -> None:
        """Test __len__ for empty profile set."""
        profileset = QuartetProfileSet()
        
        assert len(profileset) == 0
    
    def test_iter(self) -> None:
        """Test __iter__ method."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))
        profileset = QuartetProfileSet(profiles=[q1, q2])
        
        items = list(profileset)
        assert len(items) == 2
        for profile, weight in items:
            assert isinstance(profile, QuartetProfile)
            assert isinstance(weight, float)
            assert weight > 0
    
    def test_contains_existing(self) -> None:
        """Test __contains__ for existing profile."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profileset = QuartetProfileSet(profiles=[q1])
        
        assert frozenset({1, 2, 3, 4}) in profileset
    
    def test_contains_not_existing(self) -> None:
        """Test __contains__ for non-existing profile."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profileset = QuartetProfileSet(profiles=[q1])
        
        assert frozenset({5, 6, 7, 8}) not in profileset
    
    def test_repr_empty(self) -> None:
        """Test __repr__ for empty profile set."""
        profileset = QuartetProfileSet()
        
        repr_str = repr(profileset)
        assert "QuartetProfileSet" in repr_str
        assert "profiles={}" in repr_str
    
    def test_repr_with_profiles(self) -> None:
        """Test __repr__ with profiles."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profileset = QuartetProfileSet(profiles=[q1])
        
        repr_str = repr(profileset)
        assert "QuartetProfileSet" in repr_str
        assert "profiles=" in repr_str


class TestQuartetProfileSetEdgeCases:
    """Tests for QuartetProfileSet edge cases."""
    
    def test_single_profile(self) -> None:
        """Test profile set with single profile."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profileset = QuartetProfileSet(profiles=[q1])
        
        assert len(profileset) == 1
        assert len(profileset.taxa) == 4
        assert profileset.is_dense is True
    
    def test_multiple_quartets_same_taxa(self) -> None:
        """Test multiple quartets on same 4-taxon set."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        q3 = Quartet(Split({1, 4}, {2, 3}))
        
        profileset = QuartetProfileSet(profiles=[(q1, 0.5), (q2, 0.3), (q3, 0.2)])
        
        assert len(profileset) == 1
        profile = profileset.get_profile(frozenset({1, 2, 3, 4}))
        assert profile is not None
        assert len(profile) == 3
        assert profileset.get_profile_weight(frozenset({1, 2, 3, 4})) == 1.0
    
    def test_star_tree_quartets(self) -> None:
        """Test profile set with star tree quartets."""
        star1 = Quartet({1, 2, 3, 4})
        star2 = Quartet({5, 6, 7, 8})
        
        profileset = QuartetProfileSet(profiles=[star1, star2])
        
        assert len(profileset) == 2
        assert star1 in profileset.get_profile(frozenset({1, 2, 3, 4}))
        assert star2 in profileset.get_profile(frozenset({5, 6, 7, 8}))
    
    def test_mixed_resolved_and_star(self) -> None:
        """Test profile set with both resolved and star quartets."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        star = Quartet({1, 2, 3, 4})
        
        profileset = QuartetProfileSet(profiles=[q1, star])
        
        assert len(profileset) == 1
        profile = profileset.get_profile(frozenset({1, 2, 3, 4}))
        assert profile is not None
        assert len(profile) == 2
        assert q1 in profile
        assert star in profile
    
    def test_large_profile_set(self) -> None:
        """Test profile set with many profiles."""
        quartets = []
        for i in range(10):
            taxa_start = i * 4
            q = Quartet(Split(
                {taxa_start + 1, taxa_start + 2},
                {taxa_start + 3, taxa_start + 4}
            ))
            quartets.append(q)
        
        profileset = QuartetProfileSet(profiles=quartets)
        
        assert len(profileset) == 10
        assert len(profileset.taxa) == 40
    
    def test_profile_weights_sum_quartet_weights(self) -> None:
        """Test that profile weight equals sum of quartet weights."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({1, 3}, {2, 4}))
        q3 = Quartet(Split({1, 4}, {2, 3}))
        
        profileset = QuartetProfileSet(profiles=[
            (q1, 0.4),
            (q2, 0.3),
            (q3, 0.3)
        ])
        
        # Profile weight should be sum: 0.4 + 0.3 + 0.3 = 1.0
        assert profileset.get_profile_weight(frozenset({1, 2, 3, 4})) == 1.0
        
        # Individual quartet weights should be preserved
        profile = profileset.get_profile(frozenset({1, 2, 3, 4}))
        assert profile is not None
        assert profile.get_weight(q1) == 0.4
        assert profile.get_weight(q2) == 0.3
        assert profile.get_weight(q3) == 0.3


class TestQuartetProfileSetValidation:
    """Tests for QuartetProfileSet validation."""
    
    def test_invalid_type_in_profiles_error(self) -> None:
        """Test that invalid types in profiles raise error."""
        # When no QuartetProfile is detected, it tries Quartet mode
        with pytest.raises(ValueError, match="Expected Quartet"):
            QuartetProfileSet(profiles=["not a profile"])
        
        with pytest.raises(ValueError, match="Expected Quartet"):
            QuartetProfileSet(profiles=[("not a quartet", 1.0)])
        
        # Test with actual QuartetProfile to trigger QuartetProfile validation
        q1 = Quartet(Split({1, 2}, {3, 4}))
        profile1 = QuartetProfile([q1])
        
        # Mixing with invalid type should fail
        with pytest.raises(ValueError, match="Expected QuartetProfile"):
            QuartetProfileSet(profiles=[profile1, "not a profile"])
    
    def test_empty_profiles_list(self) -> None:
        """Test that empty profiles list creates empty profile set."""
        profileset = QuartetProfileSet(profiles=[])
        
        assert len(profileset) == 0
        assert len(profileset.taxa) == 0
    
    def test_taxa_superset_validation(self) -> None:
        """Test taxa superset validation with multiple profiles."""
        q1 = Quartet(Split({1, 2}, {3, 4}))
        q2 = Quartet(Split({5, 6}, {7, 8}))
        
        # Valid: taxa is superset
        profileset1 = QuartetProfileSet(
            profiles=[q1, q2],
            taxa=frozenset({1, 2, 3, 4, 5, 6, 7, 8, 9, 10})
        )
        assert len(profileset1.taxa) == 10
        
        # Invalid: taxa is not superset
        with pytest.raises(ValueError, match="Provided taxa must be a superset"):
            QuartetProfileSet(
                profiles=[q1, q2],
                taxa=frozenset({1, 2, 3, 4})
            )

