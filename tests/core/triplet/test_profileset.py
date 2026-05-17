"""
Tests for TripletProfileSet class.
"""

from typing import Mapping

import pytest

from phylozoo.core.triplet import Triplet, TripletProfile, TripletProfileSet
from phylozoo.core.split.base import Split


class TestTripletProfileSetInit:
    """Tests for TripletProfileSet initialization."""

    def test_init_empty(self) -> None:
        """Test creating an empty profile set."""
        profileset = TripletProfileSet()

        assert len(profileset) == 0
        assert len(profileset.taxa) == 0
        assert profileset.is_dense is True

    def test_init_from_triplet_profiles(self) -> None:
        """Test creating from TripletProfile objects."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        profile1 = TripletProfile({t1: 0.8, t2: 0.2})
        profile2 = TripletProfile([Triplet(Split({4}, {5, 6}))])

        profileset = TripletProfileSet(profiles=[profile1, profile2])

        assert len(profileset) == 2
        assert profileset.get_profile_weight(frozenset({1, 2, 3})) == 1.0
        assert profileset.get_profile_weight(frozenset({4, 5, 6})) == 1.0

    def test_init_from_triplet_profiles_with_weights(self) -> None:
        """Test creating from TripletProfile objects with explicit weights."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        profile1 = TripletProfile({t1: 0.8, t2: 0.2})
        profile2 = TripletProfile([Triplet(Split({4}, {5, 6}))])

        profileset = TripletProfileSet(profiles=[(profile1, 2.0), (profile2, 1.5)])

        assert len(profileset) == 2
        assert profileset.get_profile_weight(frozenset({1, 2, 3})) == 2.0
        assert profileset.get_profile_weight(frozenset({4, 5, 6})) == 1.5

    def test_init_from_triplets(self) -> None:
        """Test creating from Triplet objects (grouped into profiles)."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        t3 = Triplet(Split({4}, {5, 6}))

        profileset = TripletProfileSet(profiles=[t1, t2, t3])

        assert len(profileset) == 2
        # Each profile constructed from bare triplets gets default profile weight 1.0
        assert profileset.get_profile_weight(frozenset({1, 2, 3})) == 1.0
        assert profileset.get_profile_weight(frozenset({4, 5, 6})) == 1.0

    def test_init_from_triplets_with_weights_not_supported(self) -> None:
        """Test that passing triplets with explicit weights is not supported."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        t3 = Triplet(Split({4}, {5, 6}))

        with pytest.raises(ValueError, match="Triplet weights are not supported"):
            TripletProfileSet(profiles=[(t1, 0.8), (t2, 0.2), (t3, 1.0)])

    def test_init_mixed_triplet_profiles_error(self) -> None:
        """Test that mixing TripletProfile and Triplet raises error."""
        t1 = Triplet(Split({1}, {2, 3}))
        profile1 = TripletProfile([t1])

        with pytest.raises(ValueError, match="Cannot mix TripletProfile and Triplet"):
            TripletProfileSet(profiles=[profile1, t1])

    def test_init_non_positive_profile_weight_error(self) -> None:
        """Test that non-positive profile weights raise error."""
        t1 = Triplet(Split({1}, {2, 3}))
        profile1 = TripletProfile([t1])

        with pytest.raises(ValueError, match="Profile weight must be positive"):
            TripletProfileSet(profiles=[(profile1, 0.0)])

        with pytest.raises(ValueError, match="Profile weight must be positive"):
            TripletProfileSet(profiles=[(profile1, -0.5)])

    def test_init_duplicate_profile_taxa_error(self) -> None:
        """Test that multiple profiles with same taxa raise error."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        profile1 = TripletProfile([t1])
        profile2 = TripletProfile([t2])  # Same taxa as profile1

        with pytest.raises(ValueError, match="Multiple profiles with the same taxa set"):
            TripletProfileSet(profiles=[profile1, profile2])

        with pytest.raises(ValueError, match="Multiple profiles with the same taxa set"):
            TripletProfileSet(profiles=[(profile1, 1.0), (profile2, 2.0)])

    def test_init_duplicate_triplet_in_triplet_mode_error(self) -> None:
        """Test that duplicate triplets in triplet mode raise error."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))

        # Same triplet appears twice (same taxa)
        with pytest.raises(ValueError, match="appears multiple times in the input"):
            TripletProfileSet(profiles=[t1, t2, t1])

        # Different triplets with same taxa (this is OK - they get merged)
        t3 = Triplet(Split({3}, {1, 2}))  # Different triplet, same taxa
        profileset = TripletProfileSet(profiles=[t1, t2, t3])
        # Should work - all three triplets merged into one profile
        profile = profileset.get_profile(frozenset({1, 2, 3}))
        assert profile is not None
        assert len(profile) == 3

    def test_init_triplet_weights_not_supported(self) -> None:
        """Test that providing triplet weights is not supported at all."""
        t1 = Triplet(Split({1}, {2, 3}))

        with pytest.raises(ValueError, match="Triplet weights are not supported"):
            TripletProfileSet(profiles=[(t1, 0.0)])

        with pytest.raises(ValueError, match="Triplet weights are not supported"):
            TripletProfileSet(profiles=[(t1, -0.5)])

    def test_init_with_taxa_parameter(self) -> None:
        """Test initialization with explicit taxa parameter."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({4}, {5, 6}))

        profileset = TripletProfileSet(profiles=[t1, t2], taxa=frozenset({1, 2, 3, 4, 5, 6, 7, 8}))

        assert len(profileset.taxa) == 8
        assert 7 in profileset.taxa
        assert 8 in profileset.taxa

    def test_init_with_taxa_not_superset_error(self) -> None:
        """Test that taxa parameter must be a superset."""
        t1 = Triplet(Split({1}, {2, 3}))

        with pytest.raises(ValueError, match="Provided taxa must be a superset"):
            TripletProfileSet(profiles=[t1], taxa=frozenset({1, 2}))

    def test_init_duplicate_taxa_profiles(self) -> None:
        """Test that duplicate taxa in profiles raise an error."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        profile1 = TripletProfile([t1])
        profile2 = TripletProfile([t2])  # Same taxa as profile1

        # Should raise error for duplicate taxa sets
        with pytest.raises(ValueError, match="Multiple profiles with the same taxa set"):
            TripletProfileSet(profiles=[(profile1, 1.0), (profile2, 2.0)])


class TestTripletProfileSetProperties:
    """Tests for TripletProfileSet properties."""

    def test_profiles_property(self) -> None:
        """Test profiles property."""
        t1 = Triplet(Split({1}, {2, 3}))
        profile1 = TripletProfile([t1])

        profileset = TripletProfileSet(profiles=[profile1])

        profiles = profileset.profiles
        assert isinstance(profiles, Mapping)
        assert frozenset({1, 2, 3}) in profiles
        profile, weight = profiles[frozenset({1, 2, 3})]
        assert isinstance(profile, TripletProfile)
        assert weight == 1.0

    def test_taxa_property(self) -> None:
        """Test taxa property."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({4}, {5, 6}))

        profileset = TripletProfileSet(profiles=[t1, t2])

        assert profileset.taxa == frozenset({1, 2, 3, 4, 5, 6})
        assert isinstance(profileset.taxa, frozenset)

    def test_taxa_property_with_explicit_taxa(self) -> None:
        """Test taxa property when explicit taxa is provided."""
        t1 = Triplet(Split({1}, {2, 3}))

        profileset = TripletProfileSet(profiles=[t1], taxa=frozenset({1, 2, 3, 4, 5, 6}))

        assert len(profileset.taxa) == 6
        assert 4 in profileset.taxa
        assert 5 in profileset.taxa

    def test_is_dense_empty(self) -> None:
        """Test is_dense for empty profile set."""
        profileset = TripletProfileSet()

        assert profileset.is_dense is True

    def test_is_dense_single_profile(self) -> None:
        """Test is_dense for single profile."""
        t1 = Triplet(Split({1}, {2, 3}))
        profileset = TripletProfileSet(profiles=[t1])

        # 3 taxa, should have C(3,3) = 1 profile for dense
        assert profileset.is_dense is True

    def test_is_dense_not_dense(self) -> None:
        """Test is_dense for non-dense profile set."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({4}, {5, 6}))
        profileset = TripletProfileSet(profiles=[t1, t2])

        # 6 taxa, should have C(6,3) = 20 profiles for dense, but only 2
        assert profileset.is_dense is False

    def test_is_dense_with_explicit_taxa(self) -> None:
        """Test is_dense when explicit taxa is provided."""
        t1 = Triplet(Split({1}, {2, 3}))
        profileset = TripletProfileSet(profiles=[t1], taxa=frozenset({1, 2, 3, 4, 5, 6}))

        # 6 taxa, should have C(6,3) = 20 profiles for dense, but only 1
        assert profileset.is_dense is False

    def test_is_all_resolved_all_resolved(self) -> None:
        """Test is_all_resolved when all profiles are resolved."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        t3 = Triplet(Split({4}, {5, 6}))
        profileset = TripletProfileSet(profiles=[t1, t2, t3])

        assert profileset.is_all_resolved is True

    def test_is_all_resolved_mixed(self) -> None:
        """Test is_all_resolved when some profiles are not resolved."""
        t1 = Triplet(Split({1}, {2, 3}))
        star = Triplet({1, 2, 3})
        t2 = Triplet(Split({4}, {5, 6}))
        profileset = TripletProfileSet(profiles=[t1, star, t2])

        assert profileset.is_all_resolved is False

    def test_is_all_resolved_all_star(self) -> None:
        """Test is_all_resolved when all profiles are star trees."""
        star1 = Triplet({1, 2, 3})
        star2 = Triplet({4, 5, 6})
        profileset = TripletProfileSet(profiles=[star1, star2])

        assert profileset.is_all_resolved is False

    def test_is_all_resolved_empty(self) -> None:
        """Test is_all_resolved for empty profile set."""
        profileset = TripletProfileSet()

        assert profileset.is_all_resolved is True

    def test_is_all_resolved_profile_with_multiple_triplets(self) -> None:
        """Test is_all_resolved with profiles containing multiple resolved triplets."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        t3 = Triplet(Split({3}, {1, 2}))
        profile = TripletProfile([t1, t2, t3])
        profileset = TripletProfileSet(profiles=[profile])

        assert profile.is_resolved() is True
        assert profileset.is_all_resolved is True

    def test_max_profile_len_single_triplet(self) -> None:
        """Test max_profile_len with profiles containing single triplets."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({4}, {5, 6}))
        profileset = TripletProfileSet(profiles=[t1, t2])

        assert profileset.max_profile_len == 1

    def test_max_profile_len_multiple_triplets(self) -> None:
        """Test max_profile_len with profiles containing multiple triplets."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        t3 = Triplet(Split({3}, {1, 2}))
        profile1 = TripletProfile([t1, t2, t3])
        profile2 = TripletProfile([Triplet(Split({4}, {5, 6}))])
        profileset = TripletProfileSet(profiles=[profile1, profile2])

        assert profileset.max_profile_len == 3

    def test_max_profile_len_empty(self) -> None:
        """Test max_profile_len for empty profile set."""
        profileset = TripletProfileSet()

        assert profileset.max_profile_len == 0

    def test_max_profile_len_mixed_lengths(self) -> None:
        """Test max_profile_len with profiles of different lengths."""
        t1 = Triplet(Split({1}, {2, 3}))
        profile1 = TripletProfile([t1])  # length 1
        t3 = Triplet(Split({4}, {5, 6}))
        t4 = Triplet(Split({5}, {4, 6}))
        profile2 = TripletProfile([t3, t4])  # length 2, different taxa
        t5 = Triplet(Split({7}, {8, 9}))
        profile3 = TripletProfile([t5])  # length 1, different taxa
        profileset = TripletProfileSet(profiles=[profile1, profile2, profile3])

        assert profileset.max_profile_len == 2

    def test_max_profile_len_with_star_trees(self) -> None:
        """Test max_profile_len with profiles containing star trees."""
        t1 = Triplet(Split({1}, {2, 3}))
        star = Triplet({1, 2, 3})
        profile = TripletProfile([t1, star])
        profileset = TripletProfileSet(profiles=[profile])

        assert profileset.max_profile_len == 2


class TestTripletProfileSetMethods:
    """Tests for TripletProfileSet methods."""

    def test_get_profile_existing(self) -> None:
        """Test get_profile for existing profile."""
        t1 = Triplet(Split({1}, {2, 3}))
        profile1 = TripletProfile([t1])

        profileset = TripletProfileSet(profiles=[profile1])

        profile = profileset.get_profile(frozenset({1, 2, 3}))
        assert profile is not None
        assert isinstance(profile, TripletProfile)
        assert t1 in profile

    def test_get_profile_not_existing(self) -> None:
        """Test get_profile for non-existing profile."""
        t1 = Triplet(Split({1}, {2, 3}))
        profileset = TripletProfileSet(profiles=[t1])

        profile = profileset.get_profile(frozenset({4, 5, 6}))
        assert profile is None

    def test_get_profile_weight_existing(self) -> None:
        """Test get_profile_weight for existing profile."""
        t1 = Triplet(Split({1}, {2, 3}))
        profile1 = TripletProfile([t1])

        profileset = TripletProfileSet(profiles=[(profile1, 2.5)])

        weight = profileset.get_profile_weight(frozenset({1, 2, 3}))
        assert weight == 2.5

    def test_get_profile_weight_not_existing(self) -> None:
        """Test get_profile_weight for non-existing profile."""
        t1 = Triplet(Split({1}, {2, 3}))
        profileset = TripletProfileSet(profiles=[t1])

        weight = profileset.get_profile_weight(frozenset({4, 5, 6}))
        assert weight is None

    def test_has_profile_existing(self) -> None:
        """Test has_profile for existing profile."""
        t1 = Triplet(Split({1}, {2, 3}))
        profileset = TripletProfileSet(profiles=[t1])

        assert profileset.has_profile(frozenset({1, 2, 3})) is True

    def test_has_profile_not_existing(self) -> None:
        """Test has_profile for non-existing profile."""
        t1 = Triplet(Split({1}, {2, 3}))
        profileset = TripletProfileSet(profiles=[t1])

        assert profileset.has_profile(frozenset({4, 5, 6})) is False

    def test_all_profile_taxon_sets(self) -> None:
        """Test all_profile_taxon_sets generator."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({4}, {5, 6}))
        profileset = TripletProfileSet(profiles=[t1, t2])

        taxa_sets = list(profileset.all_profile_taxon_sets())
        assert len(taxa_sets) == 2
        assert frozenset({1, 2, 3}) in taxa_sets
        assert frozenset({4, 5, 6}) in taxa_sets

    def test_all_profile_taxon_sets_empty(self) -> None:
        """Test all_profile_taxon_sets for empty profile set."""
        profileset = TripletProfileSet()

        taxa_sets = list(profileset.all_profile_taxon_sets())
        assert len(taxa_sets) == 0


class TestTripletProfileSetMagicMethods:
    """Tests for TripletProfileSet magic methods."""

    def test_len(self) -> None:
        """Test __len__ method."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({4}, {5, 6}))
        profileset = TripletProfileSet(profiles=[t1, t2])

        assert len(profileset) == 2

    def test_len_empty(self) -> None:
        """Test __len__ for empty profile set."""
        profileset = TripletProfileSet()

        assert len(profileset) == 0

    def test_iter(self) -> None:
        """Test __iter__ method."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({4}, {5, 6}))
        profileset = TripletProfileSet(profiles=[t1, t2])

        items = list(profileset)
        assert len(items) == 2
        for profile, weight in items:
            assert isinstance(profile, TripletProfile)
            assert isinstance(weight, float)
            assert weight > 0

    def test_contains_existing(self) -> None:
        """Test __contains__ for existing profile."""
        t1 = Triplet(Split({1}, {2, 3}))
        profileset = TripletProfileSet(profiles=[t1])

        assert frozenset({1, 2, 3}) in profileset

    def test_contains_not_existing(self) -> None:
        """Test __contains__ for non-existing profile."""
        t1 = Triplet(Split({1}, {2, 3}))
        profileset = TripletProfileSet(profiles=[t1])

        assert frozenset({4, 5, 6}) not in profileset

    def test_repr_empty(self) -> None:
        """Test __repr__ for empty profile set."""
        profileset = TripletProfileSet()

        repr_str = repr(profileset)
        assert "TripletProfileSet" in repr_str
        assert "profiles={}" in repr_str

    def test_repr_with_profiles(self) -> None:
        """Test __repr__ with profiles."""
        t1 = Triplet(Split({1}, {2, 3}))
        profileset = TripletProfileSet(profiles=[t1])

        repr_str = repr(profileset)
        assert "TripletProfileSet" in repr_str
        assert "profiles=" in repr_str


class TestTripletProfileSetEdgeCases:
    """Tests for TripletProfileSet edge cases."""

    def test_single_profile(self) -> None:
        """Test profile set with single profile."""
        t1 = Triplet(Split({1}, {2, 3}))
        profileset = TripletProfileSet(profiles=[t1])

        assert len(profileset) == 1
        assert len(profileset.taxa) == 3
        assert profileset.is_dense is True

    def test_multiple_triplets_same_taxa(self) -> None:
        """Test multiple triplets on same 3-taxon set."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        t3 = Triplet(Split({3}, {1, 2}))

        profileset = TripletProfileSet(profiles=[t1, t2, t3])

        assert len(profileset) == 1
        profile = profileset.get_profile(frozenset({1, 2, 3}))
        assert profile is not None
        assert len(profile) == 3
        # Profile built from bare triplets has default weight 1.0
        assert profileset.get_profile_weight(frozenset({1, 2, 3})) == 1.0

    def test_star_tree_triplets(self) -> None:
        """Test profile set with star tree triplets."""
        star1 = Triplet({1, 2, 3})
        star2 = Triplet({4, 5, 6})

        profileset = TripletProfileSet(profiles=[star1, star2])

        assert len(profileset) == 2
        assert star1 in profileset.get_profile(frozenset({1, 2, 3}))
        assert star2 in profileset.get_profile(frozenset({4, 5, 6}))

    def test_mixed_resolved_and_star(self) -> None:
        """Test profile set with both resolved and star triplets."""
        t1 = Triplet(Split({1}, {2, 3}))
        star = Triplet({1, 2, 3})

        profileset = TripletProfileSet(profiles=[t1, star])

        assert len(profileset) == 1
        profile = profileset.get_profile(frozenset({1, 2, 3}))
        assert profile is not None
        assert len(profile) == 2
        assert t1 in profile
        assert star in profile

    def test_large_profile_set(self) -> None:
        """Test profile set with many profiles."""
        triplets = []
        for i in range(10):
            taxa_start = i * 3
            t = Triplet(Split({taxa_start + 1}, {taxa_start + 2, taxa_start + 3}))
            triplets.append(t)

        profileset = TripletProfileSet(profiles=triplets)

        assert len(profileset) == 10
        assert len(profileset.taxa) == 30

    def test_profile_weights_and_triplet_weights_from_profile(self) -> None:
        """Test profile and triplet weights when constructed from TripletProfile."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        t3 = Triplet(Split({3}, {1, 2}))

        profile = TripletProfile({t1: 0.4, t2: 0.3, t3: 0.3})
        profileset = TripletProfileSet(profiles=[profile])

        # Default profile weight is 1.0
        assert profileset.get_profile_weight(frozenset({1, 2, 3})) == 1.0

        # Individual triplet weights are preserved from the profile
        stored_profile = profileset.get_profile(frozenset({1, 2, 3}))
        assert stored_profile is not None
        assert stored_profile.get_weight(t1) == 0.4
        assert stored_profile.get_weight(t2) == 0.3
        assert stored_profile.get_weight(t3) == 0.3


class TestTripletProfileSetValidation:
    """Tests for TripletProfileSet validation."""

    def test_invalid_type_in_profiles_error(self) -> None:
        """Test that invalid types in profiles raise error."""
        # When no TripletProfile is detected, it tries Triplet mode
        with pytest.raises(ValueError, match="Expected Triplet"):
            TripletProfileSet(profiles=["not a profile"])

        with pytest.raises(ValueError, match="Expected Triplet"):
            TripletProfileSet(profiles=[("not a triplet", 1.0)])

        # Test with actual TripletProfile to trigger TripletProfile validation
        t1 = Triplet(Split({1}, {2, 3}))
        profile1 = TripletProfile([t1])

        # Mixing with invalid type should fail
        with pytest.raises(ValueError, match="Expected TripletProfile"):
            TripletProfileSet(profiles=[profile1, "not a profile"])

    def test_empty_profiles_list(self) -> None:
        """Test that empty profiles list creates empty profile set."""
        profileset = TripletProfileSet(profiles=[])

        assert len(profileset) == 0
        assert len(profileset.taxa) == 0

    def test_taxa_superset_validation(self) -> None:
        """Test taxa superset validation with multiple profiles."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({4}, {5, 6}))

        # Valid: taxa is superset
        profileset1 = TripletProfileSet(profiles=[t1, t2], taxa=frozenset({1, 2, 3, 4, 5, 6, 7, 8}))
        assert len(profileset1.taxa) == 8

        # Invalid: taxa is not superset
        with pytest.raises(ValueError, match="Provided taxa must be a superset"):
            TripletProfileSet(profiles=[t1, t2], taxa=frozenset({1, 2, 3}))
