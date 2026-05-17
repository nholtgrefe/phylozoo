"""
Tests for TripletProfile class.
"""

from typing import Mapping

import pytest

from phylozoo.core.triplet import Triplet, TripletProfile
from phylozoo.core.split.base import Split


class TestTripletProfileInit:
    """Tests for TripletProfile initialization."""

    def test_init_from_dict(self) -> None:
        """Test creating a profile from a dictionary."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        profile = TripletProfile({t1: 0.8, t2: 0.2})

        assert profile.taxa == frozenset({1, 2, 3})
        assert len(profile) == 2
        assert profile.get_weight(t1) == 0.8
        assert profile.get_weight(t2) == 0.2
        assert abs(sum(profile.triplets.values()) - 1.0) < 1e-9

    def test_init_from_list_triplets(self) -> None:
        """Test creating a profile from a list of triplets (equal weight 1/k each)."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        profile = TripletProfile([t1, t2])

        assert profile.taxa == frozenset({1, 2, 3})
        assert len(profile) == 2
        assert profile.get_weight(t1) == 0.5
        assert profile.get_weight(t2) == 0.5
        assert abs(sum(profile.triplets.values()) - 1.0) < 1e-9

    def test_init_from_list_tuples(self) -> None:
        """Test creating a profile from a list of (triplet, weight) tuples."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        profile = TripletProfile([(t1, 0.7), (t2, 0.3)])

        assert profile.taxa == frozenset({1, 2, 3})
        assert len(profile) == 2
        assert profile.get_weight(t1) == 0.7
        assert profile.get_weight(t2) == 0.3
        assert abs(sum(profile.triplets.values()) - 1.0) < 1e-9

    def test_init_empty_error(self) -> None:
        """Test that empty profile raises ValueError."""
        with pytest.raises(ValueError, match="TripletProfile must have at least one triplet"):
            TripletProfile([])

        with pytest.raises(ValueError, match="TripletProfile must have at least one triplet"):
            TripletProfile({})

    def test_init_different_taxa_error(self) -> None:
        """Test that triplets with different taxa raise ValueError."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({4}, {5, 6}))

        with pytest.raises(ValueError, match="All triplets must have the same taxa"):
            TripletProfile({t1: 0.5, t2: 0.5})

        with pytest.raises(ValueError, match="All triplets must have the same taxa"):
            TripletProfile([t1, t2])

    def test_init_non_positive_weight_error(self) -> None:
        """Test that non-positive weights raise ValueError."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))

        with pytest.raises(ValueError, match="Weight must be positive"):
            TripletProfile({t1: 0.0, t2: 1.0})

        with pytest.raises(ValueError, match="Weight must be positive"):
            TripletProfile({t1: -0.5, t2: 1.0})

        with pytest.raises(ValueError, match="Weight must be positive"):
            TripletProfile([(t1, 0.0), (t2, 1.0)])

    def test_init_duplicate_triplet_error(self) -> None:
        """Test that duplicate triplets raise ValueError."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))

        # Duplicate in list of triplets
        with pytest.raises(ValueError, match="appears multiple times in the input"):
            TripletProfile([t1, t2, t1])

        # Duplicate in list of tuples
        with pytest.raises(ValueError, match="appears multiple times in the input"):
            TripletProfile([(t1, 0.5), (t2, 0.3), (t1, 0.2)])

        # Duplicate key in dict: Python overwrites, so we get one entry per triplet.
        # Weights must sum to 1.0 (no normalization).
        triplets_dict = {t1: 0.4, t2: 0.6}
        triplets_dict[t1] = 0.4  # Overwrite same key; still {t1: 0.4, t2: 0.6}
        profile = TripletProfile(triplets_dict)
        assert profile.get_weight(t1) == pytest.approx(0.4)
        assert profile.get_weight(t2) == pytest.approx(0.6)
        assert abs(sum(profile.triplets.values()) - 1.0) < 1e-9

    def test_init_single_triplet(self) -> None:
        """Test creating a profile with a single triplet (weight 1.0)."""
        t1 = Triplet(Split({1}, {2, 3}))
        profile = TripletProfile([t1])

        assert profile.taxa == frozenset({1, 2, 3})
        assert len(profile) == 1
        assert profile.get_weight(t1) == 1.0
        assert abs(sum(profile.triplets.values()) - 1.0) < 1e-9

    def test_init_mixed_star_and_resolved(self) -> None:
        """Test creating a profile with both star and resolved triplets on same taxa."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet({1, 2, 3})  # Star tree
        profile = TripletProfile({t1: 0.6, t2: 0.4})

        assert profile.taxa == frozenset({1, 2, 3})
        assert len(profile) == 2
        assert profile.get_weight(t1) == 0.6
        assert profile.get_weight(t2) == 0.4

    def test_weights_sum_to_one(self) -> None:
        """Test that triplet weights always sum to 1.0 (within tolerance) for all input types."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        for profile in (
            TripletProfile([t1]),
            TripletProfile([t1, t2]),
            TripletProfile({t1: 0.8, t2: 0.2}),
            TripletProfile([(t1, 0.4), (t2, 0.6)]),
        ):
            assert abs(sum(profile.triplets.values()) - 1.0) < 1e-9


class TestTripletProfileProperties:
    """Tests for TripletProfile properties."""

    def test_taxa_property(self) -> None:
        """Test taxa property."""
        t1 = Triplet(Split({1}, {2, 3}))
        profile = TripletProfile([t1])

        assert profile.taxa == frozenset({1, 2, 3})
        assert isinstance(profile.taxa, frozenset)

    def test_triplets_property(self) -> None:
        """Test triplets property."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        profile = TripletProfile({t1: 0.8, t2: 0.2})

        triplets = profile.triplets
        assert isinstance(triplets, Mapping)
        assert t1 in triplets
        assert t2 in triplets
        assert triplets[t1] == 0.8
        assert triplets[t2] == 0.2

    def test_weights_sum_to_one_property(self) -> None:
        """Test that triplet weights sum to 1.0."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        profile = TripletProfile({t1: 0.3, t2: 0.7})

        assert abs(sum(profile.triplets.values()) - 1.0) < 1e-9

    def test_immutability(self) -> None:
        """Test that profile is immutable."""
        t1 = Triplet(Split({1}, {2, 3}))
        profile = TripletProfile([t1])

        with pytest.raises(AttributeError, match="Cannot modify attribute"):
            profile._taxa = frozenset({5, 6, 7})

        with pytest.raises(AttributeError, match="Cannot modify attribute"):
            profile._triplets = {}

    def test_split_single_resolved(self) -> None:
        """Test split property with single resolved triplet."""
        t1 = Triplet(Split({1}, {2, 3}))
        profile = TripletProfile([t1])

        split = profile.split
        assert split is not None
        # The singleton side is {1}; the cherry is {2, 3}.
        sides = (split.set1, split.set2)
        assert {1} in sides
        assert {2, 3} in sides

    def test_split_single_star(self) -> None:
        """Test split property with single star tree."""
        star = Triplet({1, 2, 3})
        profile = TripletProfile([star])

        assert profile.split is None

    def test_split_multiple_triplets(self) -> None:
        """Test split property with multiple triplets."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        profile = TripletProfile([t1, t2])

        assert profile.split is None

    def test_split_multiple_with_star(self) -> None:
        """Test split property with multiple triplets including star."""
        t1 = Triplet(Split({1}, {2, 3}))
        star = Triplet({1, 2, 3})
        profile = TripletProfile([t1, star])

        assert profile.split is None

    def test_split_cached(self) -> None:
        """Test that split property is cached."""
        t1 = Triplet(Split({1}, {2, 3}))
        profile = TripletProfile([t1])

        split1 = profile.split
        split2 = profile.split

        # Should return the same object (cached)
        assert split1 is split2


class TestTripletProfileMethods:
    """Tests for TripletProfile methods."""

    def test_get_weight(self) -> None:
        """Test get_weight method."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        profile = TripletProfile({t1: 0.8, t2: 0.2})

        assert profile.get_weight(t1) == 0.8
        assert profile.get_weight(t2) == 0.2

    def test_get_weight_not_found(self) -> None:
        """Test get_weight for triplet not in profile."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        t3 = Triplet(Split({3}, {1, 2}))
        profile = TripletProfile({t1: 0.8, t2: 0.2})

        assert profile.get_weight(t3) == 0.0

    def test_len(self) -> None:
        """Test __len__ method."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        t3 = Triplet(Split({3}, {1, 2}))
        profile = TripletProfile([t1, t2, t3])

        assert len(profile) == 3

    def test_iter(self) -> None:
        """Test __iter__ method."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        profile = TripletProfile([t1, t2])

        triplets_list = list(profile)
        assert len(triplets_list) == 2
        assert t1 in triplets_list
        assert t2 in triplets_list

    def test_contains(self) -> None:
        """Test __contains__ method."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        t3 = Triplet(Split({3}, {1, 2}))
        profile = TripletProfile([t1, t2])

        assert t1 in profile
        assert t2 in profile
        assert t3 not in profile

    def test_repr(self) -> None:
        """Test __repr__ method."""
        t1 = Triplet(Split({1}, {2, 3}))
        profile = TripletProfile([t1])

        repr_str = repr(profile)
        assert "TripletProfile" in repr_str
        assert "taxa" in repr_str or "1" in repr_str


class TestTripletProfileEdgeCases:
    """Tests for TripletProfile edge cases."""

    def test_duplicate_triplets_in_list(self) -> None:
        """Test that duplicate triplets in list raise an error."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        # Same triplet appears twice - should raise error
        with pytest.raises(ValueError, match="appears multiple times in the input"):
            TripletProfile([t1, t1, t2])

    def test_duplicate_triplets_in_dict(self) -> None:
        """Test that duplicate keys in dict use last value; single triplet must have weight 1.0."""
        t1 = Triplet(Split({1}, {2, 3}))
        # Same key twice: dict has only one entry (last wins). Weight must sum to 1.0.
        profile = TripletProfile({t1: 1.0})
        assert len(profile) == 1
        assert profile.get_weight(t1) == 1.0
        assert abs(sum(profile.triplets.values()) - 1.0) < 1e-9

    def test_provided_weights_stored_as_is(self) -> None:
        """Test that when weights are provided and sum to 1.0, they are stored without scaling."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        profile = TripletProfile({t1: 1.0 / 3.0, t2: 2.0 / 3.0})
        assert profile.get_weight(t1) == pytest.approx(1.0 / 3.0)
        assert profile.get_weight(t2) == pytest.approx(2.0 / 3.0)
        assert abs(sum(profile.triplets.values()) - 1.0) < 1e-9

    def test_weights_must_sum_to_one_error(self) -> None:
        """Test that provided weights that do not sum to 1.0 raise."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        with pytest.raises(ValueError, match="Weights must sum to 1.0"):
            TripletProfile({t1: 100.0, t2: 200.0})
        with pytest.raises(ValueError, match="Weights must sum to 1.0"):
            TripletProfile({t1: 0.0001, t2: 0.0002})
        with pytest.raises(ValueError, match="Weights must sum to 1.0"):
            TripletProfile({t1: 0.8})  # single triplet with weight != 1.0

    def test_many_triplets_equal_weight(self) -> None:
        """Test profile with many triplets gets equal weight 1/k each."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        t3 = Triplet(Split({3}, {1, 2}))
        star = Triplet({1, 2, 3})

        profile = TripletProfile([t1, t2, t3, star])

        assert len(profile) == 4
        assert abs(sum(profile.triplets.values()) - 1.0) < 1e-9
        assert all(profile.get_weight(t) == 0.25 for t in [t1, t2, t3, star])
        assert all(t in profile for t in [t1, t2, t3, star])


class TestTripletProfileIsResolved:
    """Tests for TripletProfile.is_resolved method."""

    def test_all_resolved_triplets(self) -> None:
        """Test is_resolved with all resolved triplets."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        profile = TripletProfile([t1, t2])

        assert profile.is_resolved() is True

    def test_single_resolved_triplet(self) -> None:
        """Test is_resolved with single resolved triplet."""
        t1 = Triplet(Split({1}, {2, 3}))
        profile = TripletProfile([t1])

        assert profile.is_resolved() is True

    def test_mixed_resolved_and_star(self) -> None:
        """Test is_resolved with mix of resolved and star triplets."""
        t1 = Triplet(Split({1}, {2, 3}))
        star = Triplet({1, 2, 3})
        profile = TripletProfile([t1, star])

        assert profile.is_resolved() is False

    def test_all_star_triplets(self) -> None:
        """Test is_resolved with all star triplets."""
        star1 = Triplet({1, 2, 3})
        # Can't have duplicate stars, so just use one
        profile = TripletProfile([star1])

        assert profile.is_resolved() is False

    def test_single_star_triplet(self) -> None:
        """Test is_resolved with single star triplet."""
        star = Triplet({1, 2, 3})
        profile = TripletProfile([star])

        assert profile.is_resolved() is False

    def test_three_resolved_triplets(self) -> None:
        """Test is_resolved with three resolved triplets."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        t3 = Triplet(Split({3}, {1, 2}))
        profile = TripletProfile([t1, t2, t3])

        assert profile.is_resolved() is True

    def test_multiple_resolved_one_star(self) -> None:
        """Test is_resolved with multiple resolved and one star."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        t3 = Triplet(Split({3}, {1, 2}))
        star = Triplet({1, 2, 3})
        profile = TripletProfile([t1, t2, t3, star])

        assert profile.is_resolved() is False

    def test_is_resolved_with_weights(self) -> None:
        """Test is_resolved works correctly with weighted triplets."""
        t1 = Triplet(Split({1}, {2, 3}))
        t2 = Triplet(Split({2}, {1, 3}))
        profile = TripletProfile({t1: 0.7, t2: 0.3})

        assert profile.is_resolved() is True

        star = Triplet({1, 2, 3})
        profile2 = TripletProfile({t1: 0.7, star: 0.3})

        assert profile2.is_resolved() is False
