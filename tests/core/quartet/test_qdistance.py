"""
Tests for the quartet distance module.
"""

import pytest

import numpy as np

from phylozoo.core.distance import DistanceMatrix
from phylozoo.core.quartet.base import Quartet
from phylozoo.core.quartet.qdistance import quartet_distance, _rho_distance
from phylozoo.core.quartet.qprofile import QuartetProfile
from phylozoo.core.quartet.qprofileset import QuartetProfileSet
from phylozoo.core.split.base import Split


class TestRhoDistance:
    """Tests for the _rho_distance helper function."""

    def test_rho_distance_single_quartet_same_side(self) -> None:
        """Test rho_distance for single quartet with leaves on same side."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profile = QuartetProfile([q1])
        rho = (0.5, 1.0, 0.5, 1.0)
        dist = _rho_distance(profile, 'A', 'B', rho)
        assert dist == 1.0  # rho_s

    def test_rho_distance_single_quartet_different_sides(self) -> None:
        """Test rho_distance for single quartet with leaves on different sides."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profile = QuartetProfile([q1])
        rho = (0.5, 1.0, 0.5, 1.0)
        dist = _rho_distance(profile, 'A', 'C', rho)
        assert dist == 0.5  # rho_c

    def test_rho_distance_two_quartets_adjacent(self) -> None:
        """Test rho_distance for two quartets with adjacent leaves."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profile = QuartetProfile([q1, q2])
        rho = (0.5, 1.0, 0.5, 1.0)
        # Check circular ordering - A and B should be adjacent
        if profile.circular_orderings:
            ordering = next(iter(profile.circular_orderings))
            if ordering.are_neighbors('A', 'B'):
                dist = _rho_distance(profile, 'A', 'B', rho)
                assert dist == 0.5  # rho_a

    def test_rho_distance_two_quartets_opposite(self) -> None:
        """Test rho_distance for two quartets with opposite leaves."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profile = QuartetProfile([q1, q2])
        rho = (0.5, 1.0, 0.5, 1.0)
        # Check circular ordering - A and D should be opposite
        if profile.circular_orderings:
            ordering = next(iter(profile.circular_orderings))
            if not ordering.are_neighbors('A', 'D'):
                dist = _rho_distance(profile, 'A', 'D', rho)
                assert dist == 1.0  # rho_o

    def test_rho_distance_unresolved_quartet_raises_error(self) -> None:
        """Test that unresolved quartet raises an error."""
        # Create a star tree quartet
        star_quartet = Quartet({'A', 'B', 'C', 'D'})
        profile = QuartetProfile([star_quartet])
        rho = (0.5, 1.0, 0.5, 1.0)
        with pytest.raises(ValueError, match="Quartet must be resolved"):
            _rho_distance(profile, 'A', 'B', rho)

    def test_rho_distance_wrong_number_of_quartets_raises_error(self) -> None:
        """Test that wrong number of quartets raises an error."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        q3 = Quartet(Split({'A', 'D'}, {'B', 'C'}))
        profile = QuartetProfile([q1, q2, q3])
        rho = (0.5, 1.0, 0.5, 1.0)
        with pytest.raises(ValueError, match="Profile must have 1 or 2 quartets"):
            _rho_distance(profile, 'A', 'B', rho)


class TestQuartetDistance:
    """Tests for the quartet_distance function."""

    def test_quartet_distance_four_taxa_dense(self) -> None:
        """Test quartet_distance with 4 taxa (dense profile set)."""
        # Create all possible quartets for 4 taxa: {A, B, C, D}
        # Only one 4-taxon combination, so one profile
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1], taxa={'A', 'B', 'C', 'D'})
        rho = (0.5, 1.0, 0.5, 1.0)
        
        dist_matrix = quartet_distance(profileset, rho)
        
        assert isinstance(dist_matrix, DistanceMatrix)
        assert len(dist_matrix) == 4
        assert dist_matrix.labels == ('A', 'B', 'C', 'D')
        # Check diagonal is 0
        for taxon in dist_matrix.labels:
            assert dist_matrix.get_distance(taxon, taxon) == 0.0
        # Check constant was added: 2*4 - 4 = 4
        # A and B are on same side, so rho_s = 1.0, delta = 2*1.0 = 2.0
        # Total = 2.0 + 4 = 6.0
        assert dist_matrix.get_distance('A', 'B') == pytest.approx(6.0)
        # A and C are on different sides, so rho_c = 0.5, delta = 2*0.5 = 1.0
        # Total = 1.0 + 4 = 5.0
        assert dist_matrix.get_distance('A', 'C') == pytest.approx(5.0)

    def test_quartet_distance_five_taxa_dense(self) -> None:
        """Test quartet_distance with 5 taxa (dense profile set)."""
        # Create all C(5,4) = 5 profiles
        taxa = {'A', 'B', 'C', 'D', 'E'}
        profiles = []
        # Create quartets for each 4-taxon combination
        for four_taxa in [
            {'A', 'B', 'C', 'D'},
            {'A', 'B', 'C', 'E'},
            {'A', 'B', 'D', 'E'},
            {'A', 'C', 'D', 'E'},
            {'B', 'C', 'D', 'E'},
        ]:
            four_list = sorted(four_taxa)
            q = Quartet(Split({four_list[0], four_list[1]}, {four_list[2], four_list[3]}))
            profiles.append(q)
        
        profileset = QuartetProfileSet(profiles=profiles, taxa=taxa)
        rho = (0.5, 1.0, 0.5, 1.0)
        
        dist_matrix = quartet_distance(profileset, rho)
        
        assert isinstance(dist_matrix, DistanceMatrix)
        assert len(dist_matrix) == 5
        # Check constant was added: 2*5 - 4 = 6
        # Diagonal should be 0
        for taxon in dist_matrix.labels:
            assert dist_matrix.get_distance(taxon, taxon) == 0.0
        # All off-diagonal entries should be >= 6 (constant)
        for i, taxon1 in enumerate(dist_matrix.labels):
            for j, taxon2 in enumerate(dist_matrix.labels):
                if i != j:
                    assert dist_matrix.get_distance(taxon1, taxon2) >= 6.0

    def test_quartet_distance_with_two_quartet_profiles(self) -> None:
        """Test quartet_distance with profiles containing 2 quartets."""
        # Create a profile with 2 quartets (four-cycle)
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profile = QuartetProfile([q1, q2])
        profileset = QuartetProfileSet(profiles=[profile], taxa={'A', 'B', 'C', 'D'})
        rho = (0.5, 1.0, 0.5, 1.0)
        
        dist_matrix = quartet_distance(profileset, rho)
        
        assert isinstance(dist_matrix, DistanceMatrix)
        assert len(dist_matrix) == 4
        # Check constant was added: 2*4 - 4 = 4
        # All distances should be >= 4
        for i, taxon1 in enumerate(dist_matrix.labels):
            for j, taxon2 in enumerate(dist_matrix.labels):
                if i != j:
                    assert dist_matrix.get_distance(taxon1, taxon2) >= 4.0

    def test_quartet_distance_not_dense_raises_error(self) -> None:
        """Test that non-dense profile set raises an error."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1], taxa={'A', 'B', 'C', 'D', 'E'})
        rho = (0.5, 1.0, 0.5, 1.0)
        
        with pytest.raises(ValueError, match="Profile set must be dense"):
            quartet_distance(profileset, rho)

    def test_quartet_distance_too_many_quartets_raises_error(self) -> None:
        """Test that profile with more than 2 quartets raises an error."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        q3 = Quartet(Split({'A', 'D'}, {'B', 'C'}))
        profile = QuartetProfile([q1, q2, q3])
        profileset = QuartetProfileSet(profiles=[profile], taxa={'A', 'B', 'C', 'D'})
        rho = (0.5, 1.0, 0.5, 1.0)
        
        with pytest.raises(ValueError, match="must contain exactly 1 or 2 quartets"):
            quartet_distance(profileset, rho)

    def test_quartet_distance_unresolved_quartet_raises_error(self) -> None:
        """Test that unresolved quartet in profile raises an error."""
        # Create a star tree quartet
        star_quartet = Quartet({'A', 'B', 'C', 'D'})
        profile = QuartetProfile([star_quartet])
        profileset = QuartetProfileSet(profiles=[profile], taxa={'A', 'B', 'C', 'D'})
        rho = (0.5, 1.0, 0.5, 1.0)
        
        with pytest.raises(ValueError, match="contains unresolved quartet"):
            quartet_distance(profileset, rho)

    def test_quartet_distance_invalid_rho_length_raises_error(self) -> None:
        """Test that invalid rho vector length raises an error."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1], taxa={'A', 'B', 'C', 'D'})
        
        with pytest.raises(ValueError, match="Rho vector must have 4 elements"):
            quartet_distance(profileset, (0.5, 1.0, 0.5))  # Only 3 elements

    def test_quartet_distance_invalid_rho_constraints_raises_error(self) -> None:
        """Test that invalid rho constraints raise an error."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1], taxa={'A', 'B', 'C', 'D'})
        
        # rho_a > rho_o should raise error
        with pytest.raises(ValueError, match="Rho vector must satisfy"):
            quartet_distance(profileset, (0.5, 1.0, 1.5, 1.0))
        
        # rho_c > rho_s should raise error
        with pytest.raises(ValueError, match="Rho vector must satisfy"):
            quartet_distance(profileset, (1.5, 1.0, 0.5, 1.0))

    def test_quartet_distance_symmetric(self) -> None:
        """Test that distance matrix is symmetric."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1], taxa={'A', 'B', 'C', 'D'})
        rho = (0.5, 1.0, 0.5, 1.0)
        
        dist_matrix = quartet_distance(profileset, rho)
        
        # Check symmetry
        for taxon1 in dist_matrix.labels:
            for taxon2 in dist_matrix.labels:
                assert dist_matrix.get_distance(taxon1, taxon2) == pytest.approx(
                    dist_matrix.get_distance(taxon2, taxon1)
                )

    def test_quartet_distance_constant_addition(self) -> None:
        """Test that constant 2*n - 4 is added correctly."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1], taxa={'A', 'B', 'C', 'D'})
        rho = (0.0, 0.0, 0.0, 0.0)  # All zeros to see just the constant
        
        dist_matrix = quartet_distance(profileset, rho)
        
        # For n=4, constant = 2*4 - 4 = 4
        # All off-diagonal should be 4 (since rho values are 0)
        for i, taxon1 in enumerate(dist_matrix.labels):
            for j, taxon2 in enumerate(dist_matrix.labels):
                if i != j:
                    assert dist_matrix.get_distance(taxon1, taxon2) == pytest.approx(4.0)
                else:
                    assert dist_matrix.get_distance(taxon1, taxon2) == 0.0

    def test_quartet_distance_nanuq_rho(self) -> None:
        """Test with NANUQ rho vector."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1], taxa={'A', 'B', 'C', 'D'})
        rho = (0.0, 1.0, 0.5, 1.0)  # NANUQ
        
        dist_matrix = quartet_distance(profileset, rho)
        
        assert isinstance(dist_matrix, DistanceMatrix)
        # A and B are on same side: rho_s = 1.0, delta = 2*1.0 = 2.0
        # Total = 2.0 + 4 = 6.0
        assert dist_matrix.get_distance('A', 'B') == pytest.approx(6.0)
        # A and C are on different sides: rho_c = 0.0, delta = 2*0.0 = 0.0
        # Total = 0.0 + 4 = 4.0
        assert dist_matrix.get_distance('A', 'C') == pytest.approx(4.0)

    def test_quartet_distance_mixed_profiles(self) -> None:
        """Test with mix of 1-quartet and 2-quartet profiles."""
        # Create profiles: some with 1 quartet, some with 2 quartets
        taxa = {'A', 'B', 'C', 'D', 'E'}
        profiles = []
        
        # Profile 1: {A, B, C, D} with 1 quartet
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profiles.append(QuartetProfile([q1]))
        
        # Profile 2: {A, B, C, E} with 2 quartets (four-cycle)
        q2 = Quartet(Split({'A', 'B'}, {'C', 'E'}))
        q3 = Quartet(Split({'A', 'C'}, {'B', 'E'}))
        profile2 = QuartetProfile([q2, q3])
        profiles.append(profile2)
        
        # Add remaining profiles to make it dense
        # {A, B, D, E}
        q4 = Quartet(Split({'A', 'B'}, {'D', 'E'}))
        profiles.append(QuartetProfile([q4]))
        # {A, C, D, E}
        q5 = Quartet(Split({'A', 'C'}, {'D', 'E'}))
        profiles.append(QuartetProfile([q5]))
        # {B, C, D, E}
        q6 = Quartet(Split({'B', 'C'}, {'D', 'E'}))
        profiles.append(QuartetProfile([q6]))
        
        profileset = QuartetProfileSet(profiles=profiles, taxa=taxa)
        rho = (0.5, 1.0, 0.5, 1.0)
        
        dist_matrix = quartet_distance(profileset, rho)
        
        assert isinstance(dist_matrix, DistanceMatrix)
        assert len(dist_matrix) == 5
        # Check constant: 2*5 - 4 = 6
        for i, taxon1 in enumerate(dist_matrix.labels):
            for j, taxon2 in enumerate(dist_matrix.labels):
                if i != j:
                    assert dist_matrix.get_distance(taxon1, taxon2) >= 6.0

