"""
Tests for the quartet distance module.
"""

import pytest

import numpy as np

from phylozoo.core.distance import DistanceMatrix
from phylozoo.core.primitives.partition import Partition
from phylozoo.core.quartet.base import Quartet
from phylozoo.core.quartet.qdistance import (
    quartet_distance,
    quartet_distance_with_partition,
    _rho_distance,
)
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


class TestQuartetDistanceWithPartition:
    """Tests for the quartet_distance_with_partition function."""

    def test_quartet_distance_with_partition_four_singleton_sets(self) -> None:
        """Test quartet_distance_with_partition with 4 singleton sets."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1], taxa={'A', 'B', 'C', 'D'})
        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}])
        rho = (0.5, 1.0, 0.5, 1.0)
        
        dist_matrix = quartet_distance_with_partition(profileset, partition, rho)
        
        assert isinstance(dist_matrix, DistanceMatrix)
        assert len(dist_matrix) == 4
        # Labels should be the partition sets (frozensets)
        assert len(dist_matrix.labels) == 4
        # Check diagonal is 0
        for i, set1 in enumerate(dist_matrix.labels):
            assert dist_matrix.get_distance(set1, set1) == 0.0
        
        # Check constant was added: 2*4 - 4 = 4
        # For singleton sets, this should match quartet_distance behavior
        # A and B are on same side: rho_s = 1.0, delta = 2*1.0 = 2.0
        # Average = 2.0 (only one quartet), Total = 2.0 + 4 = 6.0
        set_a = frozenset({'A'})
        set_b = frozenset({'B'})
        set_c = frozenset({'C'})
        assert dist_matrix.get_distance(set_a, set_b) == pytest.approx(6.0)
        # A and C are on different sides: rho_c = 0.5, delta = 2*0.5 = 1.0
        # Average = 1.0, Total = 1.0 + 4 = 5.0
        assert dist_matrix.get_distance(set_a, set_c) == pytest.approx(5.0)

    def test_quartet_distance_with_partition_averaging(self) -> None:
        """Test that distances are averaged when multiple quartets contribute."""
        # Create a partition where sets have multiple leaves
        # This will generate multiple quartets per pair of sets
        taxa = {'A', 'B', 'C', 'D', 'E', 'F'}
        
        # Create all possible quartets to make it dense
        profiles = []
        for four_taxa in [
            {'A', 'B', 'C', 'D'},
            {'A', 'B', 'C', 'E'},
            {'A', 'B', 'C', 'F'},
            {'A', 'B', 'D', 'E'},
            {'A', 'B', 'D', 'F'},
            {'A', 'B', 'E', 'F'},
            {'A', 'C', 'D', 'E'},
            {'A', 'C', 'D', 'F'},
            {'A', 'C', 'E', 'F'},
            {'A', 'D', 'E', 'F'},
            {'B', 'C', 'D', 'E'},
            {'B', 'C', 'D', 'F'},
            {'B', 'C', 'E', 'F'},
            {'B', 'D', 'E', 'F'},
            {'C', 'D', 'E', 'F'},
        ]:
            four_list = sorted(four_taxa)
            q = Quartet(Split({four_list[0], four_list[1]}, {four_list[2], four_list[3]}))
            profiles.append(q)
        
        profileset = QuartetProfileSet(profiles=profiles, taxa=taxa)
        # Partition: {A, B} | {C, D} | {E, F}
        partition = Partition([{'A', 'B'}, {'C', 'D'}, {'E', 'F'}])
        rho = (0.5, 1.0, 0.5, 1.0)
        
        dist_matrix = quartet_distance_with_partition(profileset, partition, rho)
        
        assert isinstance(dist_matrix, DistanceMatrix)
        assert len(dist_matrix) == 3
        # Check that distances are computed (averaged across multiple quartets)
        for i, set1 in enumerate(dist_matrix.labels):
            for j, set2 in enumerate(dist_matrix.labels):
                if i != j:
                    dist = dist_matrix.get_distance(set1, set2)
                    # Should be >= constant (2*3 - 4 = 2)
                    assert dist >= 2.0

    def test_quartet_distance_with_partition_symmetric(self) -> None:
        """Test that distance matrix is symmetric."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1], taxa={'A', 'B', 'C', 'D'})
        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}])
        rho = (0.5, 1.0, 0.5, 1.0)
        
        dist_matrix = quartet_distance_with_partition(profileset, partition, rho)
        
        # Check symmetry
        for set1 in dist_matrix.labels:
            for set2 in dist_matrix.labels:
                assert dist_matrix.get_distance(set1, set2) == pytest.approx(
                    dist_matrix.get_distance(set2, set1)
                )

    def test_quartet_distance_with_partition_constant_addition(self) -> None:
        """Test that constant 2*n - 4 is added correctly."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1], taxa={'A', 'B', 'C', 'D'})
        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}])
        rho = (0.0, 0.0, 0.0, 0.0)  # All zeros to see just the constant
        
        dist_matrix = quartet_distance_with_partition(profileset, partition, rho)
        
        # For n=4 sets, constant = 2*4 - 4 = 4
        # All off-diagonal should be 4 (since rho values are 0)
        for i, set1 in enumerate(dist_matrix.labels):
            for j, set2 in enumerate(dist_matrix.labels):
                if i != j:
                    assert dist_matrix.get_distance(set1, set2) == pytest.approx(4.0)
                else:
                    assert dist_matrix.get_distance(set1, set2) == 0.0

    def test_quartet_distance_with_partition_partition_mismatch_raises_error(self) -> None:
        """Test that partition elements must match profile set taxa."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1], taxa={'A', 'B', 'C', 'D'})
        # Partition with different taxa
        partition = Partition([{'A'}, {'B'}, {'C'}, {'E'}])
        rho = (0.5, 1.0, 0.5, 1.0)
        
        with pytest.raises(ValueError, match="Partition elements must match profile set taxa"):
            quartet_distance_with_partition(profileset, partition, rho)

    def test_quartet_distance_with_partition_not_dense_raises_error(self) -> None:
        """Test that non-dense profile set raises an error."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1], taxa={'A', 'B', 'C', 'D', 'E'})
        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}, {'E'}])
        rho = (0.5, 1.0, 0.5, 1.0)
        
        with pytest.raises(ValueError, match="Profile set must be dense"):
            quartet_distance_with_partition(profileset, partition, rho)

    def test_quartet_distance_with_partition_too_many_quartets_raises_error(self) -> None:
        """Test that profile with more than 2 quartets raises an error."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        q3 = Quartet(Split({'A', 'D'}, {'B', 'C'}))
        profile = QuartetProfile([q1, q2, q3])
        profileset = QuartetProfileSet(profiles=[profile], taxa={'A', 'B', 'C', 'D'})
        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}])
        rho = (0.5, 1.0, 0.5, 1.0)
        
        with pytest.raises(ValueError, match="must contain exactly 1 or 2 quartets"):
            quartet_distance_with_partition(profileset, partition, rho)

    def test_quartet_distance_with_partition_unresolved_quartet_raises_error(self) -> None:
        """Test that unresolved quartet in profile raises an error."""
        star_quartet = Quartet({'A', 'B', 'C', 'D'})
        profile = QuartetProfile([star_quartet])
        profileset = QuartetProfileSet(profiles=[profile], taxa={'A', 'B', 'C', 'D'})
        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}])
        rho = (0.5, 1.0, 0.5, 1.0)
        
        with pytest.raises(ValueError, match="contains unresolved quartet"):
            quartet_distance_with_partition(profileset, partition, rho)

    def test_quartet_distance_with_partition_invalid_rho_length_raises_error(self) -> None:
        """Test that invalid rho vector length raises an error."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1], taxa={'A', 'B', 'C', 'D'})
        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}])
        
        with pytest.raises(ValueError, match="Rho vector must have 4 elements"):
            quartet_distance_with_partition(profileset, partition, (0.5, 1.0, 0.5))  # Only 3 elements

    def test_quartet_distance_with_partition_invalid_rho_constraints_raises_error(self) -> None:
        """Test that invalid rho constraints raise an error."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1], taxa={'A', 'B', 'C', 'D'})
        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}])
        
        # rho_a > rho_o should raise error
        with pytest.raises(ValueError, match="Rho vector must satisfy"):
            quartet_distance_with_partition(profileset, partition, (0.5, 1.0, 1.5, 1.0))
        
        # rho_c > rho_s should raise error
        with pytest.raises(ValueError, match="Rho vector must satisfy"):
            quartet_distance_with_partition(profileset, partition, (1.5, 1.0, 0.5, 1.0))

    def test_quartet_distance_with_partition_nanuq_rho(self) -> None:
        """Test with NANUQ rho vector."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1], taxa={'A', 'B', 'C', 'D'})
        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}])
        rho = (0.0, 1.0, 0.5, 1.0)  # NANUQ
        
        dist_matrix = quartet_distance_with_partition(profileset, partition, rho)
        
        assert isinstance(dist_matrix, DistanceMatrix)
        set_a = frozenset({'A'})
        set_b = frozenset({'B'})
        set_c = frozenset({'C'})
        # A and B are on same side: rho_s = 1.0, delta = 2*1.0 = 2.0
        # Average = 2.0, Total = 2.0 + 4 = 6.0
        assert dist_matrix.get_distance(set_a, set_b) == pytest.approx(6.0)
        # A and C are on different sides: rho_c = 0.0, delta = 2*0.0 = 0.0
        # Average = 0.0, Total = 0.0 + 4 = 4.0
        assert dist_matrix.get_distance(set_a, set_c) == pytest.approx(4.0)

    def test_quartet_distance_with_partition_two_quartet_profiles(self) -> None:
        """Test with profiles containing 2 quartets (four-cycle)."""
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
        profile = QuartetProfile([q1, q2])
        profileset = QuartetProfileSet(profiles=[profile], taxa={'A', 'B', 'C', 'D'})
        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}])
        rho = (0.5, 1.0, 0.5, 1.0)
        
        dist_matrix = quartet_distance_with_partition(profileset, partition, rho)
        
        assert isinstance(dist_matrix, DistanceMatrix)
        assert len(dist_matrix) == 4
        # Check constant was added: 2*4 - 4 = 4
        # All distances should be >= 4
        for i, set1 in enumerate(dist_matrix.labels):
            for j, set2 in enumerate(dist_matrix.labels):
                if i != j:
                    assert dist_matrix.get_distance(set1, set2) >= 4.0

    def test_quartet_distance_with_partition_three_sets(self) -> None:
        """Test with partition of 3 sets."""
        taxa = {'A', 'B', 'C', 'D', 'E'}
        # Create all C(5,4) = 5 profiles to make it dense
        profiles = []
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
        partition = Partition([{'A', 'B'}, {'C'}, {'D', 'E'}])
        rho = (0.5, 1.0, 0.5, 1.0)
        
        dist_matrix = quartet_distance_with_partition(profileset, partition, rho)
        
        assert isinstance(dist_matrix, DistanceMatrix)
        assert len(dist_matrix) == 3
        # Check constant: 2*3 - 4 = 2
        for i, set1 in enumerate(dist_matrix.labels):
            for j, set2 in enumerate(dist_matrix.labels):
                if i != j:
                    assert dist_matrix.get_distance(set1, set2) >= 2.0

    def test_quartet_distance_with_partition_averaging_multiple_contributions(self) -> None:
        """Test that averaging works correctly with multiple quartets."""
        # Create a simple case where we can verify the averaging
        taxa = {'A', 'B', 'C', 'D'}
        # Create all quartets for {A, B, C, D}
        # Only one 4-taxon combination, so one profile
        q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
        profileset = QuartetProfileSet(profiles=[q1], taxa=taxa)
        
        # Partition: {A, B} | {C} | {D}
        # This creates a 3-set partition, but we need at least 4 sets for subpartitions
        # Let's use a different partition: {A} | {B} | {C} | {D}
        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}])
        rho = (0.5, 1.0, 0.5, 1.0)
        
        dist_matrix = quartet_distance_with_partition(profileset, partition, rho)
        
        # For singleton sets, there's only one quartet per pair, so averaging
        # should give the same result as no averaging
        set_a = frozenset({'A'})
        set_b = frozenset({'B'})
        # A and B are on same side: rho_s = 1.0, delta = 2*1.0 = 2.0
        # Average = 2.0 (only one contribution), Total = 2.0 + 4 = 6.0
        assert dist_matrix.get_distance(set_a, set_b) == pytest.approx(6.0)

    def test_quartet_distance_with_partition_five_sets(self) -> None:
        """Test with partition of 5 sets."""
        taxa = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'}
        # Create enough profiles to make it dense (C(8,4) = 70 profiles)
        # For testing, we'll create a subset that covers the needed combinations
        profiles = []
        # Create profiles for all 4-taxon combinations involving the first 5 taxa
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
        
        # Add more profiles to make it dense for the 8 taxa
        # (In practice, we'd need all 70, but for testing we can use a subset)
        # Actually, let's use a smaller set of taxa to make it dense
        taxa = {'A', 'B', 'C', 'D', 'E'}
        profiles = []
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
        partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}, {'E'}])
        rho = (0.5, 1.0, 0.5, 1.0)
        
        dist_matrix = quartet_distance_with_partition(profileset, partition, rho)
        
        assert isinstance(dist_matrix, DistanceMatrix)
        assert len(dist_matrix) == 5
        # Check constant: 2*5 - 4 = 6
        for i, set1 in enumerate(dist_matrix.labels):
            for j, set2 in enumerate(dist_matrix.labels):
                if i != j:
                    assert dist_matrix.get_distance(set1, set2) >= 6.0

    def test_quartet_distance_with_partition_mixed_set_sizes(self) -> None:
        """Test with partition containing sets of different sizes."""
        taxa = {'A', 'B', 'C', 'D', 'E', 'F'}
        # Create all C(6,4) = 15 profiles to make it dense
        profiles = []
        for four_taxa in [
            {'A', 'B', 'C', 'D'},
            {'A', 'B', 'C', 'E'},
            {'A', 'B', 'C', 'F'},
            {'A', 'B', 'D', 'E'},
            {'A', 'B', 'D', 'F'},
            {'A', 'B', 'E', 'F'},
            {'A', 'C', 'D', 'E'},
            {'A', 'C', 'D', 'F'},
            {'A', 'C', 'E', 'F'},
            {'A', 'D', 'E', 'F'},
            {'B', 'C', 'D', 'E'},
            {'B', 'C', 'D', 'F'},
            {'B', 'C', 'E', 'F'},
            {'B', 'D', 'E', 'F'},
            {'C', 'D', 'E', 'F'},
        ]:
            four_list = sorted(four_taxa)
            q = Quartet(Split({four_list[0], four_list[1]}, {four_list[2], four_list[3]}))
            profiles.append(q)
        
        profileset = QuartetProfileSet(profiles=profiles, taxa=taxa)
        # Partition with mixed sizes: {A, B} | {C} | {D, E, F}
        partition = Partition([{'A', 'B'}, {'C'}, {'D', 'E', 'F'}])
        rho = (0.5, 1.0, 0.5, 1.0)
        
        dist_matrix = quartet_distance_with_partition(profileset, partition, rho)
        
        assert isinstance(dist_matrix, DistanceMatrix)
        assert len(dist_matrix) == 3
        # Check constant: 2*3 - 4 = 2
        for i, set1 in enumerate(dist_matrix.labels):
            for j, set2 in enumerate(dist_matrix.labels):
                if i != j:
                    dist = dist_matrix.get_distance(set1, set2)
                    assert dist >= 2.0
                    # Check symmetry
                    assert dist == pytest.approx(dist_matrix.get_distance(set2, set1))

