"""
Tests for the delta heuristic module.
"""

import pytest
import numpy as np

from phylozoo.core.distance import DistanceMatrix
from phylozoo.core.split.base import Split
from phylozoo.core.quartet.base import Quartet
from phylozoo.inference.squirrel.delta_heuristic import delta_heuristic
from phylozoo.inference.squirrel.sqprofileset import SqQuartetProfileSet


class TestDeltaHeuristic:
    """Test cases for delta heuristic algorithm."""
    
    def test_fewer_than_four_taxa_returns_empty(self) -> None:
        """Test that fewer than 4 taxa returns empty profileset."""
        # 3 taxa
        matrix = np.array([[0.0, 0.1, 0.2], [0.1, 0.0, 0.3], [0.2, 0.3, 0.0]])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C'])
        profileset = delta_heuristic(dm)
        assert len(profileset) == 0
    
    def test_four_taxa_creates_one_profile(self) -> None:
        """Test that 4 taxa creates exactly one profile."""
        # Create distance matrix with clear split: AB close, CD close, AB far from CD
        matrix = np.array([
            [0.0, 0.1, 0.9, 0.9],  # A: close to B, far from C, D
            [0.1, 0.0, 0.9, 0.9],  # B: close to A, far from C, D
            [0.9, 0.9, 0.0, 0.1],  # C: close to D, far from A, B
            [0.9, 0.9, 0.1, 0.0]   # D: close to C, far from A, B
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
        profileset = delta_heuristic(dm, lam=0.3)
        
        assert len(profileset) == 1
        profile = list(profileset.profiles.values())[0][0]
        assert len(profile.quartets) == 1  # Should be a split (delta < lam)
        
        # Check that it's the correct split: AB|CD
        quartet = next(iter(profile.quartets))
        assert quartet.split == Split({'A', 'B'}, {'C', 'D'})
    
    def test_delta_threshold_split_vs_cycle(self) -> None:
        """Test that delta threshold determines split vs cycle."""
        # Create distance matrix where delta will be high (cycle case)
        # All pairs equally distant (delta = 0, creates split)
        matrix = np.array([
            [0.0, 0.5, 0.5, 0.5],
            [0.5, 0.0, 0.5, 0.5],
            [0.5, 0.5, 0.0, 0.5],
            [0.5, 0.5, 0.5, 0.0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
        
        # With equal distances, delta = 0, so should create split
        profileset = delta_heuristic(dm, lam=0.3)
        assert len(profileset) == 1
        profile = list(profileset.profiles.values())[0][0]
        assert len(profile.quartets) == 1  # Split, not cycle
    
    def test_weighted_profiles(self) -> None:
        """Test that weighted profiles have correct weights."""
        matrix = np.array([
            [0.0, 0.1, 0.9, 0.9],
            [0.1, 0.0, 0.9, 0.9],
            [0.9, 0.9, 0.0, 0.1],
            [0.9, 0.9, 0.1, 0.0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
        
        # With weights
        profileset_weighted = delta_heuristic(dm, lam=0.3, weight=True)
        profile_weighted = list(profileset_weighted.profiles.values())[0][0]
        quartet = next(iter(profile_weighted.quartets))
        weight = profile_weighted.get_weight(quartet)
        assert weight is not None
        assert weight > 0
        
        # Without weights
        profileset_unweighted = delta_heuristic(dm, lam=0.3, weight=False)
        profile_unweighted = list(profileset_unweighted.profiles.values())[0][0]
        quartet_unweighted = next(iter(profile_unweighted.quartets))
        weight_unweighted = profile_unweighted.get_weight(quartet_unweighted)
        assert weight_unweighted == 1.0
    
    def test_reticulation_leaf_assignment(self) -> None:
        """Test that cycles get reticulation leaves assigned."""
        # Create distance matrix that will produce a cycle
        # Use distances that create high delta (>= lam)
        matrix = np.array([
            [0.0, 0.2, 0.3, 0.4],
            [0.2, 0.0, 0.4, 0.3],
            [0.3, 0.4, 0.0, 0.2],
            [0.4, 0.3, 0.2, 0.0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
        
        # Use low lambda to force cycles
        profileset = delta_heuristic(dm, lam=0.01, weight=False)
        
        # Check if any profiles are cycles (2 quartets)
        has_cycle = any(len(profile.quartets) == 2 for profile, _ in profileset.profiles.values())
        
        if has_cycle:
            # Find cycle profile
            cycle_profile = next(
                (profile for profile, _ in profileset.profiles.values() if len(profile.quartets) == 2),
                None
            )
            assert cycle_profile is not None
            # Should have reticulation leaf assigned
            assert cycle_profile.reticulation_leaf is not None
            assert cycle_profile.reticulation_leaf in cycle_profile.taxa
    
    def test_lambda_validation(self) -> None:
        """Test that invalid lambda values raise errors."""
        matrix = np.array([
            [0.0, 0.1, 0.9, 0.9],
            [0.1, 0.0, 0.9, 0.9],
            [0.9, 0.9, 0.0, 0.1],
            [0.9, 0.9, 0.1, 0.0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
        
        with pytest.raises(ValueError, match="Lambda must be in"):
            delta_heuristic(dm, lam=-0.1)
        
        with pytest.raises(ValueError, match="Lambda must be in"):
            delta_heuristic(dm, lam=1.5)
    
    def test_five_taxa_creates_five_profiles(self) -> None:
        """Test that 5 taxa creates C(5,4) = 5 profiles."""
        matrix = np.array([
            [0.0, 0.1, 0.2, 0.3, 0.4],
            [0.1, 0.0, 0.2, 0.3, 0.4],
            [0.2, 0.2, 0.0, 0.3, 0.4],
            [0.3, 0.3, 0.3, 0.0, 0.4],
            [0.4, 0.4, 0.4, 0.4, 0.0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D', 'E'])
        profileset = delta_heuristic(dm, lam=0.3)
        
        # C(5,4) = 5 combinations
        assert len(profileset) == 5
    
    def test_all_profiles_have_four_taxa(self) -> None:
        """Test that all profiles contain exactly 4 taxa."""
        matrix = np.array([
            [0.0, 0.1, 0.2, 0.3, 0.4],
            [0.1, 0.0, 0.2, 0.3, 0.4],
            [0.2, 0.2, 0.0, 0.3, 0.4],
            [0.3, 0.3, 0.3, 0.0, 0.4],
            [0.4, 0.4, 0.4, 0.4, 0.0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D', 'E'])
        profileset = delta_heuristic(dm, lam=0.3)
        
        for profile, _ in profileset.profiles.values():
            assert len(profile.taxa) == 4
    
    def test_profiles_are_valid_sqprofiles(self) -> None:
        """Test that all profiles are valid SqQuartetProfile objects."""
        matrix = np.array([
            [0.0, 0.1, 0.9, 0.9],
            [0.1, 0.0, 0.9, 0.9],
            [0.9, 0.9, 0.0, 0.1],
            [0.9, 0.9, 0.1, 0.0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
        profileset = delta_heuristic(dm, lam=0.3)
        
        for profile, _ in profileset.profiles.values():
            # Should have 1 or 2 quartets
            assert len(profile.quartets) in (1, 2)
            # All quartets should be resolved
            for quartet in profile.quartets:
                assert quartet.is_resolved()

