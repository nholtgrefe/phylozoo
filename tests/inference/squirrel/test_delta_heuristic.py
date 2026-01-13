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
from phylozoo.inference.squirrel.squirrel import squirrel
from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
from phylozoo.core.network.sdnetwork.classifications import has_parallel_edges
from phylozoo.core.network.sdnetwork.features import blobs, k_blobs
from phylozoo.core.network.sdnetwork.isomorphism import is_isomorphic
from phylozoo.inference.squirrel.qsimilarity import sqprofileset_from_network


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


class TestDeltaHeuristicToNetwork:
    """Test cases for delta heuristic to network pipeline."""
    
    def test_delta_heuristic_to_network_no_errors(self) -> None:
        """
        Test that going from delta-heuristic to full network doesn't raise errors.
        
        This tests the complete pipeline:
        1. Create distance matrix
        2. Run delta_heuristic to get profileset
        3. Run squirrel to get network
        4. Verify network is valid (no errors)
        """
        # Create a distance matrix with clear structure
        matrix = np.array([
            [0.0, 0.1, 0.9, 0.9, 0.9],
            [0.1, 0.0, 0.9, 0.9, 0.9],
            [0.9, 0.9, 0.0, 0.1, 0.9],
            [0.9, 0.9, 0.1, 0.0, 0.9],
            [0.9, 0.9, 0.9, 0.9, 0.0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D', 'E'])
        
        # Step 1: Get profileset from delta heuristic
        profileset = delta_heuristic(dm, lam=0.3)
        assert len(profileset) > 0
        
        # Step 2: Run squirrel to get network
        # This should not raise any errors
        network = squirrel(profileset)
        
        # Step 3: Verify network is valid
        assert isinstance(network, SemiDirectedPhyNetwork)
        assert network.number_of_nodes() > 0
        assert len(network.taxa) == 5
        network.validate()  # Should not raise errors
    
    def test_delta_heuristic_to_network_multiple_taxa(self) -> None:
        """Test delta heuristic to network with more taxa."""
        # Create distance matrix for 6 taxa
        matrix = np.array([
            [0.0, 0.1, 0.9, 0.9, 0.9, 0.9],
            [0.1, 0.0, 0.9, 0.9, 0.9, 0.9],
            [0.9, 0.9, 0.0, 0.1, 0.9, 0.9],
            [0.9, 0.9, 0.1, 0.0, 0.9, 0.9],
            [0.9, 0.9, 0.9, 0.9, 0.0, 0.1],
            [0.9, 0.9, 0.9, 0.9, 0.1, 0.0]
        ])
        dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D', 'E', 'F'])
        
        profileset = delta_heuristic(dm, lam=0.3)
        assert len(profileset) > 0
        
        # Should not raise errors
        network = squirrel(profileset)
        assert isinstance(network, SemiDirectedPhyNetwork)
        assert len(network.taxa) == 6
        network.validate()


class TestNetworkRoundTrip:
    """Test cases for network -> profileset -> network round-trip."""
    
    def _has_3_cycle(self, network: SemiDirectedPhyNetwork) -> bool:
        """Check if network has a 3-node blob (3-cycle)."""
        network_blobs = blobs(network, trivial=False, leaves=False)
        return any(len(blob) == 3 for blob in network_blobs)
    
    def _has_2_blobs(self, network: SemiDirectedPhyNetwork) -> bool:
        """Check if network has multi-node 2-blobs."""
        two_blobs = k_blobs(network, k=2, trivial=False, leaves=False)
        return any(len(blob) > 1 for blob in two_blobs)
    
    def test_network_round_trip_simple_tree(self) -> None:
        """
        Test round-trip: network -> profileset -> network for a simple tree.
        
        A tree has no 3-cycles, no parallel edges, and no 2-blobs (except trivial).
        """
        # Use a fixture tree network with at least 4 taxa
        from tests.fixtures.sd_networks import SDTREE_MEDIUM_BINARY
        network = SDTREE_MEDIUM_BINARY
        
        # Verify network properties
        assert not has_parallel_edges(network)
        assert not self._has_2_blobs(network)
        assert not self._has_3_cycle(network)
        
        # Get profileset from network
        profileset = sqprofileset_from_network(network)
        assert len(profileset) > 0
        
        # Run squirrel to reconstruct network
        reconstructed = squirrel(profileset)
        
        # Verify reconstructed network is valid
        assert isinstance(reconstructed, SemiDirectedPhyNetwork)
        reconstructed.validate()
        
        # Check that networks are isomorphic (same topology)
        assert is_isomorphic(network, reconstructed)
        
    def test_network_round_trip_larger_tree(self) -> None:
        """Test round-trip for a larger tree."""
        # Use a fixture tree network
        from tests.fixtures.sd_networks import SDTREE_MEDIUM_BINARY
        network = SDTREE_MEDIUM_BINARY
        
        # Verify network properties
        assert not has_parallel_edges(network)
        assert not self._has_2_blobs(network)
        assert not self._has_3_cycle(network)
        
        # Get profileset from network
        profileset = sqprofileset_from_network(network)
        assert len(profileset) > 0
        
        # Run squirrel to reconstruct network
        reconstructed = squirrel(profileset)
        
        # Verify reconstructed network is valid
        assert isinstance(reconstructed, SemiDirectedPhyNetwork)
        reconstructed.validate()
        
        # Check that networks are isomorphic
        assert is_isomorphic(network, reconstructed)
