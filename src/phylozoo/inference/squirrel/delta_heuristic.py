"""
Delta heuristic module.

This module provides the delta heuristic algorithm for inferring quartet profiles
from distance matrices. The delta heuristic uses pairwise distances to determine
whether quartets should be resolved as splits (single quartet) or cycles (two quartets).
"""

from __future__ import annotations

import itertools

from ...core.distance import DistanceMatrix
from ...core.quartet.base import Quartet
from ...core.split.base import Split
from ...core.primitives.circular_ordering import CircularOrdering
from ...utils.exceptions import PhyloZooValueError
from .sqprofile import SqQuartetProfile
from .sqprofileset import SqQuartetProfileSet


def delta_heuristic(
    distance_matrix: DistanceMatrix,
    lam: float = 0.3,
    weight: bool = True,
) -> SqQuartetProfileSet:
    """
    Infer quartet profiles from a distance matrix using the delta heuristic.
    
    The delta heuristic algorithm:
    1. For each combination of 4 taxa, computes distances for the 3 possible splits
    2. Sorts splits by their total distance (sum of pairwise distances within each side)
    3. Computes a delta value: (d2 - d1) / (d2 - d0) where d0 <= d1 <= d2
    4. If delta < lambda: creates a single quartet (split) profile
    5. If delta >= lambda: creates a cycle (two quartets) profile with circular ordering
    6. Assigns reticulation leaves to cycles based on delta sum per taxon
    
    Parameters
    ----------
    distance_matrix : DistanceMatrix
        The distance matrix containing pairwise distances between taxa.
    lam : float, optional
        Lambda threshold for determining split vs cycle. Values < lam create splits,
        values >= lam create cycles. By default 0.3.
    weight : bool, optional
        If True, quartets are weighted based on delta value. By default True.
    
    Returns
    -------
    SqQuartetProfileSet
        A set of squirrel quartet profiles inferred from the distance matrix.
    
    Raises
    ------
    PhyloZooValueError
        If lambda is not in [0, 1].
    
    Examples
    --------
    >>> import numpy as np
    >>> from phylozoo.core.distance import DistanceMatrix
    >>> from phylozoo.inference.squirrel.delta_heuristic import delta_heuristic
    >>> 
    >>> # Create a simple distance matrix
    >>> matrix = np.array([
    ...     [0.0, 0.1, 0.9, 0.9],
    ...     [0.1, 0.0, 0.9, 0.9],
    ...     [0.9, 0.9, 0.0, 0.1],
    ...     [0.9, 0.9, 0.1, 0.0]
    ... ])
    >>> dm = DistanceMatrix(matrix, labels=['A', 'B', 'C', 'D'])
    >>> profileset = delta_heuristic(dm, lam=0.3)
    >>> len(profileset) == 1  # Only one quartet for 4 taxa
    True
    
    Notes
    -----
    The algorithm is optimized for speed using numpy vectorization where possible.
    Delta values are computed as: delta = (d2 - d1) / (d2 - d0) where:
    - d0 = smallest split distance
    - d1 = middle split distance  
    - d2 = largest split distance
    
    If all distances are equal, delta = 0 (creates a split).
    """
    if lam < 0 or lam > 1:
        raise PhyloZooValueError(f"Lambda must be in [0, 1], got {lam}")
    
    taxa_list = list(distance_matrix.labels)
    n = len(taxa_list)
    
    if n < 4:
        # Not enough taxa for quartets
        return SqQuartetProfileSet()
    
    # Build taxon to index mapping for fast lookups
    taxon_to_idx = {taxon: idx for idx, taxon in enumerate(taxa_list)}
    
    # Get distance matrix as numpy array for fast access
    dist_array = distance_matrix.np_array
    
    # Track delta sums per taxon for reticulation assignment
    delta_sum: dict[str, float] = {taxon: 0.0 for taxon in taxa_list}
    
    # Collect profiles
    profiles: list[SqQuartetProfile] = []
    
    # Iterate through all combinations of 4 taxa
    # Use numpy for faster distance lookups
    for i, j, k, l in itertools.combinations(range(n), 4):
        a, b, c, d = taxa_list[i], taxa_list[j], taxa_list[k], taxa_list[l]
        
        # Compute distances for the 3 possible splits (vectorized)
        # Split 1: {a,b} | {c,d}
        split1_dist = dist_array[i, j] + dist_array[k, l]
        
        # Split 2: {a,c} | {b,d}
        split2_dist = dist_array[i, k] + dist_array[j, l]
        
        # Split 3: {a,d} | {b,c}
        split3_dist = dist_array[i, l] + dist_array[j, k]
        
        # Store distances and split indices (avoid creating Split objects until needed)
        dists = [split1_dist, split2_dist, split3_dist]
        split_configs = [
            ({a, b}, {c, d}),  # Split 1
            ({a, c}, {b, d}),  # Split 2
            ({a, d}, {b, c}),  # Split 3
        ]
        
        # Sort by distance (ascending) - use indices to avoid recreating splits
        sorted_indices = sorted(range(3), key=lambda idx: dists[idx])
        
        # Extract distances
        d0, d1, d2 = dists[sorted_indices[0]], dists[sorted_indices[1]], dists[sorted_indices[2]]
        
        # Compute delta value
        if d0 == d1 == d2:
            delta = 0.0
        else:
            numerator = d2 - d1
            denominator = d2 - d0
            if denominator == 0:
                delta = 0.0
            else:
                delta = numerator / denominator
        
        # Update delta sums
        delta_sum[a] += delta
        delta_sum[b] += delta
        delta_sum[c] += delta
        delta_sum[d] += delta
        
        # Determine quartet type based on delta
        if delta < lam:
            # Create single quartet (split) profile
            best_idx = sorted_indices[0]
            best_set1, best_set2 = split_configs[best_idx]
            best_split = Split(best_set1, best_set2)
            quartet = Quartet(best_split)
            
            if weight:
                w = abs(lam - delta) / lam if lam > 0 else 1.0
                profile = SqQuartetProfile({quartet: w})
            else:
                profile = SqQuartetProfile([quartet])
        else:
            # Create cycle (two quartets) profile
            # The "wrong" split is the one with largest distance
            wrong_idx = sorted_indices[2]
            wrong_set1, wrong_set2 = split_configs[wrong_idx]
            
            # Determine circular ordering based on wrong split
            # The circular ordering should avoid the wrong split
            if wrong_set1 == {a, b} and wrong_set2 == {c, d}:
                # Wrong split is ab|cd, so order is a,c,b,d (avoiding ab together)
                circ = CircularOrdering([a, c, b, d])
            elif wrong_set1 == {a, c} and wrong_set2 == {b, d}:
                # Wrong split is ac|bd, so order is a,b,c,d (avoiding ac together)
                circ = CircularOrdering([a, b, c, d])
            else:  # wrong_set1 == {a, d} and wrong_set2 == {b, c}
                # Wrong split is ad|bc, so order is a,b,d,c (avoiding ad together)
                circ = CircularOrdering([a, b, d, c])
            
            # Create the two quartets from circular ordering
            # For circular ordering [a, b, c, d], quartets are:
            # (a b | c d) and (a d | b c)
            order_list = list(circ.order)
            a_ord, b_ord, c_ord, d_ord = order_list
            
            q1 = Quartet(Split({a_ord, b_ord}, {c_ord, d_ord}))
            q2 = Quartet(Split({a_ord, d_ord}, {b_ord, c_ord}))
            
            if weight:
                w = abs(lam - delta) / (1.0 - lam) if lam < 1.0 else 1.0
                quartets_dict = {q1: w, q2: w}
            else:
                quartets_dict = {q1: 1.0, q2: 1.0}
            
            profile = SqQuartetProfile(quartets_dict)
        
        profiles.append(profile)
    
    # Assign reticulation leaves to cycles based on delta sum
    # Order taxa by delta sum (descending - highest delta sum = most likely reticulation)
    reticulation_order = sorted(delta_sum.items(), key=lambda x: x[1], reverse=True)
    reticulation_order_taxa = [taxon for taxon, _ in reticulation_order]
    
    # Update profiles: assign reticulation leaves to cycles
    final_profiles: list[SqQuartetProfile] = []
    for profile in profiles:
        if len(profile.quartets) == 2:
            # This is a cycle - assign reticulation leaf
            # Find first taxon in reticulation order that is in this profile's taxa
            ret_leaf = next(
                (taxon for taxon in reticulation_order_taxa if taxon in profile.taxa),
                None
            )
            
            if ret_leaf is not None:
                # Create new profile with reticulation leaf
                final_profiles.append(
                    SqQuartetProfile(
                        dict(profile.quartets),
                        reticulation_leaf=ret_leaf
                    )
                )
            else:
                # Should not happen, but keep profile without reticulation leaf
                final_profiles.append(profile)
        else:
            # Single quartet - no reticulation leaf
            final_profiles.append(profile)
    
    return SqQuartetProfileSet(profiles=final_profiles)

