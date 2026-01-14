"""
Quartet distance module.

This module provides functions for computing distance matrices from quartet profiles
based on the quartet distance metric.

"""

import itertools
from typing import TYPE_CHECKING

import numpy as np

from ...utils.exceptions import PhyloZooValueError
from ..distance import DistanceMatrix

if TYPE_CHECKING:
    from .qprofile import QuartetProfile
    from .qprofileset import QuartetProfileSet
    from ..primitives.partition import Partition


def quartet_distance(
    profileset: 'QuartetProfileSet',
    rho: tuple[float, float, float, float],
) -> DistanceMatrix:
    """
    Compute a distance matrix between taxa based on quartet profiles.
    
    This function computes pairwise distances between taxa by aggregating contributions
    from all quartet profiles. The distance is based on the quartet topology and uses
    a rho vector to weight different quartet types.
    
    Parameters
    ----------
    profileset : QuartetProfileSet
        The quartet profile set to compute distances from. Must be dense.
        Each profile must contain exactly 1 or 2 resolved quartets.
    rho : tuple[float, float, float, float]
        Rho vector (rho_c, rho_s, rho_a, rho_o) specifying distance contributions
        for different quartet topologies:
        - rho_c: contribution when two leaves are on different sides of a split
          (used for profiles with 1 quartet)
        - rho_s: contribution when two leaves are on the same side of a split
          (used for profiles with 1 quartet)
        - rho_a: contribution when two leaves are adjacent in a circular ordering
          (used for profiles with 2 quartets, which induce a circular ordering)
        - rho_o: contribution when two leaves are opposite in a circular ordering
          (used for profiles with 2 quartets, which induce a circular ordering)
        Note: For each quartet profile containing a pair of leaves, 2*rho_xy is
        added to the distance between leaves x and y, where rho_xy is the
        appropriate rho value based on the profile type and leaf positions.
    
    Returns
    -------
    DistanceMatrix
        A distance matrix with pairwise distances between taxa.
        The matrix is symmetric with zero diagonal.
    
    Raises
    ------
    ValueError
        If rho vector has invalid length or values.
        If the profile set is not dense.
        If any profile contains more than 2 quartets or contains unresolved quartets.
    
    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> from phylozoo.core.quartet.base import Quartet
    >>> from phylozoo.core.quartet.qprofileset import QuartetProfileSet
    >>> # Create a simple profile set
    >>> q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
    >>> q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
    >>> profileset = QuartetProfileSet(profiles=[q1, q2])
    >>> # Compute distance matrix with Squirrel rho vector
    >>> dist_matrix = quartet_distance(profileset, rho=(0.5, 1.0, 0.5, 1.0))
    >>> # Access distances
    >>> dist_matrix.get_distance('A', 'B')  # Distance between A and B
    1.0
    
    Notes
    -----
    Common rho vector values:
    - Squirrel/MONAD: (0.5, 1.0, 0.5, 1.0)
    - NANUQ: (0.0, 1.0, 0.5, 1.0)
    """
    # Validate rho vector
    if len(rho) != 4:
        raise PhyloZooValueError(f"Rho vector must have 4 elements, got {len(rho)}")
    rho_c, rho_s, rho_a, rho_o = rho
    
    # Validate rho constraints
    if rho_a > rho_o or rho_c > rho_s:
        raise PhyloZooValueError(
            "Rho vector must satisfy: rho_a <= rho_o and rho_c <= rho_s"
        )
    
    # Check if profile set is dense
    if not profileset.is_dense:
        raise PhyloZooValueError("Profile set must be dense (have a profile for every 4-taxon combination)")
    
    # Get all taxa and initialize distance matrix as numpy array
    taxa = sorted(profileset.taxa)
    n = len(taxa)
    D = np.zeros((n, n), dtype=np.float64)
    
    # Create mapping from taxon to index
    taxon_to_index = {taxon: i for i, taxon in enumerate(taxa)}
    
    # Iterate over all 4-taxon combinations
    for four_taxa in itertools.combinations(taxa, 4):
        four_taxa_set = frozenset(four_taxa)
        profile = profileset.get_profile(four_taxa_set)
        
        if profile is None:
            raise PhyloZooValueError(
                f"No profile found for 4-taxon set {four_taxa_set}. "
                "Profile set must be dense."
            )
        
        # Validate profile: must have 1 or 2 quartets, all resolved
        num_quartets = len(profile.quartets)
        if num_quartets == 0:
            raise PhyloZooValueError(
                f"Profile for {four_taxa_set} has no quartets"
            )
        if num_quartets > 2:
            raise PhyloZooValueError(
                f"Profile for {four_taxa_set} has {num_quartets} quartets. "
                "Each profile must contain exactly 1 or 2 quartets."
            )
        
        # Check all quartets are resolved
        for quartet in profile.quartets:
            if not quartet.is_resolved():
                raise PhyloZooValueError(
                    f"Profile for {four_taxa_set} contains unresolved quartet (star tree). "
                    "All quartets must be resolved."
                )
        
        # Compute rho-distance for each pair of leaves using the profile
        for leaf1, leaf2 in itertools.combinations(four_taxa, 2):
            rho_dist = _rho_distance(profile, leaf1, leaf2, rho)
            
            # Get indices
            i = taxon_to_index[leaf1]
            j = taxon_to_index[leaf2]
            
            delta = 2 * rho_dist
            D[i, j] += delta
            D[j, i] += delta  # Symmetric matrix
    
    # Add constant 2*n - 4 to all off-diagonal entries
    constant = 2 * n - 4
    D += constant
    np.fill_diagonal(D, 0)  # Keep diagonal at 0
    
    # Create and return DistanceMatrix
    return DistanceMatrix(D, labels=taxa)


def _rho_distance(
    profile: 'QuartetProfile',
    leaf1: str,
    leaf2: str,
    rho: tuple[float, float, float, float],
) -> float:
    """
    Compute rho-distance between two leaves in a quartet profile.
    
    If the profile has 1 quartet (split quarnet), uses split-based logic:
    - Same side of split: rho_s
    - Different sides: rho_c
    
    If the profile has 2 quartets (four-cycle), uses circular ordering logic:
    - Adjacent in ordering: rho_a
    - Opposite in ordering: rho_o
    
    Parameters
    ----------
    profile : QuartetProfile
        The quartet profile. Must have 1 or 2 resolved quartets.
    leaf1 : str
        First leaf.
    leaf2 : str
        Second leaf.
    rho : tuple[float, float, float, float]
        Rho vector (rho_c, rho_s, rho_a, rho_o).
    
    Returns
    -------
    float
        The rho-distance.
    """
    rho_c, rho_s, rho_a, rho_o = rho
    num_quartets = len(profile.quartets)
    
    if num_quartets == 1:
        # Single quartet: use split-based logic
        quartet = next(iter(profile.quartets))
        split = quartet.split
        if split is None:
            raise PhyloZooValueError("Quartet must be resolved (have a split)")
        
        # Check if both leaves are on the same side of the split
        same_side = (leaf1 in split.set1 and leaf2 in split.set1) or (
            leaf1 in split.set2 and leaf2 in split.set2
        )
        
        if same_side:
            return float(rho_s)
        else:
            return float(rho_c)
    
    elif num_quartets == 2:
        # Two quartets: use circular ordering logic
        circular_orderings = profile.circular_orderings
        if circular_orderings is None or len(circular_orderings) == 0:
            raise PhyloZooValueError(
                "Profile with 2 quartets must have a circular ordering"
            )
        
        # Use the first circular ordering (they should all be equivalent)
        ordering = next(iter(circular_orderings))
        
        # Check if leaves are adjacent or opposite in the circular ordering
        if ordering.are_neighbors(leaf1, leaf2):
            # Adjacent: use rho_a
            return float(rho_a)
        else:
            # Opposite: use rho_o
            return float(rho_o)
    
    else:
        raise PhyloZooValueError(
            f"Profile must have 1 or 2 quartets, got {num_quartets}"
        )


def quartet_distance_with_partition(
    profileset: 'QuartetProfileSet',
    partition: 'Partition',
    rho: tuple[float, float, float, float] = (0.5, 1.0, 0.5, 1.0),
) -> DistanceMatrix:
    """
    Compute a distance matrix between partition sets based on quartet profiles.
    
    This function computes pairwise distances between sets in a partition by
    aggregating contributions from quartet profiles. For each 4-subpartition,
    it considers all quartets with one leaf from each set and aggregates
    their contributions.
    
    Parameters
    ----------
    profileset : QuartetProfileSet
        The quartet profile set to compute distances from. Must be dense.
        Each profile must contain exactly 1 or 2 resolved quartets.
    partition : Partition
        A partition of the taxa. The distance matrix will be computed
        between the sets (parts) of this partition.
    rho : tuple[float, float, float, float], optional
        Rho vector (rho_c, rho_s, rho_a, rho_o) specifying distance contributions.
        By default (0.5, 1.0, 0.5, 1.0) (Squirrel/MONAD).
    
    Returns
    -------
    DistanceMatrix
        A distance matrix with pairwise distances between partition sets.
        The matrix is symmetric with zero diagonal.
        Labels are the partition sets (frozensets).
    
    Raises
    ------
    PhyloZooValueError
        If rho vector has invalid length or values.
        If partition elements don't match profile set taxa.
        If profile set is not dense.
        If any profile contains more than 2 quartets or contains unresolved quartets.
    
    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> from phylozoo.core.quartet.base import Quartet
    >>> from phylozoo.core.quartet.qprofileset import QuartetProfileSet
    >>> from phylozoo.core.primitives.partition import Partition
    >>> 
    >>> # Create a simple profile set
    >>> q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
    >>> profileset = QuartetProfileSet(profiles=[q1])
    >>> 
    >>> # Create a partition
    >>> partition = Partition([{'A'}, {'B'}, {'C'}, {'D'}])
    >>> 
    >>> # Compute distance matrix
    >>> dist_matrix = quartet_distance_with_partition(profileset, partition)
    >>> len(dist_matrix)
    4
    """
    # Validate rho vector
    if len(rho) != 4:
        raise PhyloZooValueError(f"Rho vector must have 4 elements, got {len(rho)}")
    rho_c, rho_s, rho_a, rho_o = rho
    
    # Validate rho constraints
    if rho_a > rho_o or rho_c > rho_s:
        raise PhyloZooValueError(
            "Rho vector must satisfy: rho_a <= rho_o and rho_c <= rho_s"
        )
    
    # Validate partition matches profile set taxa
    if partition.elements != profileset.taxa:
        raise PhyloZooValueError(
            "Partition elements must match profile set taxa. "
            f"Partition has {partition.elements}, profile set has {profileset.taxa}"
        )
    
    # Check if profile set is dense
    if not profileset.is_dense:
        raise PhyloZooValueError("Profile set must be dense (have a profile for every 4-taxon combination)")
    
    # Get partition size and set order
    n = len(partition)  # Number of sets in partition
    set_order = list(partition.parts)
    
    # Pre-compute set-to-index mapping (O(n) once instead of O(n) per lookup)
    set_to_index: dict[frozenset, int] = {part: idx for idx, part in enumerate(set_order)}
    
    # Initialize distance matrix
    D = np.zeros((n, n), dtype=np.float64)
    
    # Iterate over subpartitions containing four sets of the partition
    for four_sub_partition in partition.subpartitions(4):
        # Map sets to indices using pre-computed dictionary (O(1) lookup)
        X_indices: dict[frozenset, int] = {
            X: set_to_index[X] for X in four_sub_partition.parts
        }
        
        # Pre-compute leaf-to-set mapping for this 4-subpartition (O(1) lookup instead of O(4))
        leaf_to_set: dict[str, frozenset] = {}
        for part in four_sub_partition.parts:
            for leaf in part:
                leaf_to_set[leaf] = part
        
        # Dictionary to collect distances for each pair of sets
        # Key: tuple of two set indices (faster than frozenset), Value: list of distance contributions
        d_lists: dict[tuple[int, int], list[float]] = {}
        indices_list = list(X_indices.values())
        for i, j in itertools.combinations(indices_list, 2):
            # Use tuple instead of frozenset for slightly faster hashing
            d_lists[(i, j)] = []
        
        # Iterate over all quartets on the four sets of the partition
        # For each 4-subpartition, get all representative partitions (one leaf per set)
        for four_leaf_partition in four_sub_partition.representative_partitions():
            four_taxa_set = four_leaf_partition.elements
            
            # Get profile for this 4-taxon set
            profile = profileset.get_profile(four_taxa_set)
            if profile is None:
                raise PhyloZooValueError(
                    f"No profile found for 4-taxon set {four_taxa_set}. "
                    "Profile set must be dense."
                )
            
            # Validate profile: must have 1 or 2 quartets, all resolved
            num_quartets = len(profile.quartets)
            if num_quartets == 0:
                raise PhyloZooValueError(
                    f"Profile for {four_taxa_set} has no quartets"
                )
            if num_quartets > 2:
                raise PhyloZooValueError(
                    f"Profile for {four_taxa_set} has {num_quartets} quartets. "
                    "Each profile must contain exactly 1 or 2 quartets."
                )
            
            # Check all quartets are resolved
            for quartet in profile.quartets:
                if not quartet.is_resolved():
                    raise PhyloZooValueError(
                        f"Profile for {four_taxa_set} contains unresolved quartet. "
                        "All quartets must be resolved."
                    )
            
            # Compute rho-distance for each pair of leaves
            for leaf1, leaf2 in itertools.combinations(four_taxa_set, 2):
                # Get rho-distance (same for all quartets in the profile if they share structure)
                rho_dist = _rho_distance(profile, leaf1, leaf2, rho)
                
                # Find which sets these leaves belong to (O(1) lookup)
                Xi = leaf_to_set[leaf1]
                Xj = leaf_to_set[leaf2]
                
                # Get indices (O(1) lookup)
                i = X_indices[Xi]
                j = X_indices[Xj]
                
                # Ensure i < j for consistent tuple key ordering
                if i > j:
                    i, j = j, i
                
                # Add to distance list for this pair of sets
                delta = 2 * rho_dist
                d_lists[(i, j)].append(delta)
        
        # Compute the average distance between the four sets
        for (i, j), value in d_lists.items():
            if len(value) == 0:
                continue  # No contributions for this pair
            
            # Compute mean (average) of all contributions
            delta = sum(value) / len(value)
            D[i, j] += delta
            D[j, i] += delta  # Symmetric matrix
    
    # Add constant 2*n - 4 to all off-diagonal entries (same as quartet_distance)
    constant = 2 * n - 4
    D += constant
    np.fill_diagonal(D, 0)
    
    # Create DistanceMatrix with partition sets as labels
    return DistanceMatrix(D, labels=set_order)

