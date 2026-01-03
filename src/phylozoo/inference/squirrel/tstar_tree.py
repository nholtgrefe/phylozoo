"""
T* tree algorithm module.

This module provides functions for computing the B*-set of splits and the T* tree
from a quartet profile set.
"""

import itertools
from typing import TYPE_CHECKING

from ...core.network.sdnetwork import SemiDirectedPhyNetwork
from ...core.quartet.qprofileset import QuartetProfileSet
from ...core.split.algorithms import splitsystem_to_tree
from ...core.split.base import Split
from ...core.split.splitsystem import SplitSystem

if TYPE_CHECKING:
    pass


def bstar(profileset: QuartetProfileSet) -> SplitSystem:
    """
    Compute the B*-set of splits from a quartet profile set.
    
    This function extracts trivial (single-quartet) resolved profiles from the
    profile set and uses the incremental O(n^5) B* algorithm to compute a
    compatible set of full splits.
    
    The algorithm only considers profiles that are:
    - Trivial (contain exactly one quartet)
    - Resolved (the quartet has a non-trivial split)
    
    This essentially mimics a set of quartets with 2|2 splits.
    
    Parameters
    ----------
    profileset : QuartetProfileSet
        The quartet profile set to extract splits from.
    
    Returns
    -------
    SplitSystem
        A split system containing the B*-set of splits.
    
    Raises
    ------
    ValueError
        If there are fewer than 4 taxa in the profile set.
    
    Notes
    -----
    Uses the incremental O(n^5) algorithm from:
    Berry, Vincent, and Olivier Gascuel. "Inferring evolutionary trees with
    strong combinatorial evidence." Theoretical computer science 240.2 (2000): 271-298.
    
    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> from phylozoo.core.quartet.base import Quartet
    >>> from phylozoo.core.quartet.qprofileset import QuartetProfileSet
    >>> q1 = Quartet(Split({1, 2}, {3, 4}))
    >>> q2 = Quartet(Split({1, 3}, {2, 4}))
    >>> profileset = QuartetProfileSet(profiles=[q1, q2])
    >>> bstar_splits = bstar(profileset)
    >>> isinstance(bstar_splits, SplitSystem)
    True
    """
    # Extract all taxa
    all_taxa = list(profileset.taxa)
    
    if len(all_taxa) < 4:
        raise ValueError("B* algorithm requires at least 4 taxa")
    
    # Build a dictionary mapping 4-taxon sets to their splits
    # Only include trivial, resolved profiles
    quartet_splits: dict[frozenset[str], Split] = {}
    
    for taxa_set, (profile, _) in profileset.profiles.items():
        # Only consider trivial profiles (single quartet)
        if not profile.is_trivial():
            continue
        
        # profile.split returns None for star trees or non-trivial profiles
        split = profile.split
        if split is not None:
            quartet_splits[taxa_set] = split
    
    # Helper function to obtain split for a 4-taxon set
    def obtain_split(four_elements: frozenset[str]) -> Split | None:
        """Get the split for a given 4-taxon set, if it exists."""
        return quartet_splits.get(four_elements)
    
    # Initialize with first 4 taxa
    order = all_taxa.copy()
    a, b, c, d = order[0:4]
    
    # Start with trivial splits for first 4 taxa
    bstar_splits = [
        Split({a}, {b, c, d}),
        Split({b}, {a, c, d}),
        Split({c}, {a, b, d}),
        Split({d}, {a, b, c})
    ]
    
    # Add the quartet split for first 4 taxa if it exists
    abcd_split = obtain_split(frozenset({a, b, c, d}))
    if abcd_split is not None:
        bstar_splits.append(abcd_split)
    
    # Incrementally add each remaining taxon
    for i, element in enumerate(order):
        if i < 4:
            continue
        
        # Build set incrementally for efficiency
        previous_taxa = frozenset(order[0:i])
        new_bstar = [Split({element}, previous_taxa)]
        
        for split in bstar_splits:
            # Check compatibility before creating candidate splits
            # Check if candidate_split1 is compatible (element added to set1)
            add1 = True
            for x in split.set1:
                for y, z in itertools.combinations(split.set2, 2):
                    quartet_split = obtain_split(frozenset({element, x, y, z}))
                    if quartet_split is not None:
                        expected_split = Split({x, element}, {y, z})
                        if quartet_split != expected_split:
                            add1 = False
                            break
                if not add1:
                    break
            
            # Check if candidate_split2 is compatible (element added to set2)
            add2 = True
            for x, y in itertools.combinations(split.set1, 2):
                for z in split.set2:
                    quartet_split = obtain_split(frozenset({element, x, y, z}))
                    if quartet_split is not None:
                        expected_split = Split({x, y}, {z, element})
                        if quartet_split != expected_split:
                            add2 = False
                            break
                if not add2:
                    break
            
            # Only create candidate splits if they're compatible
            if add1:
                new_bstar.append(Split(split.set1 | {element}, split.set2))
            if add2:
                new_bstar.append(Split(split.set1, split.set2 | {element}))
        
        # Update bstar_splits after processing all splits for this element
        bstar_splits = new_bstar
    
    return SplitSystem(bstar_splits)


def tstar_tree(profileset: QuartetProfileSet) -> SemiDirectedPhyNetwork:
    """
    Compute the T* tree from a quartet profile set.
    
    The T* tree is constructed by:
    1. Computing the B*-set of splits using the bstar algorithm
    2. Converting the split system to a tree using splitsystem_to_tree
    
    Parameters
    ----------
    profileset : QuartetProfileSet
        The quartet profile set to construct the tree from.
    
    Returns
    -------
    SemiDirectedPhyNetwork
        The T* tree as a semi-directed phylogenetic network.
    
    Raises
    ------
    ValueError
        If there are fewer than 4 taxa in the profile set, or if the B*-set
        of splits is not compatible with a tree.
    
    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> from phylozoo.core.quartet.base import Quartet
    >>> from phylozoo.core.quartet.qprofileset import QuartetProfileSet
    >>> from phylozoo.core.network.sdnetwork.classifications import is_tree
    >>> q1 = Quartet(Split({1, 2}, {3, 4}))
    >>> q2 = Quartet(Split({1, 3}, {2, 4}))
    >>> profileset = QuartetProfileSet(profiles=[q1, q2])
    >>> tree = tstar_tree(profileset)
    >>> is_tree(tree)
    True
    """
    # Compute B*-set of splits, compatible by construction
    bstar_splits = bstar(profileset)
    
    # Convert split system to tree
    return splitsystem_to_tree(bstar_splits, check_compatibility=True)

