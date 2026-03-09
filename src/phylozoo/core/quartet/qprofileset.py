"""
Quartet profile set module.

A quartet profile set is a collection of quartet profiles covering multiple four-taxon sets.
This module provides the QuartetProfileSet class.
"""

from functools import cached_property
from math import comb
from types import MappingProxyType
from typing import Iterator, Mapping, TYPE_CHECKING

from ...utils.exceptions import PhyloZooValueError
from .base import Quartet
from .qprofile import QuartetProfile

if TYPE_CHECKING:
    pass


class QuartetProfileSet:
    """
    Immutable collection of quartet profiles with two-level weights.
    
    A QuartetProfileSet groups quartets by their 4-taxon sets into profiles.
    Each profile has a weight (profile weight), and each quartet within a
    profile also has a weight (quartet weight).
    
    This allows representing uncertainty or multiple hypotheses about quartet
    topologies for the same set of 4 taxa, with different weights assigned
    to each hypothesis.
    
    Parameters
    ----------
    profiles : list[QuartetProfile | Quartet | tuple[QuartetProfile, float]] | None, optional
        List of QuartetProfile objects, Quartet objects, or tuples with profile weights.
        
        - If QuartetProfile: used directly (optionally with a profile-weight tuple).
        - If Quartet: automatically grouped by taxa into profiles. For each 4-taxon
          set, all quartets on that set are collected into a :class:`QuartetProfile`
          with equal weights :math:`1/k` (where :math:`k` is the number of quartets
          for that taxa set). Each resulting profile in the set receives default
          profile weight 1.0.
        - If tuple: (profile, weight) where weight is the profile weight.
        
        Passing quartets together with explicit weights (e.g. ``(Quartet, weight)``)
        is not supported. To use non-uniform quartet weights within a profile,
        construct a :class:`QuartetProfile` explicitly and pass that (optionally
        with a profile weight).
        
        By default None.
    taxa : frozenset[str] | None, optional
        Total set of taxa. If provided, must be a superset of all taxa in the
        profiles. Allows specifying taxa for which no profile exists.
        By default None (computed from profiles).
    
    Raises
    ------
    PhyloZooValueError
        If any profile would be empty, if any weight is non-positive, if profiles/quartets
        are mixed incorrectly, or if provided taxa is not a superset of profile taxa.
    
    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> q1 = Quartet(Split({1, 2}, {3, 4}))
    >>> q2 = Quartet(Split({1, 3}, {2, 4}))
    >>> q3 = Quartet(Split({5, 6}, {7, 8}))
    
    >>> # From quartets (grouped into profiles, equal weights per profile)
    >>> profileset = QuartetProfileSet(profiles=[q1, q2, q3])
    >>> len(profileset)
    2
    >>> profileset.get_profile_weight(frozenset({1, 2, 3, 4}))
    1.0
    
    >>> # From QuartetProfile objects (better control)
    >>> profile1 = QuartetProfile({q1: 0.8, q2: 0.2})
    >>> profile2 = QuartetProfile([q3])
    >>> profileset2 = QuartetProfileSet(profiles=[(profile1, 2.0), (profile2, 1.5)])
    >>> profileset2.get_profile_weight(frozenset({1, 2, 3, 4}))
    2.0
    
    Attributes
    ----------
    profiles : Mapping[frozenset[str], tuple[QuartetProfile, float]]
        Read-only mapping of 4-taxon sets to (profile, profile_weight) tuples.
        Each profile contains multiple quartets with their individual weights.
    """
    
    def __init__(
        self,
        profiles: list[QuartetProfile | Quartet | tuple[QuartetProfile, float]] | None = None,
        taxa: frozenset[str] | None = None,
    ) -> None:
        """
        Initialize a quartet profile set.
        
        Parameters
        ----------
        profiles : list[QuartetProfile | Quartet | tuple[QuartetProfile, float]] | None, optional
            List of QuartetProfile objects, Quartet objects, or tuples with profile weights.
            
            - If QuartetProfile: used directly (optionally with profile-weight tuple)
            - If Quartet: automatically grouped by taxa into profiles. Each resulting
              profile receives default weight 1.0 and uses :class:`QuartetProfile`
              to handle internal quartet weights (equal weights 1/k when created
              from a list of quartets).
            - If tuple: (profile, weight) where weight is the profile weight and
              must be positive.
            
            Passing quartets with explicit weights (e.g. ``(Quartet, weight)``) is
            not supported. To specify non-uniform quartet weights, construct a
            :class:`QuartetProfile` explicitly and pass that instead.
            
            By default None.
        taxa : frozenset[str] | None, optional
            Total set of taxa. If provided, must be a superset of all taxa in the
            profiles. Allows specifying taxa for which no profile exists.
            By default None (computed from profiles).
        
        Raises
        ------
        ValueError
            If any profile would be empty, if any weight is non-positive, or if
            the provided taxa is not a superset of all taxa in profiles.
        """
        if profiles is None:
            profiles = []
        
        all_taxa_from_input: set[str] = set()
        profiles_dict: dict[frozenset[str], tuple[QuartetProfile, float]] = {}
        
        # Process items in a single pass, determining mode from first item
        mode: str | None = None  # 'profile' or 'quartet'
        profile_data: dict[frozenset[str], dict[Quartet, float]] = {}
        
        for item in profiles:
            # Extract object and weight (weights only allowed for QuartetProfile)
            if isinstance(item, tuple):
                obj, weight = item
                if isinstance(obj, Quartet):
                    raise PhyloZooValueError(
                        "Quartet weights are not supported in QuartetProfileSet. "
                        "Construct a QuartetProfile with quartet weights and pass "
                        "that instead."
                    )
            else:
                obj = item
                weight = 1.0
            
            # Determine mode from first item, validate consistency for subsequent items
            if mode is None:
                if isinstance(obj, QuartetProfile):
                    mode = 'profile'
                elif isinstance(obj, Quartet):
                    mode = 'quartet'
                else:
                    raise PhyloZooValueError(f"Expected QuartetProfile or Quartet, got {type(obj)}")
            else:
                # Validate that subsequent items match the mode
                if mode == 'profile' and not isinstance(obj, QuartetProfile):
                    if isinstance(obj, Quartet):
                        raise PhyloZooValueError("Cannot mix QuartetProfile and Quartet objects in profiles list")
                    raise PhyloZooValueError(f"Expected QuartetProfile, got {type(obj)}")
                elif mode == 'quartet' and not isinstance(obj, Quartet):
                    if isinstance(obj, QuartetProfile):
                        raise PhyloZooValueError("Cannot mix QuartetProfile and Quartet objects in profiles list")
                    raise PhyloZooValueError(f"Expected Quartet, got {type(obj)}")
            
            # Validate weight
            if weight <= 0:
                obj_type = "profile" if mode == 'profile' else "quartet"
                raise PhyloZooValueError(
                    f"{obj_type.capitalize()} weight must be positive, got {weight} for {obj_type} {obj}"
                )
            
            if mode == 'profile':
                # Mode 1: From QuartetProfile objects
                profile = obj
                profile_taxa = profile.taxa
                all_taxa_from_input.update(profile_taxa)
                # Check for duplicate taxa sets
                if profile_taxa in profiles_dict:
                    raise PhyloZooValueError(
                        f"Multiple profiles with the same taxa set {profile_taxa} are not allowed. "
                        "Each 4-taxon set can only have one profile."
                    )
                profiles_dict[profile_taxa] = (profile, weight)
            
            else:
                # Mode 2: From Quartet objects (group into profiles, no quartet weights)
                quartet = obj
                quartet_taxa = quartet.taxa
                all_taxa_from_input.update(quartet_taxa)
                
                if quartet_taxa not in profile_data:
                    profile_data[quartet_taxa] = {}
                # Check for duplicate quartets during iteration
                if quartet in profile_data[quartet_taxa]:
                    raise PhyloZooValueError(
                        f"Quartet {quartet} appears multiple times in the input. "
                        "Each quartet can only appear once per taxa set."
                    )
                # Value is unused beyond duplicate detection; all quartets are unweighted here.
                profile_data[quartet_taxa][quartet] = 1.0
        
        # If we processed quartets, create profiles from grouped data.
        # Each taxa set becomes a QuartetProfile built from its quartets, with
        # default profile weight 1.0 in the set.
        if mode == 'quartet':
            for taxa_set, quartets_dict in profile_data.items():
                # Validate no empty profiles
                if len(quartets_dict) == 0:
                    raise PhyloZooValueError(f"Cannot have empty profile for taxa {taxa_set}")
                
                profile = QuartetProfile(list(quartets_dict.keys()))
                profile_weight = 1.0
                profiles_dict[taxa_set] = (profile, profile_weight)
        
        # Store as immutable
        self._profiles = MappingProxyType(profiles_dict)
        
        # Handle taxa parameter
        if taxa is not None:
            # Validate that provided taxa is a superset
            if not all_taxa_from_input.issubset(taxa):
                missing = all_taxa_from_input - taxa
                raise PhyloZooValueError(
                    f"Provided taxa must be a superset of all taxa in profiles/quartets. "
                    f"Missing taxa: {missing}"
                )
            self._taxa = frozenset(taxa)
        else:
            # Compute all taxa from profiles
            self._taxa = frozenset(all_taxa_from_input)
    
    @property
    def profiles(self) -> Mapping[frozenset[str], tuple[QuartetProfile, float]]:
        """
        Get the profiles and their weights (read-only).
        
        Returns
        -------
        Mapping[frozenset[str], tuple[QuartetProfile, float]]
            Read-only mapping of 4-taxon sets to (profile, profile_weight) tuples.
        """
        return self._profiles
    
    @property
    def taxa(self) -> frozenset[str]:
        """
        Get all taxa in the profile set.
        
        Returns
        -------
        frozenset[str]
            Set of all taxon labels. If taxa was specified during initialization,
            this includes all specified taxa (even those without profiles).
            Otherwise, returns taxa that appear in at least one profile.
        """
        return self._taxa
    
    @cached_property
    def is_dense(self) -> bool:
        """
        Check if the quartet profile set is dense.
        
        A dense quartet profile set has a profile for every possible combination
        of 4 taxa from the total set of taxa.
        
        Returns
        -------
        bool
            True if the profile set is dense (has C(n, 4) profiles where n is
            the number of taxa), False otherwise.
        """
        n = len(self._taxa)
        if n < 4:
            return len(self._profiles) == 0
        return len(self._profiles) == comb(n, 4)
    
    @cached_property
    def is_all_resolved(self) -> bool:
        """
        Check if all profiles in the set are resolved.
        
        A profile is resolved if all its quartets are resolved (i.e., not star trees).
        This property returns True only if every profile in the set is resolved.
        
        Returns
        -------
        bool
            True if all profiles are resolved, False otherwise.
        
        Examples
        --------
        >>> from phylozoo.core.split.base import Split
        >>> from phylozoo.core.quartet.base import Quartet
        >>> q1 = Quartet(Split({1, 2}, {3, 4}))
        >>> q2 = Quartet(Split({1, 3}, {2, 4}))
        >>> profileset = QuartetProfileSet(profiles=[q1, q2])
        >>> profileset.is_all_resolved
        True
        >>> star = Quartet({1, 2, 3, 4})
        >>> profileset2 = QuartetProfileSet(profiles=[q1, star])
        >>> profileset2.is_all_resolved
        False
        """
        return all(profile.is_resolved() for profile, _ in self._profiles.values())
    
    @cached_property
    def max_profile_len(self) -> int:
        """
        Get the maximum number of quartets in any profile.
        
        Returns
        -------
        int
            The largest number of quartets in any profile in the set.
            Returns 0 if the set is empty.
        
        Examples
        --------
        >>> from phylozoo.core.split.base import Split
        >>> from phylozoo.core.quartet.base import Quartet
        >>> q1 = Quartet(Split({1, 2}, {3, 4}))
        >>> q2 = Quartet(Split({1, 3}, {2, 4}))
        >>> q3 = Quartet(Split({1, 4}, {2, 3}))
        >>> profileset = QuartetProfileSet(profiles=[q1, q2, q3])
        >>> profileset.max_profile_len
        3
        """
        if len(self._profiles) == 0:
            return 0
        return max(len(profile) for profile, _ in self._profiles.values())
    
    def get_profile(self, taxa: frozenset[str]) -> QuartetProfile | None:
        """
        Get the profile for a 4-taxon set.
        
        Parameters
        ----------
        taxa : frozenset[str]
            The 4-taxon set.
        
        Returns
        -------
        QuartetProfile | None
            The profile for the taxa, or None if not found.
        """
        result = self._profiles.get(taxa)
        return result[0] if result else None
    
    def get_profile_weight(self, taxa: frozenset[str]) -> float | None:
        """
        Get the profile weight for a 4-taxon set.
        
        Parameters
        ----------
        taxa : frozenset[str]
            The 4-taxon set.
        
        Returns
        -------
        float | None
            The profile weight, or None if not found.
        """
        result = self._profiles.get(taxa)
        return result[1] if result else None
    
    def has_profile(self, taxa: frozenset[str]) -> bool:
        """
        Check if a profile exists for the given 4-taxon set.
        
        Parameters
        ----------
        taxa : frozenset[str]
            The 4-taxon set.
        
        Returns
        -------
        bool
            True if a profile exists, False otherwise.
        """
        return taxa in self._profiles
    
    def all_profile_taxon_sets(self) -> Iterator[frozenset[str]]:
        """
        Get all 4-taxon sets that have profiles.
        
        Returns
        -------
        Iterator[frozenset[str]]
            Iterator over all 4-taxon sets that have profiles.
        """
        return iter(self._profiles.keys())
    
    def __iter__(self) -> Iterator[tuple[QuartetProfile, float]]:
        """
        Return an iterator over (profile, profile_weight) pairs.
        
        Returns
        -------
        Iterator[tuple[QuartetProfile, float]]
            Iterator over (profile, weight) tuples.
        """
        return iter(self._profiles.values())
    
    def __len__(self) -> int:
        """
        Return the number of profiles.
        
        Returns
        -------
        int
            Number of profiles.
        """
        return len(self._profiles)
    
    def __contains__(self, taxa: frozenset[str]) -> bool:
        """
        Check if a profile exists for the given 4-taxon set.
        
        Parameters
        ----------
        taxa : frozenset[str]
            The 4-taxon set.
        
        Returns
        -------
        bool
            True if a profile exists, False otherwise.
        """
        return taxa in self._profiles
    
    def __repr__(self) -> str:
        """
        Return string representation of the profile set that can be used to initialize it.
        
        Returns
        -------
        str
            String representation that can be used to recreate the object.
        """
        if len(self._profiles) == 0:
            return "QuartetProfileSet(profiles={})"
        
        # Build list of (profile, weight) tuples for initialization
        # Use the quartets dict directly since that's what __init__ accepts
        profile_items = []
        for taxa, (profile, weight) in self._profiles.items():
            # Get the quartets dict from the profile
            quartets_dict = dict(profile.quartets)
            profile_items.append(f"(QuartetProfile({repr(quartets_dict)}), {weight})")
        
        profiles_str = ", ".join(profile_items)
        return f"QuartetProfileSet(profiles=[{profiles_str}])"
    
    def __str__(self) -> str:
        """
        Return human-readable string representation of the quartet profile set.
        
        Displays one line per profile, showing the profile (using its __str__ method)
        and its profile weight. Aligns with QuartetProfile's __str__ format.
        
        Returns
        -------
        str
            Human-readable string representation.
        
        Examples
        --------
        >>> from phylozoo.core.split.base import Split
        >>> q1 = Quartet(Split({1, 2}, {3, 4}))
        >>> q2 = Quartet(Split({1, 3}, {2, 4}))
        >>> profile1 = QuartetProfile({q1: 0.8, q2: 0.2})
        >>> profileset = QuartetProfileSet(profiles=[(profile1, 1.0)])
        >>> str(profileset)
        'QuartetProfileSet({\\n  QuartetProfile({...}) [weight: 1.0]\\n})'
        """
        if len(self._profiles) == 0:
            return "QuartetProfileSet({})"
        
        # Sort profiles by taxa for consistent display
        sorted_profiles = sorted(
            self._profiles.items(),
            key=lambda item: sorted(item[0])
        )
        
        # Show all profiles with weights, one per line
        profile_lines = []
        for taxa, (profile, weight) in sorted_profiles:
            # Use profile's __str__ directly, followed by weight
            profile_str = str(profile)
            profile_lines.append(f"  {profile_str} [weight: {weight}],")
        
        # Remove trailing comma from last line
        if profile_lines:
            profile_lines[-1] = profile_lines[-1].rstrip(',')
        
        return f"QuartetProfileSet({{\n" + "\n".join(profile_lines) + "\n})"



