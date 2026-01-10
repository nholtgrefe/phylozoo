"""
Squirrel quartet profile set module.

This module provides the SqQuartetProfileSet class for working with
collections of squirrel quartet profiles.
"""

from typing import TYPE_CHECKING, Mapping

from ...core.quartet.base import Quartet
from ...core.quartet.qprofileset import QuartetProfileSet

from .sqprofile import SqQuartetProfile

if TYPE_CHECKING:
    pass


class SqQuartetProfileSet(QuartetProfileSet):
    """
    Squirrel quartet profile set.
    
    A SqQuartetProfileSet is a collection of SqQuartetProfile objects,
    which are quartet profiles containing either 1 or 2 resolved quartets.
    
    Parameters
    ----------
    profiles : list[SqQuartetProfile | Quartet | tuple[SqQuartetProfile, float] | tuple[Quartet, float]] | None, optional
        List of SqQuartetProfile objects, Quartet objects, or tuples with weights.
        - If SqQuartetProfile: used directly (optionally with weight tuple)
        - If Quartet: automatically grouped by taxa into SqQuartetProfile objects
        - If tuple: (profile/quartet, weight) where weight defaults to 1.0
        When quartets are provided, each profile's weight is the sum of weights
        of quartets in that profile. Each quartet within the profile keeps its input weight.
        By default None.
    taxa : frozenset[str] | None, optional
        Total set of taxa. If provided, must be a superset of all taxa in the
        profiles. Allows specifying taxa for which no profile exists.
        By default None (computed from profiles).
    
    Attributes
    ----------
    profiles : Mapping[frozenset[str], tuple[SqQuartetProfile, float]]
        Read-only mapping of 4-taxon sets to (profile, profile_weight) tuples.
        Each profile contains 1 or 2 resolved quartets with their individual weights.
    
    Raises
    ------
    ValueError
        If any profile would be empty, if any weight is non-positive, or if
        any profile is not a valid SqQuartetProfile (e.g., has unresolved quartets
        or more than 2 quartets).
    
    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> from phylozoo.core.quartet.base import Quartet
    >>> q1 = Quartet(Split({1, 2}, {3, 4}))
    >>> q2 = Quartet(Split({1, 3}, {2, 4}))
    >>> q3 = Quartet(Split({5, 6}, {7, 8}))
    >>> # From quartets (automatically creates SqQuartetProfile, merges by taxa)
    >>> sq_profileset = SqQuartetProfileSet(profiles=[q1, q2, q3])
    >>> len(sq_profileset)
    2  # q1 and q2 merged into one profile, q3 is separate
    >>> # From SqQuartetProfile objects
    >>> sq_profile1 = SqQuartetProfile([q1, q2], reticulation_leaf=1)
    >>> sq_profile2 = SqQuartetProfile([q3])
    >>> sq_profileset2 = SqQuartetProfileSet(profiles=[sq_profile1, sq_profile2])
    >>> len(sq_profileset2)
    2
    """
    
    def __init__(
        self,
        profiles: list[SqQuartetProfile | Quartet | tuple[SqQuartetProfile, float] | tuple[Quartet, float]] | None = None,
        taxa: frozenset[str] | None = None,
    ) -> None:
        """
        Initialize a squirrel quartet profile set.
        
        Parameters
        ----------
        profiles : list[SqQuartetProfile | Quartet | tuple[SqQuartetProfile, float] | tuple[Quartet, float]] | None, optional
            List of SqQuartetProfile objects, Quartet objects, or tuples with weights.
            - If SqQuartetProfile: used directly (optionally with weight tuple)
            - If Quartet: automatically grouped by taxa into SqQuartetProfile objects
            - If tuple: (profile/quartet, weight) where weight defaults to 1.0
            When quartets are provided, each profile's weight is the sum of weights
            of quartets in that profile. Each quartet within the profile keeps its input weight.
            By default None.
        taxa : frozenset[str] | None, optional
            Total set of taxa. If provided, must be a superset of all taxa in the
            profiles. Allows specifying taxa for which no profile exists.
            By default None (computed from profiles).
        
        Raises
        ------
        ValueError
            If any profile would be empty, if any weight is non-positive, or if
            any profile is not a valid SqQuartetProfile (e.g., has unresolved quartets
            or more than 2 quartets).
        """
        if profiles is None:
            profiles = []
        
        # Separate SqQuartetProfile objects from Quartet objects
        sq_profiles: list[tuple[SqQuartetProfile, float]] = []
        quartet_data: dict[frozenset[str], dict[Quartet, float]] = {}
        
        for item in profiles:
            if isinstance(item, tuple):
                obj, weight = item
            else:
                obj = item
                weight = 1.0
            
            if isinstance(obj, SqQuartetProfile):
                # Already a SqQuartetProfile, use as-is
                sq_profiles.append((obj, weight))
            elif isinstance(obj, Quartet):
                # Group quartets by taxa for merging
                quartet_taxa = obj.taxa
                if quartet_taxa not in quartet_data:
                    quartet_data[quartet_taxa] = {}
                # Check for duplicate quartets during iteration
                if obj in quartet_data[quartet_taxa]:
                    raise ValueError(
                        f"Quartet {obj} appears multiple times in the input. "
                        "Each quartet can only appear once per taxa set."
                    )
                quartet_data[quartet_taxa][obj] = weight
            else:
                raise ValueError(
                    f"Expected SqQuartetProfile or Quartet, got {type(obj)}"
                )
        
        # Convert grouped quartets to SqQuartetProfile objects
        converted_profiles: list[tuple[SqQuartetProfile, float]] = []
        
        # Add pre-existing SqQuartetProfile objects
        converted_profiles.extend(sq_profiles)
        
        # Process grouped quartets: merge quartets with same taxa into profiles
        for taxa_set, quartets_dict in quartet_data.items():
            if len(quartets_dict) == 0:
                raise ValueError(f"Cannot have empty profile for taxa {taxa_set}")
            
            # Create SqQuartetProfile from all quartets with this taxa
            # This will validate that we have 1-2 resolved quartets
            sq_profile = SqQuartetProfile(quartets_dict)
            
            # Profile weight = sum of quartet weights
            profile_weight = sum(quartets_dict.values())
            converted_profiles.append((sq_profile, profile_weight))
        
        # Initialize parent with converted profiles
        super().__init__(profiles=converted_profiles, taxa=taxa)
        
        # Validate that all profiles are SqQuartetProfile instances
        for profile, _ in self._profiles.values():
            if not isinstance(profile, SqQuartetProfile):
                raise ValueError(
                    f"All profiles must be SqQuartetProfile instances, "
                    f"got {type(profile)}"
                )
    
    @property
    def profiles(self) -> Mapping[frozenset[str], tuple[SqQuartetProfile, float]]:
        """
        Get the profiles and their weights (read-only).
        
        Returns
        -------
        Mapping[frozenset[str], tuple[SqQuartetProfile, float]]
            Read-only mapping of 4-taxon sets to (profile, profile_weight) tuples.
        """
        # Return the parent's profiles property (already validated to be SqQuartetProfile)
        return self._profiles  # type: ignore
    
    def get_profile(self, taxa: frozenset[str]) -> SqQuartetProfile | None:
        """
        Get the SqQuartetProfile for a 4-taxon set.
        
        Parameters
        ----------
        taxa : frozenset[str]
            The 4-taxon set.
        
        Returns
        -------
        SqQuartetProfile | None
            The profile for the taxa, or None if not found.
        """
        result = self._profiles.get(taxa)
        return result[0] if result else None  # type: ignore
    
    def __repr__(self) -> str:
        """
        Return string representation of the profile set that can be used to initialize it.
        
        Returns
        -------
        str
            String representation that can be used to recreate the object.
        """
        if len(self._profiles) == 0:
            return "SqQuartetProfileSet(profiles={})"
        
        # Build list of (profile, weight) tuples for initialization
        # Use the quartets dict directly since that's what __init__ accepts
        profile_items = []
        for taxa, (profile, weight) in self._profiles.items():
            # Get the quartets dict from the profile
            quartets_dict = dict(profile.quartets)
            ret_leaf_str = f", reticulation_leaf='{profile.reticulation_leaf}'" if profile.reticulation_leaf else ""
            profile_items.append(
                f"(SqQuartetProfile({repr(quartets_dict)}{ret_leaf_str}), {weight})"
            )
        
        profiles_str = ", ".join(profile_items)
        return f"SqQuartetProfileSet(profiles=[{profiles_str}])"

