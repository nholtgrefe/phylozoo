"""
Triplet profile set module.

A triplet profile set is a collection of triplet profiles covering multiple three-taxon sets.
This module provides the TripletProfileSet class.
"""

from functools import cached_property
from math import comb
from types import MappingProxyType
from typing import Iterator, Mapping, TYPE_CHECKING

from ...utils.exceptions import PhyloZooValueError
from .base import Triplet
from .tprofile import TripletProfile

if TYPE_CHECKING:
    pass


class TripletProfileSet:
    """
    Immutable collection of triplet profiles with two-level weights.

    A TripletProfileSet groups triplets by their 3-taxon sets into profiles.
    Each profile has a weight (profile weight), and each triplet within a
    profile also has a weight (triplet weight).

    This allows representing uncertainty or multiple hypotheses about triplet
    topologies for the same set of 3 taxa, with different weights assigned
    to each hypothesis.

    Parameters
    ----------
    profiles : list[TripletProfile | Triplet | tuple[TripletProfile, float]] | None, optional
        List of TripletProfile objects, Triplet objects, or tuples with profile weights.

        - If TripletProfile: used directly (optionally with a profile-weight tuple).
        - If Triplet: automatically grouped by taxa into profiles. For each 3-taxon
          set, all triplets on that set are collected into a :class:`TripletProfile`
          with equal weights :math:`1/k` (where :math:`k` is the number of triplets
          for that taxa set). Each resulting profile in the set receives default
          profile weight 1.0.
        - If tuple: (profile, weight) where weight is the profile weight.

        Passing triplets together with explicit weights (e.g. ``(Triplet, weight)``)
        is not supported. To use non-uniform triplet weights within a profile,
        construct a :class:`TripletProfile` explicitly and pass that (optionally
        with a profile weight).

        By default None.
    taxa : frozenset[str] | None, optional
        Total set of taxa. If provided, must be a superset of all taxa in the
        profiles. Allows specifying taxa for which no profile exists.
        By default None (computed from profiles).

    Raises
    ------
    PhyloZooValueError
        If any profile would be empty, if any weight is non-positive, if profiles/triplets
        are mixed incorrectly, or if provided taxa is not a superset of profile taxa.

    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> t1 = Triplet(Split({1}, {2, 3}))
    >>> t2 = Triplet(Split({2}, {1, 3}))
    >>> t3 = Triplet(Split({4}, {5, 6}))

    >>> # From triplets (grouped into profiles, equal weights per profile)
    >>> profileset = TripletProfileSet(profiles=[t1, t2, t3])
    >>> len(profileset)
    2
    >>> profileset.get_profile_weight(frozenset({1, 2, 3}))
    1.0

    >>> # From TripletProfile objects (better control)
    >>> profile1 = TripletProfile({t1: 0.8, t2: 0.2})
    >>> profile2 = TripletProfile([t3])
    >>> profileset2 = TripletProfileSet(profiles=[(profile1, 2.0), (profile2, 1.5)])
    >>> profileset2.get_profile_weight(frozenset({1, 2, 3}))
    2.0
    """

    def __init__(
        self,
        profiles: list[TripletProfile | Triplet | tuple[TripletProfile, float]] | None = None,
        taxa: frozenset[str] | None = None,
    ) -> None:
        """
        Initialize a triplet profile set.

        Parameters
        ----------
        profiles : list[TripletProfile | Triplet | tuple[TripletProfile, float]] | None, optional
            List of TripletProfile objects, Triplet objects, or tuples with profile weights.

            - If TripletProfile: used directly (optionally with profile-weight tuple)
            - If Triplet: automatically grouped by taxa into profiles. Each resulting
              profile receives default weight 1.0 and uses :class:`TripletProfile`
              to handle internal triplet weights (equal weights 1/k when created
              from a list of triplets).
            - If tuple: (profile, weight) where weight is the profile weight and
              must be positive.

            Passing triplets with explicit weights (e.g. ``(Triplet, weight)``) is
            not supported. To specify non-uniform triplet weights, construct a
            :class:`TripletProfile` explicitly and pass that instead.

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
        profiles_dict: dict[frozenset[str], tuple[TripletProfile, float]] = {}

        # Process items in a single pass, determining mode from first item
        mode: str | None = None  # 'profile' or 'triplet'
        profile_data: dict[frozenset[str], dict[Triplet, float]] = {}

        for item in profiles:
            # Extract object and weight (weights only allowed for TripletProfile)
            if isinstance(item, tuple):
                obj, weight = item
                if isinstance(obj, Triplet):
                    raise PhyloZooValueError(
                        "Triplet weights are not supported in TripletProfileSet. "
                        "Construct a TripletProfile with triplet weights and pass "
                        "that instead."
                    )
            else:
                obj = item
                weight = 1.0

            # Determine mode from first item, validate consistency for subsequent items
            if mode is None:
                if isinstance(obj, TripletProfile):
                    mode = 'profile'
                elif isinstance(obj, Triplet):
                    mode = 'triplet'
                else:
                    raise PhyloZooValueError(f"Expected TripletProfile or Triplet, got {type(obj)}")
            else:
                # Validate that subsequent items match the mode
                if mode == 'profile' and not isinstance(obj, TripletProfile):
                    if isinstance(obj, Triplet):
                        raise PhyloZooValueError("Cannot mix TripletProfile and Triplet objects in profiles list")
                    raise PhyloZooValueError(f"Expected TripletProfile, got {type(obj)}")
                elif mode == 'triplet' and not isinstance(obj, Triplet):
                    if isinstance(obj, TripletProfile):
                        raise PhyloZooValueError("Cannot mix TripletProfile and Triplet objects in profiles list")
                    raise PhyloZooValueError(f"Expected Triplet, got {type(obj)}")

            # Validate weight
            if weight <= 0:
                obj_type = "profile" if mode == 'profile' else "triplet"
                raise PhyloZooValueError(
                    f"{obj_type.capitalize()} weight must be positive, got {weight} for {obj_type} {obj}"
                )

            if mode == 'profile':
                # Mode 1: From TripletProfile objects
                profile = obj
                profile_taxa = profile.taxa
                all_taxa_from_input.update(profile_taxa)
                # Check for duplicate taxa sets
                if profile_taxa in profiles_dict:
                    raise PhyloZooValueError(
                        f"Multiple profiles with the same taxa set {profile_taxa} are not allowed. "
                        "Each 3-taxon set can only have one profile."
                    )
                profiles_dict[profile_taxa] = (profile, weight)

            else:
                # Mode 2: From Triplet objects (group into profiles, no triplet weights)
                triplet = obj
                triplet_taxa = triplet.taxa
                all_taxa_from_input.update(triplet_taxa)

                if triplet_taxa not in profile_data:
                    profile_data[triplet_taxa] = {}
                # Check for duplicate triplets during iteration
                if triplet in profile_data[triplet_taxa]:
                    raise PhyloZooValueError(
                        f"Triplet {triplet} appears multiple times in the input. "
                        "Each triplet can only appear once per taxa set."
                    )
                # Value is unused beyond duplicate detection; all triplets are unweighted here.
                profile_data[triplet_taxa][triplet] = 1.0

        # If we processed triplets, create profiles from grouped data.
        # Each taxa set becomes a TripletProfile built from its triplets, with
        # default profile weight 1.0 in the set.
        if mode == 'triplet':
            for taxa_set, triplets_dict in profile_data.items():
                # Validate no empty profiles
                if len(triplets_dict) == 0:
                    raise PhyloZooValueError(f"Cannot have empty profile for taxa {taxa_set}")

                profile = TripletProfile(list(triplets_dict.keys()))
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
                    f"Provided taxa must be a superset of all taxa in profiles/triplets. "
                    f"Missing taxa: {missing}"
                )
            self._taxa = frozenset(taxa)
        else:
            # Compute all taxa from profiles
            self._taxa = frozenset(all_taxa_from_input)

    @property
    def profiles(self) -> Mapping[frozenset[str], tuple[TripletProfile, float]]:
        """
        Get the profiles and their weights (read-only).

        Returns
        -------
        Mapping[frozenset[str], tuple[TripletProfile, float]]
            Read-only mapping of 3-taxon sets to (profile, profile_weight) tuples.
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
        Check if the triplet profile set is dense.

        A dense triplet profile set has a profile for every possible combination
        of 3 taxa from the total set of taxa.

        Returns
        -------
        bool
            True if the profile set is dense (has C(n, 3) profiles where n is
            the number of taxa), False otherwise.
        """
        n = len(self._taxa)
        if n < 3:
            return len(self._profiles) == 0
        return len(self._profiles) == comb(n, 3)

    @cached_property
    def is_all_resolved(self) -> bool:
        """
        Check if all profiles in the set are resolved.

        A profile is resolved if all its triplets are resolved (i.e., not star trees).
        This property returns True only if every profile in the set is resolved.

        Returns
        -------
        bool
            True if all profiles are resolved, False otherwise.

        Examples
        --------
        >>> from phylozoo.core.split.base import Split
        >>> from phylozoo.core.triplet.base import Triplet
        >>> t1 = Triplet(Split({1}, {2, 3}))
        >>> t2 = Triplet(Split({4}, {5, 6}))
        >>> profileset = TripletProfileSet(profiles=[t1, t2])
        >>> profileset.is_all_resolved
        True
        >>> star = Triplet({1, 2, 3})
        >>> profileset2 = TripletProfileSet(profiles=[t2, star])
        >>> profileset2.is_all_resolved
        False
        """
        return all(profile.is_resolved() for profile, _ in self._profiles.values())

    @cached_property
    def max_profile_len(self) -> int:
        """
        Get the maximum number of triplets in any profile.

        Returns
        -------
        int
            The largest number of triplets in any profile in the set.
            Returns 0 if the set is empty.

        Examples
        --------
        >>> from phylozoo.core.split.base import Split
        >>> from phylozoo.core.triplet.base import Triplet
        >>> t1 = Triplet(Split({1}, {2, 3}))
        >>> t2 = Triplet(Split({2}, {1, 3}))
        >>> t3 = Triplet(Split({3}, {1, 2}))
        >>> profileset = TripletProfileSet(profiles=[t1, t2, t3])
        >>> profileset.max_profile_len
        3
        """
        if len(self._profiles) == 0:
            return 0
        return max(len(profile) for profile, _ in self._profiles.values())

    def get_profile(self, taxa: frozenset[str]) -> TripletProfile | None:
        """
        Get the profile for a 3-taxon set.

        Parameters
        ----------
        taxa : frozenset[str]
            The 3-taxon set.

        Returns
        -------
        TripletProfile | None
            The profile for the taxa, or None if not found.
        """
        result = self._profiles.get(taxa)
        return result[0] if result else None

    def get_profile_weight(self, taxa: frozenset[str]) -> float | None:
        """
        Get the profile weight for a 3-taxon set.

        Parameters
        ----------
        taxa : frozenset[str]
            The 3-taxon set.

        Returns
        -------
        float | None
            The profile weight, or None if not found.
        """
        result = self._profiles.get(taxa)
        return result[1] if result else None

    def has_profile(self, taxa: frozenset[str]) -> bool:
        """
        Check if a profile exists for the given 3-taxon set.

        Parameters
        ----------
        taxa : frozenset[str]
            The 3-taxon set.

        Returns
        -------
        bool
            True if a profile exists, False otherwise.
        """
        return taxa in self._profiles

    def all_profile_taxon_sets(self) -> Iterator[frozenset[str]]:
        """
        Get all 3-taxon sets that have profiles.

        Returns
        -------
        Iterator[frozenset[str]]
            Iterator over all 3-taxon sets that have profiles.
        """
        return iter(self._profiles.keys())

    def __iter__(self) -> Iterator[tuple[TripletProfile, float]]:
        """
        Return an iterator over (profile, profile_weight) pairs.

        Returns
        -------
        Iterator[tuple[TripletProfile, float]]
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
        Check if a profile exists for the given 3-taxon set.

        Parameters
        ----------
        taxa : frozenset[str]
            The 3-taxon set.

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
            return "TripletProfileSet(profiles={})"

        # Build list of (profile, weight) tuples for initialization
        # Use the triplets dict directly since that's what __init__ accepts
        profile_items = []
        for taxa, (profile, weight) in self._profiles.items():
            # Get the triplets dict from the profile
            triplets_dict = dict(profile.triplets)
            profile_items.append(f"(TripletProfile({repr(triplets_dict)}), {weight})")

        profiles_str = ", ".join(profile_items)
        return f"TripletProfileSet(profiles=[{profiles_str}])"

    def __str__(self) -> str:
        """
        Return human-readable string representation of the triplet profile set.

        Displays one line per profile, showing the profile (using its __str__ method)
        and its profile weight. Aligns with TripletProfile's __str__ format.

        Returns
        -------
        str
            Human-readable string representation.

        Examples
        --------
        >>> from phylozoo.core.split.base import Split
        >>> t1 = Triplet(Split({1}, {2, 3}))
        >>> t2 = Triplet(Split({2}, {1, 3}))
        >>> profile1 = TripletProfile({t1: 0.8, t2: 0.2})
        >>> profileset = TripletProfileSet(profiles=[(profile1, 1.0)])
        >>> str(profileset)
        'TripletProfileSet({\\n  TripletProfile({...}) [weight: 1.0]\\n})'
        """
        if len(self._profiles) == 0:
            return "TripletProfileSet({})"

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

        return f"TripletProfileSet({{\n" + "\n".join(profile_lines) + "\n})"
