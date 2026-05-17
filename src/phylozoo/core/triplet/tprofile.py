"""
Triplet profile module.

A triplet profile groups multiple triplets on the same 3-taxon set, each with an
associated weight representing the relative importance or frequency of each triplet
topology. This module provides the TripletProfile class; total weight is always 1.0.
"""

from types import MappingProxyType
from typing import Iterator, Mapping

from ...utils.exceptions import PhyloZooValueError
from ..split.base import Split
from .base import Triplet

# Tolerance for checking that stored weights sum to 1.0.
_WEIGHT_SUM_TOLERANCE = 1e-9


class TripletProfile:
    """
    Immutable profile for triplets on the same 3-taxon set.

    A TripletProfile groups multiple triplets that share the same 3 taxa,
    each with an associated weight. The weights always sum to 1.0 (within a
    small tolerance), so the profile represents a probability distribution
    over triplet topologies.

    - If no weights are provided (list of triplets), each triplet is assigned
      equal weight 1/k, where k is the number of triplets.
    - If weights are provided (dict or list of (triplet, weight) tuples),
      they must sum to 1.0 (within tolerance); they are not scaled.

    Parameters
    ----------
    triplets : dict[Triplet, float] | Mapping[Triplet, float] | list[Triplet] | list[tuple[Triplet, float]]
        Input triplets. Can be:

        - A dictionary mapping triplets to weights (must sum to 1.0)
        - A list of triplets (each assigned weight 1/k)
        - A list of (triplet, weight) tuples (weights must sum to 1.0)

        Taxa are automatically extracted from the triplets.

    Raises
    ------
    PhyloZooValueError
        If triplets is empty, if triplets have different taxa, if any
        weight is non-positive, or if provided weights do not sum to 1.0.

    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> t1 = Triplet(Split({1}, {2, 3}))
    >>> t2 = Triplet(Split({2}, {1, 3}))

    >>> # From dictionary with weights (already sum to 1.0)
    >>> profile = TripletProfile({t1: 0.8, t2: 0.2})
    >>> profile.taxa
    frozenset({1, 2, 3})
    >>> profile.get_weight(t1)
    0.8

    >>> # From list of triplets (equal weight 1/k each)
    >>> profile2 = TripletProfile([t1, t2])
    >>> profile2.get_weight(t1)
    0.5
    >>> profile2.get_weight(t2)
    0.5
    """

    __slots__ = ("_taxa", "_triplets", "_initialized", "_split_cache")

    # Slot type annotations (set via object.__setattr__ in __init__).
    _taxa: frozenset[str]
    _triplets: Mapping[Triplet, float]
    _initialized: bool
    _split_cache: Split | None

    def __init__(
        self,
        triplets: (
            dict[Triplet, float]
            | Mapping[Triplet, float]
            | list[Triplet]
            | list[tuple[Triplet, float]]
        ),
    ) -> None:
        """
        Initialize a triplet profile.

        Total weight is always 1.0:
        - List of triplets: each gets weight 1/k (k = number of triplets).
        - Dict or list of (triplet, weight) tuples: weights must sum to 1.0 (within
          tolerance). If they do not sum to 1.0, initialization
          is invalid.

        Parameters
        ----------
        triplets : dict[Triplet, float] | Mapping[Triplet, float] | list[Triplet] | list[tuple[Triplet, float]]
            Input triplets. Can be:

            - A dictionary mapping triplets to weights (must sum to 1.0)
            - A list of triplets (each assigned weight 1/k)
            - A list of (triplet, weight) tuples (weights must sum to 1.0)
            Taxa are automatically extracted from triplets.

        Raises
        ------
        PhyloZooValueError
            If triplets is empty, if triplets have different taxa, if any triplet appears
            multiple times, if any weight is non-positive, or if provided weights do not
            sum to 1.0 (within tolerance).
        """
        if isinstance(triplets, list):
            # Check if it's a list of triplets or list of tuples
            if len(triplets) == 0:
                raw_dict: dict[Triplet, float] = {}
            elif isinstance(triplets[0], Triplet):
                # List of triplets: assign equal weight 1/k to each
                raw_dict = {}
                for t in triplets:
                    if t in raw_dict:
                        raise PhyloZooValueError(
                            f"Triplet {t} appears multiple times in the input. "
                            "Each triplet can only appear once in a profile."
                        )
                    raw_dict[t] = 1.0
                k = len(raw_dict)
                raw_dict = {t: 1.0 / k for t in raw_dict}
            else:
                # List of tuples: convert to dict (validate after)
                raw_dict = {}
                for t, weight in triplets:
                    if t in raw_dict:
                        raise PhyloZooValueError(
                            f"Triplet {t} appears multiple times in the input. "
                            "Each triplet can only appear once in a profile."
                        )
                    raw_dict[t] = weight
        else:
            raw_dict = dict(triplets)

        # Set triplets before validation
        object.__setattr__(self, "_triplets", raw_dict)
        self._validate_triplets()

        # When weights are provided they must sum to 1.0 (no scaling). List-of-triplets
        # case already has 1/k so sum is 1.0.
        total = sum(self._triplets.values())
        if abs(total - 1.0) > _WEIGHT_SUM_TOLERANCE:
            raise PhyloZooValueError(
                f"Weights must sum to 1.0 (within tolerance {_WEIGHT_SUM_TOLERANCE}), got {total}"
            )

        # Extract taxa and store as immutable
        first_triplet = next(iter(self._triplets.keys()))
        taxa_set = first_triplet.taxa
        object.__setattr__(self, "_taxa", taxa_set)
        object.__setattr__(self, "_triplets", MappingProxyType(self._triplets))
        object.__setattr__(self, "_initialized", True)

    def _validate_triplets(self) -> None:
        """
        Validate triplets.

        Raises PhyloZooValueError if validation fails. Does nothing if validation passes.
        Uses self._triplets for validation.

        Raises
        ------
        PhyloZooValueError
            If triplets is empty, if triplets have different taxa, or if any
            weight is non-positive.
        """
        if len(self._triplets) == 0:
            raise PhyloZooValueError("TripletProfile must have at least one triplet")

        # Get taxa from first triplet
        first_triplet = next(iter(self._triplets.keys()))
        taxa_set = first_triplet.taxa

        # Validate all triplets have same taxa
        for triplet in self._triplets:
            if triplet.taxa != taxa_set:
                raise PhyloZooValueError(
                    f"All triplets must have the same taxa. "
                    f"Expected {taxa_set}, got {triplet.taxa}"
                )

        # Validate weights are positive
        for triplet, weight in self._triplets.items():
            if weight <= 0:
                raise PhyloZooValueError(
                    f"Weight must be positive, got {weight} for triplet {triplet}"
                )

    def __setattr__(self, name: str, value: object) -> None:
        """
        Prevent modification of attributes after initialization.

        Raises
        ------
        AttributeError
            If attempting to modify any attribute after initialization.
        """
        if not hasattr(self, "_initialized") or not self._initialized:
            super().__setattr__(name, value)
            return

        raise AttributeError(f"Cannot modify attribute '{name}'. TripletProfile is immutable.")

    @property
    def taxa(self) -> frozenset[str]:
        """
        Get the taxa of the profile.

        Returns
        -------
        frozenset[str]
            The 3 taxon labels.
        """
        return self._taxa

    @property
    def triplets(self) -> Mapping[Triplet, float]:
        """
        Get the triplets and their weights (read-only).

        Returns
        -------
        Mapping[Triplet, float]
            Read-only mapping of triplets to weights.
        """
        return self._triplets

    def get_weight(self, triplet: Triplet) -> float:
        """
        Get the weight of a triplet.

        Parameters
        ----------
        triplet : Triplet
            The triplet to get the weight for.

        Returns
        -------
        float
            The weight of the triplet, or 0.0 if not found.
        """
        return self._triplets.get(triplet, 0.0)

    @property
    def split(self) -> Split | None:
        """
        Get the split if the profile has a single triplet.

        Returns the split of the triplet if the profile contains exactly one
        triplet. If that triplet is a star tree, returns None. If the profile
        has multiple triplets, returns None (since there would be multiple splits).

        The result is cached after first computation.

        Returns
        -------
        Split | None
            The split of the single triplet, or None if multiple triplets or star tree.
        """
        if hasattr(self, "_split_cache"):
            return self._split_cache

        if len(self._triplets) == 1:
            # Single triplet: return its split (None if star tree)
            triplet = next(iter(self._triplets))
            result = triplet.split
        else:
            # Multiple triplets: return None
            result = None

        # Cache the result
        object.__setattr__(self, "_split_cache", result)
        return result

    def __len__(self) -> int:
        """
        Return the number of triplets in the profile.

        Returns
        -------
        int
            Number of triplets.
        """
        return len(self._triplets)

    def is_trivial(self) -> bool:
        """
        Check if the profile is trivial (contains exactly one triplet).

        A trivial profile is essentially a single triplet, meaning it represents
        a single topology rather than a distribution over multiple topologies.

        Returns
        -------
        bool
            True if the profile contains exactly one triplet, False otherwise.

        Examples
        --------
        >>> from phylozoo.core.split.base import Split
        >>> t1 = Triplet(Split({1}, {2, 3}))
        >>> profile = TripletProfile([t1])
        >>> profile.is_trivial()
        True
        >>> t2 = Triplet(Split({2}, {1, 3}))
        >>> profile2 = TripletProfile({t1: 0.8, t2: 0.2})
        >>> profile2.is_trivial()
        False
        """
        return len(self._triplets) == 1

    def is_resolved(self) -> bool:
        """
        Check if the profile is resolved.

        A profile is resolved if all its triplets are resolved. A triplet is resolved
        if it has a trivial split (i.e., it's not a star tree).

        Returns
        -------
        bool
            True if all triplets are resolved (profile is resolved), False otherwise.

        Examples
        --------
        >>> from phylozoo.core.split.base import Split
        >>> t1 = Triplet(Split({1}, {2, 3}))
        >>> t2 = Triplet(Split({2}, {1, 3}))
        >>> profile = TripletProfile([t1, t2])
        >>> profile.is_resolved()
        True
        >>> t3 = Triplet({1, 2, 3})  # Star tree
        >>> profile2 = TripletProfile([t1, t3])
        >>> profile2.is_resolved()
        False
        """
        return all(t.is_resolved() for t in self._triplets)

    def __iter__(self) -> Iterator[Triplet]:
        """
        Return an iterator over the triplets.

        Returns
        -------
        Iterator[Triplet]
            Iterator over triplets.
        """
        return iter(self._triplets)

    def __contains__(self, triplet: Triplet) -> bool:
        """
        Check if a triplet is in the profile.

        Parameters
        ----------
        triplet : Triplet
            Triplet to check.

        Returns
        -------
        bool
            True if the triplet is in the profile, False otherwise.
        """
        return triplet in self._triplets

    def __repr__(self) -> str:
        """
        Return string representation of the profile.

        Returns
        -------
        str
            String representation.
        """
        return f"TripletProfile(taxa={set(self._taxa)}, triplets={dict(self._triplets)})"

    def __str__(self) -> str:
        """
        Return human-readable string representation of the triplet profile.

        Displays one line per triplet, showing the triplet (using its __str__ method)
        and its weight.

        Returns
        -------
        str
            Human-readable string representation.

        Examples
        --------
        >>> from phylozoo.core.split.base import Split
        >>> t1 = Triplet(Split({1}, {2, 3}))
        >>> t2 = Triplet(Split({2}, {1, 3}))
        >>> profile = TripletProfile({t1: 0.8, t2: 0.2})
        >>> str(profile)
        'TripletProfile({\\n  Triplet(1 | 2 3): 0.8,\\n  Triplet(2 | 1 3): 0.2\\n})'
        """
        if len(self._triplets) == 0:
            return "TripletProfile({})"

        # Sort triplets for consistent display
        # Sort by triplet string representation for deterministic ordering
        sorted_triplets = sorted(self._triplets.items(), key=lambda item: str(item[0]))

        # Show all triplets with weights, one per line
        triplet_lines = [f"  {triplet}: {weight}," for triplet, weight in sorted_triplets]
        # Remove trailing comma from last line
        if triplet_lines:
            triplet_lines[-1] = triplet_lines[-1].rstrip(",")

        return "TripletProfile({\n" + "\n".join(triplet_lines) + "\n})"
