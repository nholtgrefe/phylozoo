"""
Quartet profile module.

This module provides the QuartetProfile class for representing multiple quartets
on the same 4-taxon set with weights. A QuartetProfile always has total weight 1.0:
if no weights are provided, each quartet is assigned 1/k; if weights are provided,
they must sum to 1.0 (within a small tolerance) and are not scaled.
"""

from types import MappingProxyType
from typing import Iterator, Mapping

from ...utils.exceptions import PhyloZooValueError
from ..primitives.circular_ordering import CircularOrdering
from ..split.base import Split
from .base import Quartet

# Tolerance for checking that stored weights sum to 1.0.
_WEIGHT_SUM_TOLERANCE = 1e-9


class QuartetProfile:
    """
    Immutable profile for quartets on the same 4-taxon set.
    
    A QuartetProfile groups multiple quartets that share the same 4 taxa,
    each with an associated weight. The weights always sum to 1.0 (within a
    small tolerance), so the profile represents a probability distribution
    over quartet topologies.
    
    - If no weights are provided (list of quartets), each quartet is assigned
      equal weight 1/k, where k is the number of quartets.
    - If weights are provided (dict or list of (quartet, weight) tuples),
      they must sum to 1.0 (within tolerance); they are not scaled.
    
    Parameters
    ----------
    quartets : dict[Quartet, float] | Mapping[Quartet, float] | list[Quartet] | list[tuple[Quartet, float]]
        Input quartets. Can be:
        - A dictionary mapping quartets to weights (must sum to 1.0)
        - A list of quartets (each assigned weight 1/k)
        - A list of (quartet, weight) tuples (weights must sum to 1.0)
        Taxa are automatically extracted from the quartets.
    
    Attributes
    ----------
    taxa : frozenset[str]
        The 4 taxon labels (extracted from quartets).
    quartets : Mapping[Quartet, float]
        Read-only mapping of quartets to their weights.
    
    Raises
    ------
    PhyloZooValueError
        If quartets is empty, if quartets have different taxa, if any
        weight is non-positive, or if provided weights do not sum to 1.0.
    
    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> q1 = Quartet(Split({1, 2}, {3, 4}))
    >>> q2 = Quartet(Split({1, 3}, {2, 4}))
    >>> # From dictionary with weights (already sum to 1.0)
    >>> profile = QuartetProfile({q1: 0.8, q2: 0.2})
    >>> profile.taxa
    frozenset({1, 2, 3, 4})
    >>> profile.get_weight(q1)
    0.8
    >>> # From list of quartets (equal weight 1/k each)
    >>> profile2 = QuartetProfile([q1, q2])
    >>> profile2.get_weight(q1)
    0.5
    >>> profile2.get_weight(q2)
    0.5
    """
    
    __slots__ = ('_taxa', '_quartets', '_initialized', '_split_cache', '_circular_orderings_cache')
    
    def __init__(
        self,
        quartets: (
            dict[Quartet, float]
            | Mapping[Quartet, float]
            | list[Quartet]
            | list[tuple[Quartet, float]]
        ),
    ) -> None:
        """
        Initialize a quartet profile.
        
        Total weight is always 1.0:
        - List of quartets: each gets weight 1/k (k = number of quartets).
        - Dict or list of (quartet, weight) tuples: weights must sum to 1.0 (within
          tolerance). If they do not sum to 1.0, initialization
          is invalid.
        
        Parameters
        ----------
        quartets : dict[Quartet, float] | Mapping[Quartet, float] | list[Quartet] | list[tuple[Quartet, float]]
            Input quartets. Can be:
            - A dictionary mapping quartets to weights (must sum to 1.0)
            - A list of quartets (each assigned weight 1/k)
            - A list of (quartet, weight) tuples (weights must sum to 1.0)
            Taxa are automatically extracted from quartets.
        
        Raises
        ------
        PhyloZooValueError
            If quartets is empty, if quartets have different taxa, if any quartet appears
            multiple times, if any weight is non-positive, or if provided weights do not
            sum to 1.0 (within tolerance).
        """
        if isinstance(quartets, list):
            # Check if it's a list of quartets or list of tuples
            if len(quartets) == 0:
                raw_dict: dict[Quartet, float] = {}
            elif isinstance(quartets[0], Quartet):
                # List of quartets: assign equal weight 1/k to each
                raw_dict = {}
                for q in quartets:
                    if q in raw_dict:
                        raise PhyloZooValueError(
                            f"Quartet {q} appears multiple times in the input. "
                            "Each quartet can only appear once in a profile."
                        )
                    raw_dict[q] = 1.0
                k = len(raw_dict)
                raw_dict = {q: 1.0 / k for q in raw_dict}
            else:
                # List of tuples: convert to dict (validate after)
                raw_dict = {}
                for q, weight in quartets:
                    if q in raw_dict:
                        raise PhyloZooValueError(
                            f"Quartet {q} appears multiple times in the input. "
                            "Each quartet can only appear once in a profile."
                        )
                    raw_dict[q] = weight
        else:
            raw_dict = dict(quartets)
        
        # Set quartets before validation
        object.__setattr__(self, '_quartets', raw_dict)
        self._validate_quartets()
        
        # When weights are provided they must sum to 1.0 (no scaling). List-of-quartets
        # case already has 1/k so sum is 1.0.
        total = sum(self._quartets.values())
        if abs(total - 1.0) > _WEIGHT_SUM_TOLERANCE:
            raise PhyloZooValueError(
                f"Weights must sum to 1.0 (within tolerance {_WEIGHT_SUM_TOLERANCE}), got {total}"
            )
        
        # Extract taxa and store as immutable
        first_quartet = next(iter(self._quartets.keys()))
        taxa_set = first_quartet.taxa
        object.__setattr__(self, '_taxa', taxa_set)
        object.__setattr__(self, '_quartets', MappingProxyType(self._quartets))
        object.__setattr__(self, '_initialized', True)
    
    def _validate_quartets(self) -> None:
        """
        Validate quartets.
        
        Raises PhyloZooValueError if validation fails. Does nothing if validation passes.
        Uses self._quartets for validation.
        
        Raises
        ------
        PhyloZooValueError
            If quartets is empty, if quartets have different taxa, or if any
            weight is non-positive.
        """
        if len(self._quartets) == 0:
            raise PhyloZooValueError("QuartetProfile must have at least one quartet")
        
        # Get taxa from first quartet
        first_quartet = next(iter(self._quartets.keys()))
        taxa_set = first_quartet.taxa
        
        # Validate all quartets have same taxa
        for quartet in self._quartets:
            if quartet.taxa != taxa_set:
                raise PhyloZooValueError(
                    f"All quartets must have the same taxa. "
                    f"Expected {taxa_set}, got {quartet.taxa}"
                )
        
        # Validate weights are positive
        for quartet, weight in self._quartets.items():
            if weight <= 0:
                raise PhyloZooValueError(
                    f"Weight must be positive, got {weight} for quartet {quartet}"
                )
    
    def __setattr__(self, name: str, value: any) -> None:
        """
        Prevent modification of attributes after initialization.
        
        Raises
        ------
        AttributeError
            If attempting to modify any attribute after initialization.
        """
        if not hasattr(self, '_initialized') or not self._initialized:
            super().__setattr__(name, value)
            return
        
        raise AttributeError(
            f"Cannot modify attribute '{name}'. QuartetProfile is immutable."
        )
    
    @property
    def taxa(self) -> frozenset[str]:
        """
        Get the taxa of the profile.
        
        Returns
        -------
        frozenset[str]
            The 4 taxon labels.
        """
        return self._taxa
    
    @property
    def quartets(self) -> Mapping[Quartet, float]:
        """
        Get the quartets and their weights (read-only).
        
        Returns
        -------
        Mapping[Quartet, float]
            Read-only mapping of quartets to weights.
        """
        return self._quartets
    
    def get_weight(self, quartet: Quartet) -> float:
        """
        Get the weight of a quartet.
        
        Parameters
        ----------
        quartet : Quartet
            The quartet to get the weight for.
        
        Returns
        -------
        float
            The weight of the quartet, or 0.0 if not found.
        """
        return self._quartets.get(quartet, 0.0)
    
    @property
    def split(self) -> Split | None:
        """
        Get the split if the profile has a single quartet.
        
        Returns the split of the quartet if the profile contains exactly one
        quartet. If that quartet is a star tree, returns None. If the profile
        has multiple quartets, returns None (since there would be multiple splits).
        
        The result is cached after first computation.
        
        Returns
        -------
        Split | None
            The split of the single quartet, or None if multiple quartets or star tree.
        """
        if hasattr(self, '_split_cache'):
            return self._split_cache
        
        if len(self._quartets) == 1:
            # Single quartet: return its split (None if star tree)
            quartet = next(iter(self._quartets))
            result = quartet.split
        else:
            # Multiple quartets: return None
            result = None
        
        # Cache the result
        object.__setattr__(self, '_split_cache', result)
        return result
    
    @property
    def circular_orderings(self) -> frozenset[CircularOrdering] | None:
        """
        Get circular orderings that are congruent with every quartet in the profile.
        
        The result is the intersection of all circular orderings of the quartets.
        Since star trees are congruent with all orderings, only resolved quartets
        constrain the result.
        
        The result is cached after first computation.
        
        Returns
        -------
        frozenset[CircularOrdering] | None
            Set of circular orderings congruent with all quartets, or None if
            no such ordering exists (e.g., if all three possible resolved quartets
            are in the profile).
        """
        if hasattr(self, '_circular_orderings_cache'):
            return self._circular_orderings_cache
        
        # Consider only resolved quartets (star trees are congruent with all orderings)
        resolved_quartets = [q for q in self._quartets if q.is_resolved()]
        
        if len(resolved_quartets) == 0:
            # Only star trees: return all 3 circular orderings
            taxa_list = sorted(self._taxa)
            a, b, c, d = taxa_list
            result = frozenset([
                CircularOrdering([a, b, c, d]),
                CircularOrdering([a, c, b, d]),
                CircularOrdering([a, b, d, c]),
            ])
        elif len(resolved_quartets) == 1:
            # Single resolved quartet: return its circular orderings
            result = resolved_quartets[0].circular_orderings
        elif len(resolved_quartets) == 2:
            # Two resolved quartets: return intersection of their circular orderings
            orderings1 = resolved_quartets[0].circular_orderings
            orderings2 = resolved_quartets[1].circular_orderings
            intersection = orderings1 & orderings2
            result = intersection if len(intersection) > 0 else None
        else:
            # Three or more resolved quartets: no ordering is congruent with all
            result = None
        
        # Cache the result
        object.__setattr__(self, '_circular_orderings_cache', result)
        return result
    
    def __len__(self) -> int:
        """
        Return the number of quartets in the profile.
        
        Returns
        -------
        int
            Number of quartets.
        """
        return len(self._quartets)
    
    def is_trivial(self) -> bool:
        """
        Check if the profile is trivial (contains exactly one quartet).
        
        A trivial profile is essentially a single quartet, meaning it represents
        a single topology rather than a distribution over multiple topologies.
        
        Returns
        -------
        bool
            True if the profile contains exactly one quartet, False otherwise.
        
        Examples
        --------
        >>> from phylozoo.core.split.base import Split
        >>> q1 = Quartet(Split({1, 2}, {3, 4}))
        >>> profile = QuartetProfile([q1])
        >>> profile.is_trivial()
        True
        >>> q2 = Quartet(Split({1, 3}, {2, 4}))
        >>> profile2 = QuartetProfile({q1: 0.8, q2: 0.2})
        >>> profile2.is_trivial()
        False
        """
        return len(self._quartets) == 1
    
    def is_resolved(self) -> bool:
        """
        Check if the profile is resolved.
        
        A profile is resolved if all its quartets are resolved. A quartet is resolved
        if it has a non-trivial split (i.e., it's not a star tree).
        
        Returns
        -------
        bool
            True if all quartets are resolved (profile is resolved), False otherwise.
        
        Examples
        --------
        >>> from phylozoo.core.split.base import Split
        >>> q1 = Quartet(Split({1, 2}, {3, 4}))
        >>> q2 = Quartet(Split({1, 3}, {2, 4}))
        >>> profile = QuartetProfile([q1, q2])
        >>> profile.is_resolved()
        True
        >>> q3 = Quartet({1, 2, 3, 4})  # Star tree
        >>> profile2 = QuartetProfile([q1, q3])
        >>> profile2.is_resolved()
        False
        """
        return all(q.is_resolved() for q in self._quartets)
    
    def __iter__(self) -> Iterator[Quartet]:
        """
        Return an iterator over the quartets.
        
        Returns
        -------
        Iterator[Quartet]
            Iterator over quartets.
        """
        return iter(self._quartets)
    
    def __contains__(self, quartet: Quartet) -> bool:
        """
        Check if a quartet is in the profile.
        
        Parameters
        ----------
        quartet : Quartet
            Quartet to check.
        
        Returns
        -------
        bool
            True if the quartet is in the profile, False otherwise.
        """
        return quartet in self._quartets
    
    def __repr__(self) -> str:
        """
        Return string representation of the profile.
        
        Returns
        -------
        str
            String representation.
        """
        return f"QuartetProfile(taxa={set(self._taxa)}, quartets={dict(self._quartets)})"
    
    def __str__(self) -> str:
        """
        Return human-readable string representation of the quartet profile.
        
        Displays one line per quartet, showing the quartet (using its __str__ method)
        and its weight.
        
        Returns
        -------
        str
            Human-readable string representation.
        
        Examples
        --------
        >>> from phylozoo.core.split.base import Split
        >>> q1 = Quartet(Split({1, 2}, {3, 4}))
        >>> q2 = Quartet(Split({1, 3}, {2, 4}))
        >>> profile = QuartetProfile({q1: 0.8, q2: 0.2})
        >>> str(profile)
        'QuartetProfile({\\n  Quartet(1 2 | 3 4): 0.8,\\n  Quartet(1 3 | 2 4): 0.2\\n})'
        """
        if len(self._quartets) == 0:
            return "QuartetProfile({})"
        
        # Sort quartets for consistent display
        # Sort by quartet string representation for deterministic ordering
        sorted_quartets = sorted(self._quartets.items(), key=lambda item: str(item[0]))
        
        # Show all quartets with weights, one per line
        quartet_lines = [
            f"  {quartet}: {weight}," for quartet, weight in sorted_quartets
        ]
        # Remove trailing comma from last line
        if quartet_lines:
            quartet_lines[-1] = quartet_lines[-1].rstrip(',')
        
        return f"QuartetProfile({{\n" + "\n".join(quartet_lines) + "\n})"

