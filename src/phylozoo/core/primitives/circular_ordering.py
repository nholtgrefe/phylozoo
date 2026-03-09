"""
Circular ordering module.

A circular ordering represents a cyclic arrangement of elements, such as the order
of taxa around a circular phylogenetic tree or network. Circular orderings are
used in network reconstruction algorithms and visualization. This module provides
classes for working with circular orderings of sets and elements.
"""

import itertools
import warnings
from typing import Any, Dict, Iterator, List, Set, TypeVar, Generator

from .partition import Partition
from phylozoo.utils.exceptions import PhyloZooWarning, PhyloZooValueError, PhyloZooAttributeError

T = TypeVar('T')


class CircularSetOrdering(Partition):
    """
    A circular ordering of sets that forms a partition.
    
    This class represents a partition with an additional circular ordering structure.
    The sets must be disjoint (as required by Partition), and the ordering is
    considered equivalent under cyclic permutations and reversals.
    
    The ordering is stored in canonical form (lexicographically smallest rotation,
    considering both forward and reverse directions), which makes equality and
    hashing operations efficient.
    
    Parameters
    ----------
    setorder : list[set[T]]
        List of sets in circular order.         The sets must be disjoint.
    
    Notes
    -----
    The circular ordering is immutable and stored in canonical form. Two
    CircularSetOrderings are considered equal if they are the same up to cyclic
    permutation or reversal. Since both are canonical, equality is a simple
    tuple comparison.
    
    Examples
    --------
    >>> cso = CircularSetOrdering([{1, 2}, {3}, {4, 5}])
    >>> len(cso)
    3
    >>> {1, 2} in cso
    True
    >>> cso.are_neighbors({1, 2}, {3})
    True
    """
    
    __slots__ = ('_setorder',)
    
    def __setattr__(self, name: str, value: Any) -> None:
        """
        Override to allow setting _setorder during initialization.
        
        Parameters
        ----------
        name : str
            Attribute name.
        value : Any
            Attribute value.

        Raises
        ------
        PhyloZooAttributeError
            If attempting to modify any attribute after initialization.
        PhyloZooValueError
            If the sets overlap.
        PhyloZooWarning
            If an empty set is added to the circular set ordering.
        """
        # Allow setting _setorder during initialization (before _initialized exists)
        if name == '_setorder' and not hasattr(self, '_initialized'):
            object.__setattr__(self, name, value)
            return
        
        # Use parent's __setattr__ for everything else
        super().__setattr__(name, value)
    
    def __init__(self, setorder: list[set[T]]) -> None:
        # Convert to frozensets, compute elements, and validate in one pass
        # (using optimized Partition approach)
        parts_frozen: list[frozenset] = []
        elements_set: set[T] = set()
        total_size: int = 0
        
        for part in setorder:
            part_frozen = frozenset(part)
            
            # Warn if empty set is added
            if len(part_frozen) == 0:
                warnings.warn(
                    "Empty set added to circular set ordering. This may cause unexpected behavior.",
                    PhyloZooWarning,
                    stacklevel=2
                )
            
            parts_frozen.append(part_frozen)
            part_size = len(part_frozen)
            total_size += part_size
            
            # Check for overlaps during element collection (early validation)
            for elt in part_frozen:
                if elt in elements_set:
                    raise PhyloZooValueError("Invalid partition: sets overlap")
                elements_set.add(elt)
        
        # Validate total size matches
        if total_size != len(elements_set):
            raise PhyloZooValueError("Invalid partition: sets overlap")
        
        # Store parts in original order (not sorted, unlike Partition)
        # We need to preserve order for circular ordering
        self._initialized = False
        object.__setattr__(self, '_parts', tuple(parts_frozen))
        object.__setattr__(self, '_elements', frozenset(elements_set))
        object.__setattr__(self, '_initialized', True)
        
        # Compute and store canonical ordering (lexicographically smallest rotation)
        canonical = self._canonical_form(tuple(parts_frozen))
        object.__setattr__(self, '_setorder', canonical)
    
    @property
    def setorder(self) -> tuple:
        """
        Get the circular ordering of sets (read-only).
        
        Returns
        -------
        tuple
            Tuple of frozensets representing the circular ordering.
        """
        return self._setorder
    
    def __repr__(self) -> str:
        """
        Return string representation of the circular set ordering.
        
        Returns
        -------
        str
            String representation showing the sets in order.
        """
        unfrozen_sets = [set(s) for s in self._setorder]
        return f"CircularSetOrdering({unfrozen_sets})"
    
    def __iter__(self) -> Iterator[frozenset]:
        """
        Return an iterator over the sets in circular order.
        
        Returns
        -------
        Iterator[frozenset]
            Iterator over the sets in the circular ordering.
        """
        return iter(self._setorder)
    
    def _canonical_form(self, ordering: tuple) -> tuple:
        """
        Compute canonical form of a circular ordering.
        
        The canonical form is the lexicographically smallest rotation,
        considering both forward and reverse directions.
        
        Uses an efficient O(n) algorithm: finds the minimum entry, then
        determines direction by comparing its two neighbors.
        
        Parameters
        ----------
        ordering : tuple
            The ordering to canonicalize.
        
        Returns
        -------
        tuple
            The canonical form of the ordering.
        """
        if not ordering:
            return ordering
        
        n = len(ordering)
        if n == 1:
            return ordering
        
        # Convert to comparable form for lexicographic comparison
        def _to_comparable_form(elt) -> tuple:
            """Convert element to comparable form for sorting."""
            if isinstance(elt, frozenset):
                # For frozensets, sort elements (handle mixed types like Partition does)
                try:
                    sorted_elements = sorted(elt)
                    return tuple(sorted_elements)
                except TypeError:
                    # Fall back to complex sorting for mixed types
                    element_keys = []
                    for item in elt:
                        item_type = type(item).__name__
                        try:
                            element_keys.append((item_type, item))
                        except TypeError:
                            element_keys.append((item_type, str(item)))
                    element_keys.sort()
                    return tuple(element_keys)
            # For non-frozensets (shouldn't happen in CircularSetOrdering, but handle gracefully)
            return elt
        
        # Find the lexicographically smallest entry
        min_idx = 0
        min_comp = _to_comparable_form(ordering[0])
        
        for i in range(1, n):
            comp = _to_comparable_form(ordering[i])
            if comp < min_comp:
                min_comp = comp
                min_idx = i
        
        # Get the two neighbors of the minimum entry (circular)
        left_neighbor_idx = (min_idx - 1) % n
        right_neighbor_idx = (min_idx + 1) % n
        
        left_neighbor_comp = _to_comparable_form(ordering[left_neighbor_idx])
        right_neighbor_comp = _to_comparable_form(ordering[right_neighbor_idx])
        
        # Determine direction: if right neighbor is smaller, we go forward;
        # if left neighbor is smaller, we go reverse
        if right_neighbor_comp < left_neighbor_comp:
            # Forward: start at min_idx, go right
            return ordering[min_idx:] + ordering[:min_idx]
        else:
            # Reverse: start at min_idx, go left (which is reverse)
            # Build reverse ordering starting from min_idx
            result = [ordering[min_idx]]
            for i in range(1, n):
                idx = (min_idx - i) % n
                result.append(ordering[idx])
            return tuple(result)
    
    def __eq__(self, other: Any) -> bool:
        """
        Check if two circular set orderings are equal.
        
        Two circular set orderings are equal if they have the same sets and
        the same circular ordering (up to cyclic permutation or reversal).
        Since both are stored in canonical form, we can compare directly.
        
        Parameters
        ----------
        other : Any
            Object to compare with.
        
        Returns
        -------
        bool
            True if circular set orderings are equal, False otherwise.
        
        Examples
        --------
        >>> cso1 = CircularSetOrdering([{1, 2}, {3}, {4}])
        >>> cso2 = CircularSetOrdering([{3}, {4}, {1, 2}])  # Cyclic permutation
        >>> cso1 == cso2
        True
        """
        if not isinstance(other, CircularSetOrdering):
            return False
        
        # Check partition equality (compare sets of parts, not tuples, since order differs)
        if set(self._parts) != set(other._parts):
            return False
        
        # Since both are in canonical form, direct comparison works
        return self._setorder == other._setorder
    
    def __hash__(self) -> int:
        """
        Return hash of the circular set ordering.
        
        Since the ordering is stored in canonical form, we can hash it directly.
        
        Returns
        -------
        int
            Hash value based on canonical circular ordering.
        """
        # Since _setorder is already in canonical form, hash it directly
        return hash(self._setorder)
    
    def are_singletons(self) -> bool:
        """
        Check if all sets in the ordering are singletons.
        
        Returns
        -------
        bool
            True if all sets have exactly one element, False otherwise.
        
        Examples
        --------
        >>> cso = CircularSetOrdering([{1}, {2}, {3}])
        >>> cso.are_singletons()
        True
        >>> cso2 = CircularSetOrdering([{1, 2}, {3}])
        >>> cso2.are_singletons()
        False
        """
        return all(len(s) == 1 for s in self._setorder)
    
    def to_circular_ordering(self) -> 'CircularOrdering':
        """
        Convert to CircularOrdering if all sets are singletons.
        
        Returns
        -------
        CircularOrdering
            A CircularOrdering with the elements from singleton sets.
        
        Raises
        ------
        PhyloZooValueError
            If not all sets are singletons.
        
        Examples
        --------
        >>> cso = CircularSetOrdering([{1}, {2}, {3}])
        >>> co = cso.to_circular_ordering()
        >>> co
        CircularOrdering([1, 2, 3])
        """
        if not self.are_singletons():
            raise PhyloZooValueError("Not all sets are singletons")
        elements = [next(iter(s)) for s in self._setorder]
        return CircularOrdering(elements)
    
    def are_neighbors(
        self, set1: set[T, frozenset], set2: set[T, frozenset]
    ) -> bool:
        """
        Check if two sets are neighbors in the circular ordering.
        
        Two sets are neighbors if they appear consecutively in the circular
        ordering (including the wrap-around from last to first).
        
        Parameters
        ----------
        set1 : set[T, frozenset]
            First set to check.
        set2 : set[T, frozenset]
            Second set to check.
        
        Returns
        -------
        bool
            True if sets are neighbors, False otherwise.
        
        Raises
        ------
        PhyloZooValueError
            If set1 equals set2, or if either set is not in the setorder.
        
        Examples
        --------
        >>> cso = CircularSetOrdering([{1, 2}, {3}, {4}])
        >>> cso.are_neighbors({1, 2}, {3})
        True
        >>> cso.are_neighbors({1, 2}, {4})
        True  # Wrap-around
        >>> cso.are_neighbors({1, 2}, {4})
        True
        """
        set1_frozen = frozenset(set1) if isinstance(set1, set) else set1
        set2_frozen = frozenset(set2) if isinstance(set2, set) else set2
        
        if set1_frozen == set2_frozen:
            raise PhyloZooValueError("Cannot check if a set is neighbour of itself")
        if set1_frozen not in self._setorder:
            raise PhyloZooValueError(f"Set {set1} not found in setorder")
        if set2_frozen not in self._setorder:
            raise PhyloZooValueError(f"Set {set2} not found in setorder")
        
        i = self._setorder.index(set1_frozen)
        j = self._setorder.index(set2_frozen)
        n = len(self._setorder)
        diff = abs(i - j)
        return diff == 1 or diff == n - 1
    
    def suborderings(self, size: int = 4) -> Iterator['CircularSetOrdering']:
        """
        Generate all suborderings of a specified size.
        
        A subordering is a circular set ordering formed by selecting a subset
        of sets from the current ordering (all combinations, not just consecutive).
        
        Parameters
        ----------
        size : int, optional
            Number of sets to include in each subordering, by default 4.
        
        Yields
        ------
        CircularSetOrdering
            Suborderings of the specified size.
        
        Examples
        --------
        >>> cso = CircularSetOrdering([{1}, {2}, {3}, {4}, {5}])
        >>> subs = list(cso.suborderings(size=2))
        >>> len(subs)
        10  # C(5,2) = 10 combinations
        """
        for comb in itertools.combinations(self._setorder, size):
            yield CircularSetOrdering([set(s) for s in comb])
    
    def representative_orderings(self) -> Iterator['CircularOrdering']:
        """
        Generate all CircularOrderings with exactly one element per set.
        
        For each set in the circular set ordering, selects exactly one element,
        creating all possible circular orderings.
        
        Yields
        ------
        CircularOrdering
            Representative circular orderings.
        
        Examples
        --------
        >>> cso = CircularSetOrdering([{1, 2}, {3}])
        >>> reps = list(cso.representative_orderings())
        >>> len(reps)
        2  # 2 choices from first set * 1 choice from second set
        
        Notes
        -----
        The number of representative orderings grows exponentially with the number
        of sets and their sizes. For large orderings, this method may be slow
        or consume significant memory. This method is a generator, so it yields
        results lazily rather than building a list in memory.
        """
        for combination in itertools.product(*self._setorder):
            yield CircularOrdering(list(combination))


class CircularOrdering(CircularSetOrdering):
    """
    A circular ordering of elements.
    
    This class represents a circular ordering where each element is in its own
    singleton set. It is a special case of CircularSetOrdering.
    
    The ordering is stored in canonical form (lexicographically smallest rotation,
    considering both forward and reverse directions), which makes equality and
    hashing operations efficient.
    
    Parameters
    ----------
    order : list[T]
        List of elements in circular order.
    
    Notes
    -----
    The circular ordering is immutable and stored in canonical form. Two
    CircularOrderings are considered equal if they are the same up to cyclic
    permutation or reversal. Since both are canonical, equality is a simple
    tuple comparison.
    
    Examples
    --------
    >>> co = CircularOrdering([1, 2, 3, 4])
    >>> len(co)
    4
    >>> 1 in co
    True
    >>> co.are_neighbors(1, 2)
    True
    """
    
    __slots__ = ('_order',)
    
    def __setattr__(self, name: str, value: Any) -> None:
        """
        Override to allow setting _order during initialization.
        
        Parameters
        ----------
        name : str
            Attribute name.
        value : Any
            Attribute value.

        Raises
        ------
        PhyloZooAttributeError
            If attempting to modify any attribute after initialization.
        PhyloZooValueError
            If the elements are not unique.
        """
        # Allow setting _order during initialization (before _initialized exists)
        if name == '_order' and not hasattr(self, '_initialized'):
            object.__setattr__(self, name, value)
            return
        
        # Use parent's __setattr__ for everything else
        super().__setattr__(name, value)
    
    def __init__(self, order: list[T]) -> None:
        # Check for duplicates before creating partition
        if len(order) != len(set(order)):
            raise PhyloZooValueError("Elements in CircularOrdering must be unique")
        
        # Convert elements to singleton sets for the partition
        parts = [{elt} for elt in order]
        
        # Initialize CircularSetOrdering first (computes canonical setorder)
        super().__init__(parts)
        
        # Extract canonical element order from canonical setorder
        # Since setorder is already in canonical form, extract elements in that order
        canonical_order = tuple(next(iter(s)) for s in self._setorder)
        object.__setattr__(self, '_order', canonical_order)
    
    @property
    def order(self) -> tuple:
        """
        Get the circular ordering of elements (read-only).
        
        Returns
        -------
        tuple
            Tuple of elements in the circular ordering.
        """
        return self._order
    
    def __repr__(self) -> str:
        """
        Return string representation of the circular ordering.
        
        Returns
        -------
        str
            String representation showing the elements in order.
        """
        return f"CircularOrdering({list(self._order)})"

    def __iter__(self) -> Iterator[T]:
        """
        Return an iterator over the elements in circular order.
        
        Returns
        -------
        Iterator[T]
            Iterator over the elements in the circular ordering.
        """
        return iter(self._order)
    
    def __contains__(self, element: T) -> bool:
        """
        Check if an element is in the circular ordering.
        
        Parameters
        ----------
        element : T
            Element to check.
        
        Returns
        -------
        bool
            True if element is in the ordering, False otherwise.
        
        Examples
        --------
        >>> co = CircularOrdering([1, 2, 3])
        >>> 1 in co
        True
        >>> 4 in co
        False
        """
        return element in self._order
    
    def __eq__(self, other: Any) -> bool:
        """
        Check if two circular orderings are equal.
        
        Two circular orderings are equal if they have the same elements and
        the same circular ordering (up to cyclic permutation or reversal).
        
        Parameters
        ----------
        other : Any
            Object to compare with.
        
        Returns
        -------
        bool
            True if circular orderings are equal, False otherwise.
        
        Examples
        --------
        >>> co1 = CircularOrdering([1, 2, 3])
        >>> co2 = CircularOrdering([2, 3, 1])  # Cyclic permutation
        >>> co1 == co2
        True
        """
        if not isinstance(other, CircularOrdering):
            return False
        return super().__eq__(other)
    
    def __hash__(self) -> int:
        """
        Return hash of the circular ordering.
        
        Returns
        -------
        int
            Hash value based on canonical circular ordering.
        """
        return super().__hash__()
    
    def are_neighbors(self, elt1: T, elt2: T) -> bool:
        """
        Check if two elements are neighbors in the circular ordering.
        
        Two elements are neighbors if they appear consecutively in the circular
        ordering (including the wrap-around from last to first).
        
        Parameters
        ----------
        elt1 : T
            First element to check.
        elt2 : T
            Second element to check.
        
        Returns
        -------
        bool
            True if elements are neighbors, False otherwise.
        
        Raises
        ------
        PhyloZooValueError
            If elt1 equals elt2, or if either element is not in the order.
        
        Examples
        --------
        >>> co = CircularOrdering([1, 2, 3, 4])
        >>> co.are_neighbors(1, 2)
        True
        >>> co.are_neighbors(1, 4)
        True  # Wrap-around
        """
        return super().are_neighbors({elt1}, {elt2})
    
    def suborderings(self, size: int = 4) -> Iterator['CircularOrdering']:
        """
        Generate all suborderings of a specified size.
        
        A subordering is a circular ordering formed by selecting a subset
        of elements from the current ordering (all combinations, not just consecutive).
        
        Parameters
        ----------
        size : int, optional
            Number of elements to include in each subordering, by default 4.
        
        Yields
        ------
        CircularOrdering
            Suborderings of the specified size.
        
        Examples
        --------
        >>> co = CircularOrdering([1, 2, 3, 4, 5])
        >>> subs = list(co.suborderings(size=2))
        >>> len(subs)
        10  # C(5,2) = 10 combinations
        """
        for comb in itertools.combinations(self._order, size):
            yield CircularOrdering(list(comb))
    
    def to_circular_setordering(
        self, mapping: dict[T, set[T]] | None = None
    ) -> CircularSetOrdering:
        """
        Convert to CircularSetOrdering with optional mapping.
        
        Maps each element to a set. If no mapping is provided, each element
        is mapped to its singleton set.
        
        Parameters
        ----------
        mapping : dict[T, set[T]] | None, optional
            Dictionary mapping elements to sets. If None, elements are mapped
            to singleton sets, by default None.
        
        Returns
        -------
        CircularSetOrdering
            A CircularSetOrdering with the mapped sets.
        
        Examples
        --------
        >>> co = CircularOrdering([1, 2, 3])
        >>> cso = co.to_circular_setordering()
        >>> cso
        CircularSetOrdering([{1}, {2}, {3}])
        >>> # With custom mapping
        >>> mapping = {1: {1, 10}, 2: {2}, 3: {3, 30}}
        >>> cso2 = co.to_circular_setordering(mapping)
        >>> cso2
        CircularSetOrdering([{1, 10}, {2}, {3, 30}])
        """
        if mapping is None:
            mapping = {elt: frozenset({elt}) for elt in self._order}
        
        new_setorder = [set(mapping[elt]) for elt in self._order]
        return CircularSetOrdering(new_setorder)
