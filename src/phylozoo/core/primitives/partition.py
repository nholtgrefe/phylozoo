"""
Partition module.

This module provides a class for working with partitions of sets.
"""

import itertools
import warnings
from typing import Any, Iterator, TypeVar

T = TypeVar('T')


class Partition:
    """
    General class for partitions of sets.
    
    A partition is a collection of disjoint sets (called parts or blocks) whose
    union equals the original set. Each element appears in exactly one part.
    
    This class is immutable - once created, the partition structure cannot be
    modified. All parts are stored as frozensets.
    
    Parameters
    ----------
    parts : list[set[T]]
        List of sets of elements. The sets must be disjoint (no overlapping
        elements).
    
    Attributes
    ----------
    parts : tuple[frozenset, ...]
        Tuple of frozen sets representing the partition blocks in canonical form (read-only).
    elements : frozenset
        Frozen set containing the union of all elements in the partition (read-only).
    
    Raises
    ------
    ValueError
        If the sets overlap (i.e., the partition is invalid).
    
    Examples
    --------
    >>> partition = Partition([{1, 2}, {3, 4}, {5}])
    >>> len(partition)
    3
    >>> partition.size()
    5
    >>> {1, 2} in partition
    True
    >>> partition.get_part(3)
    frozenset({3, 4})
    >>> partition.parts  # Read-only tuple
    (frozenset({1, 2}), frozenset({3, 4}), frozenset({5}))
    
    Notes
    -----
    The partition is immutable after initialization. Attempts to modify attributes
    will raise AttributeError.
    """
    
    __slots__ = ('_parts', '_elements', '_initialized')
    
    def __init__(self, parts: list[set[T]]) -> None:
        # Convert to frozensets, compute elements, and validate in one pass
        parts_frozen: list[frozenset] = []
        elements_set: set[T] = set()
        total_size: int = 0
        
        for part in parts:
            part_frozen = frozenset(part)
            
            # Warn if empty set is added
            if len(part_frozen) == 0:
                warnings.warn(
                    "Empty set added to partition. This may cause unexpected behavior.",
                    UserWarning,
                    stacklevel=2
                )
            
            parts_frozen.append(part_frozen)
            part_size = len(part_frozen)
            total_size += part_size
            
            # Check for overlaps during element collection (early validation)
            for elt in part_frozen:
                if elt in elements_set:
                    raise ValueError("Invalid partition: sets overlap")
                elements_set.add(elt)
        
        # Validate total size matches (catches any remaining edge cases)
        if total_size != len(elements_set):
            raise ValueError("Invalid partition: sets overlap")
        
        # Store in canonical form
        self._parts: tuple = self._canonical_form(tuple(parts_frozen))
        self._elements: frozenset = frozenset(elements_set)
        self._initialized: bool = True
    
    @staticmethod
    def _canonical_form(parts: tuple) -> tuple:
        """
        Compute canonical form of a partition.
        
        The canonical form sorts parts by size first, then by their elements
        in sorted order. This ensures deterministic ordering regardless of
        input order.
        
        Parameters
        ----------
        parts : tuple
            Tuple of frozensets representing the partition parts.
        
        Returns
        -------
        tuple
            Tuple of frozensets in canonical (sorted) order.
        """
        def sort_key(part: frozenset) -> tuple:
            """Sort key: (size, sorted_elements). Uses simple sorting when possible."""
            try:
                # Try simple sorting first (works for same/comparable types)
                sorted_elements = sorted(part)
                return (len(part), tuple(sorted_elements))
            except TypeError:
                # Fall back to complex sorting for mixed types
                element_keys = []
                for elt in part:
                    elt_type = type(elt).__name__
                    try:
                        element_keys.append((elt_type, elt))
                    except TypeError:
                        element_keys.append((elt_type, str(elt)))
                element_keys.sort()
                return (len(part), tuple(element_keys))
        
        return tuple(sorted(parts, key=sort_key))
    
    def __setattr__(self, name: str, value: Any) -> None:
        """
        Prevent modification of attributes after initialization.
        
        Raises
        ------
        AttributeError
            If attempting to modify any attribute after initialization.
        """
        # Allow setting during initialization
        if not hasattr(self, '_initialized') or not self._initialized:
            super().__setattr__(name, value)
            return
        
        # Prevent modification after initialization
        if name in self.__slots__:
            raise AttributeError(
                f"Cannot modify attribute '{name}'. Partition is immutable."
            )
        raise AttributeError(
            f"Cannot set attribute '{name}'. Partition is immutable."
        )
    
    @property
    def parts(self) -> tuple:
        """
        Get the parts of the partition (read-only).
        
        Returns
        -------
        tuple
            Tuple of frozensets representing the partition blocks.
        """
        return self._parts
    
    @property
    def elements(self) -> frozenset:
        """
        Get the elements of the partition (read-only).
        
        Returns
        -------
        frozenset
            Frozen set containing all elements in the partition.
        """
        return self._elements
    
    def __repr__(self) -> str:
        """
        Return string representation of the partition.
        
        Returns
        -------
        str
            String representation showing the parts as sets.
        """
        unfrozen_set = [set(part) for part in self._parts]
        return f"Partition({unfrozen_set})"
    
    def __eq__(self, other: Any) -> bool:
        """
        Check if two partitions are equal.
        
        Two partitions are equal if they have the same parts (order doesn't matter).
        Since parts are stored in canonical (sorted) order, we can compare tuples directly.
        
        Parameters
        ----------
        other : Any
            Object to compare with.
        
        Returns
        -------
        bool
            True if partitions are equal, False otherwise.
        
        Examples
        --------
        >>> p1 = Partition([{1, 2}, {3, 4}])
        >>> p2 = Partition([{3, 4}, {1, 2}])
        >>> p1 == p2
        True
        """
        if not isinstance(other, Partition):
            return False
        # Since parts are stored in canonical order, we can compare directly
        return self._parts == other._parts
    
    def __hash__(self) -> int:
        """
        Return hash of the partition.
        
        Returns
        -------
        int
            Hash value based on the parts.
        
        Notes
        -----
        Partitions are hashable because parts are stored as frozensets.
        Since parts are stored in canonical (sorted) order, we can hash
        the tuple directly.
        """
        return hash(self._parts)
    
    def __contains__(self, subset: set[T, frozenset]) -> bool:
        """
        Check if a subset is one of the parts in the partition.
        
        Parameters
        ----------
        subset : set[T, frozenset]
            Subset to check.
        
        Returns
        -------
        bool
            True if the subset is a part of the partition, False otherwise.
        
        Examples
        --------
        >>> partition = Partition([{1, 2}, {3, 4}])
        >>> {1, 2} in partition
        True
        >>> {1, 3} in partition
        False
        """
        subset_frozen = frozenset(subset) if isinstance(subset, set) else subset
        return any(subset_frozen == part for part in self._parts)
    
    def __iter__(self) -> Iterator[frozenset]:
        """
        Return an iterator over the parts of the partition.
        
        The iteration order is deterministic and consistent regardless of
        the input order. Parts are sorted by size first, then by their
        elements in sorted order.
        
        Returns
        -------
        Iterator[frozenset]
            Iterator over the parts in sorted order.
        
        Examples
        --------
        >>> partition = Partition([{3, 4}, {1, 2}])
        >>> list(partition)
        [frozenset({1, 2}), frozenset({3, 4})]  # Always sorted, regardless of input order
        """
        # Parts are already stored in sorted order, so just iterate
        return iter(self._parts)
    
    def __len__(self) -> int:
        """
        Return the number of parts in the partition.
        
        Returns
        -------
        int
            Number of parts (blocks) in the partition.
        
        Examples
        --------
        >>> partition = Partition([{1, 2}, {3, 4}, {5}])
        >>> len(partition)
        3
        """
        return len(self._parts)
    
    def size(self) -> int:
        """
        Return the total number of elements the partition covers.
        
        Returns
        -------
        int
            Total number of elements in all parts.
        
        Examples
        --------
        >>> partition = Partition([{1, 2}, {3, 4}, {5}])
        >>> partition.size()
        5
        """
        return len(self._elements)
    
    def _is_valid(self) -> bool:
        """
        Check if the partition is valid.
        
        A partition is valid if each element appears in exactly one part
        of the partition (i.e., parts are disjoint and cover all elements).
        
        Returns
        -------
        bool
            True if partition is valid, False otherwise.
        
        Notes
        -----
        The implementation checks if the sum of sizes of all parts equals
        the size of the union of all parts, which is equivalent to ensuring
        no overlapping elements.
        """
        return sum(len(part) for part in self._parts) == len(self._elements)
    
    def get_part(self, element: T) -> frozenset:
        """
        Get the part that contains the given element.
        
        Parameters
        ----------
        element : T
            Element to find.
        
        Returns
        -------
        frozenset
            The part (frozenset) containing the element.
        
        Raises
        ------
        ValueError
            If the element is not found in any part of the partition.
        
        Examples
        --------
        >>> partition = Partition([{1, 2}, {3, 4}])
        >>> partition.get_part(1)
        frozenset({1, 2})
        >>> partition.get_part(5)
        ValueError: Element 5 not found in partition
        """
        for part in self._parts:
            if element in part:
                return part
        raise ValueError(f"Element {element} not found in partition")
    
    def subpartitions(self, size: int = 4) -> Iterator['Partition']:
        """
        Generate all subpartitions of a specified size.
        
        A subpartition is a partition formed by selecting a subset of the
        parts from the current partition.
        
        Parameters
        ----------
        size : int, optional
            Number of parts to include in each subpartition, by default 4.
        
        Yields
        ------
        Partition
            Subpartitions of the specified size.
        
        Examples
        --------
        >>> partition = Partition([{1}, {2}, {3}, {4}, {5}])
        >>> subparts = list(partition.subpartitions(size=2))
        >>> len(subparts)
        10  # C(5,2) = 10
        """
        for comb in itertools.combinations(self._parts, size):
            yield Partition([set(part) for part in comb])
    
    def representative_partitions(self) -> Iterator['Partition']:
        """
        Generate all partitions with exactly one element per part.
        
        For each part in the current partition, selects exactly one element,
        creating a new partition where each part is a singleton set.
        
        Yields
        ------
        Partition
            Representative partitions with exactly one element per part.
        
        Examples
        --------
        >>> partition = Partition([{1, 2}, {3, 4}])
        >>> reps = list(partition.representative_partitions())
        >>> len(reps)
        4  # 2 choices from first part * 2 choices from second part
        >>> reps[0]
        Partition([{1}, {3}])
        """
        for comb in itertools.product(*self._parts):
            yield Partition([{elt} for elt in comb])
    
    def is_refinement(self, other: 'Partition') -> bool:
        """
        Check if this partition is a refinement of another partition.
        
        A partition P is a refinement of partition Q if every part of P is
        a subset of some part of Q.
        
        Parameters
        ----------
        other : Partition
            The partition to check against.
        
        Returns
        -------
        bool
            True if this partition is a refinement of 'other', False otherwise.
        
        Raises
        ------
        ValueError
            If 'other' is not a Partition instance or covers different elements.
        
        Examples
        --------
        >>> p1 = Partition([{1}, {2}, {3}])  # Fine partition
        >>> p2 = Partition([{1, 2}, {3}])    # Coarser partition
        >>> p1.is_refinement(p2)
        True
        >>> p2.is_refinement(p1)
        False
        """
        if not isinstance(other, Partition):
            raise ValueError("The argument must be an instance of Partition")
        if not self._elements == other._elements:
            raise ValueError("Other partition covers different elements.")
        
        for part in self._parts:
            if not any(part.issubset(other_part) for other_part in other._parts):
                return False
        return True
    