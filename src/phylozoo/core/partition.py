"""
Partition module.

This module provides a class for working with partitions of sets.
"""

import itertools
from typing import Any, Iterator, List, Set, TypeVar, Union

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
    parts : List[Set[T]]
        List of sets of elements. The sets must be disjoint (no overlapping
        elements).
    
    Attributes
    ----------
    parts : Tuple[frozenset, ...]
        Tuple of frozen sets representing the partition blocks (read-only).
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
    
    def __init__(self, parts: List[Set[T]]) -> None:
        # Store in private attributes
        self._parts: tuple = tuple(frozenset(part) for part in parts)
        self._elements: frozenset = (
            frozenset().union(*self._parts) if self._parts else frozenset()
        )
        self._initialized: bool = True
        
        # Always validate
        if not self._is_valid():
            raise ValueError("Invalid partition: sets overlap")
    
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
        return set(self.parts) == set(other.parts)
    
    def __hash__(self) -> int:
        """
        Return hash of the partition.
        
        Returns
        -------
        int
            Hash value based on the parts (order-independent).
        
        Notes
        -----
        Partitions are hashable because parts are stored as frozensets.
        The hash is order-independent to match __eq__ behavior.
        """
        return hash(frozenset(self._parts))
    
    def __contains__(self, subset: Union[Set[T], frozenset]) -> bool:
        """
        Check if a subset is one of the parts in the partition.
        
        Parameters
        ----------
        subset : Union[Set[T], frozenset]
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
        # Sort by size first, then by sorted elements for deterministic order
        sorted_parts = sorted(
            self._parts,
            key=lambda part: (len(part), sorted(part))
        )
        return iter(sorted_parts)
    
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
    
    def subpartitions(self, size: int = 4) -> List['Partition']:
        """
        Generate all subpartitions of a specified size.
        
        A subpartition is a partition formed by selecting a subset of the
        parts from the current partition.
        
        Parameters
        ----------
        size : int, optional
            Number of parts to include in each subpartition, by default 4.
        
        Returns
        -------
        List[Partition]
            List of all subpartitions of the specified size.
        
        Examples
        --------
        >>> partition = Partition([{1}, {2}, {3}, {4}, {5}])
        >>> subparts = partition.subpartitions(size=2)
        >>> len(subparts)
        10  # C(5,2) = 10
        """
        subparts = self._k_combinations(list(self._parts), size)
        return [Partition([set(part) for part in subpart]) for subpart in subparts]
    
    def representative_partitions(self) -> List['Partition']:
        """
        Generate all partitions with exactly one element per part.
        
        For each part in the current partition, selects exactly one element,
        creating a new partition where each part is a singleton set.
        
        Returns
        -------
        List[Partition]
            List of all representative partitions.
        
        Examples
        --------
        >>> partition = Partition([{1, 2}, {3, 4}])
        >>> reps = partition.representative_partitions()
        >>> len(reps)
        4  # 2 choices from first part * 2 choices from second part
        >>> reps[0]
        Partition([{1}, {3}])
        """
        combinations = list(itertools.product(*self._parts))
        partitions = [[{elt} for elt in comb] for comb in combinations]
        return [Partition(part) for part in partitions]
    
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
    
    @staticmethod
    def _k_combinations(structure: List[Any], k: int) -> List[List[Any]]:
        """
        Helper function: returns all combinations of size k from a given list.
        
        Parameters
        ----------
        structure : List[Any]
            List to generate combinations from.
        k : int
            Size of combinations to generate.
        
        Returns
        -------
        List[List[Any]]
            List of all combinations of size k.
        
        Examples
        --------
        >>> Partition._k_combinations([1, 2, 3], 2)
        [[1, 2], [1, 3], [2, 3]]
        """
        return [list(comb) for comb in itertools.combinations(structure, k)]
