"""
Splits module.

This module provides classes for working with phylogenetic splits. A split is a
2-partition {A, B} of a set of elements where A ∪ B equals the full set and
A ∩ B = ∅.
"""

from functools import cached_property
from typing import Set, TypeVar

from ...utils.exceptions import PhyloZooValueError
from ..primitives.partition import Partition

T = TypeVar('T')


class Split(Partition):
    """
    Class for 2-partitions of sets, child-class of the general Partition class.
    
    A split is a 2-partition of a set of elements. It takes as input two sets
    of elements that form the split.
    
    Parameters
    ----------
    set1 : set[T]
        First set of elements in the split.
    set2 : set[T]
        Second set of elements in the split.
    
    Attributes
    ----------
    set1 : set[T]
        First set of elements in the split.
    set2 : set[T]
        Second set of elements in the split.
    elements : frozenset
        Set containing all elements from both sides of the split (inherited from Partition).
    
    Raises
    ------
    PhyloZooValueError
        If the sets overlap (i.e., the split is invalid).
    
    Examples
    --------
    >>> split = Split({1, 2}, {3, 4})
    >>> split.is_trivial()
    False
    >>> split.elements
    {1, 2, 3, 4}
    >>> split2 = Split({1}, {2, 3, 4})
    >>> split2.is_trivial()
    True
    """
    
    def __init__(self, set1: set[T], set2: set[T]) -> None:
        """
        Initialize a split.
        
        Parameters
        ----------
        set1 : set[T]
            First set of elements in the split.
        set2 : set[T]
            Second set of elements in the split.
        
        Raises
        ------
        PhyloZooValueError
            If the sets overlap (i.e., the split is invalid) or if either set is empty.
        """
        # Validate that neither set is empty
        if len(set1) == 0 or len(set2) == 0:
            raise PhyloZooValueError("Split sets cannot be empty")
        
        # Initialize parent Partition with the two parts first
        super().__init__([set1, set2])
        
        # Use canonical ordering from Partition to determine set1 and set2
        # Partition stores parts in canonical order, so we use that order
        parts_iter = iter(self._parts)
        set1_canonical = set(next(parts_iter))
        set2_canonical = set(next(parts_iter))
        
        # Set additional attributes using object.__setattr__ to bypass immutability
        object.__setattr__(self, 'set1', set1_canonical)
        object.__setattr__(self, 'set2', set2_canonical)
    
    def __repr__(self) -> str:
        """
        Return string representation of the split.
        
        Returns
        -------
        str
            String representation.
        """
        return f"Split({self.set1}, {self.set2})"

    @cached_property
    def is_trivial(self) -> bool:
        """
        Check if this is a trivial split.

        A trivial split is one where one of the sets has size 1.

        Returns
        -------
        bool
            True if the split is trivial, False otherwise.
        """
        return len(self.set1) == 1 or len(self.set2) == 1

    
