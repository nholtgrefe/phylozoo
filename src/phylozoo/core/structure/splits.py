"""
Splits module.

This module provides classes for working with phylogenetic splits.
"""

import itertools
from typing import Set, TypeVar

from .partition import Partition

T = TypeVar('T')


class Split(Partition):
    """
    Class for 2-partitions of sets, child-class of the general Partition class.
    
    A split is a 2-partition of a set of elements. It takes as input two sets
    of elements that form the split.
    
    Parameters
    ----------
    set1 : Set[T]
        First set of elements in the split.
    set2 : Set[T]
        Second set of elements in the split.
    
    Attributes
    ----------
    set1 : Set[T]
        First set of elements in the split.
    set2 : Set[T]
        Second set of elements in the split.
    elements : Set[T]
        Set containing all elements from both sides of the split.
    
    Raises
    ------
    ValueError
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
    
    def __init__(self, set1: Set[T], set2: Set[T]) -> None:
        """
        Initialize a split.
        
        Parameters
        ----------
        set1 : Set[T]
            First set of elements in the split.
        set2 : Set[T]
            Second set of elements in the split.
        
        Raises
        ------
        ValueError
            If the sets overlap (i.e., the split is invalid) or if either set is empty.
        """
        # Validate that neither set is empty
        if len(set1) == 0 or len(set2) == 0:
            raise ValueError("Split sets cannot be empty")
        
        # Initialize parent Partition with the two parts first
        super().__init__([set1, set2])
        
        # Use canonical ordering from Partition to determine set1 and set2
        # Partition stores parts in canonical order, so we use that order
        canonical_parts = list(self._parts)
        set1_canonical = set(canonical_parts[0])
        set2_canonical = set(canonical_parts[1])
        
        # Set additional attributes using object.__setattr__ to bypass immutability
        object.__setattr__(self, 'set1', set1_canonical)
        object.__setattr__(self, 'set2', set2_canonical)
    
    @property
    def elements(self) -> Set[T]:
        """
        Get the set of all elements in the split.
        
        Returns
        -------
        Set[T]
            Set containing all elements from both sides of the split.
        """
        return self.set1 | self.set2
    
    def __repr__(self) -> str:
        """
        Return string representation of the split.
        
        Returns
        -------
        str
            String representation.
        """
        return f"Split({set(self.set1)}, {set(self.set2)})"
    
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
    
    def is_compatible(self, other: 'Split') -> bool:
        """
        Check if this split is compatible with another split.
        
        Two splits are compatible if they have the same set of elements, and one
        of the sets of this split is a subset of one of the sets of the other split
        (and hence the other set is a superset).
        
        Parameters
        ----------
        other : Split
            The split to check against.
        
        Returns
        -------
        bool
            True if this split is compatible with the other, False otherwise.
        
        Raises
        ------
        ValueError
            If 'other' is not a Split instance.
        
        Examples
        --------
        >>> split1 = Split({1, 2}, {3, 4})
        >>> split2 = Split({1}, {2, 3, 4})
        >>> split1.is_compatible(split2)
        True
        >>> split3 = Split({1, 2, 3}, {4})
        >>> split1.is_compatible(split3)
        True
        >>> split4 = Split({1, 3}, {2, 4})
        >>> split1.is_compatible(split4)
        False
        """
        if not isinstance(other, Split):
            raise ValueError("Not a Split instance")
        
        # First check that they have the same elements
        if self.elements != other.elements:
            return False
        
        # Check if one set of this split is a subset of one set of the other split
        # Case 1: self.set1 is subset of other.set1, so self.set2 is superset of other.set2
        if self.set1.issubset(other.set1) and other.set2.issubset(self.set2):
            return True
        
        # Case 2: self.set1 is subset of other.set2, so self.set2 is superset of other.set1
        if self.set1.issubset(other.set2) and other.set1.issubset(self.set2):
            return True
        
        # Case 3: other.set1 is subset of self.set1, so other.set2 is superset of self.set2
        if other.set1.issubset(self.set1) and self.set2.issubset(other.set2):
            return True
        
        # Case 4: other.set1 is subset of self.set2, so other.set2 is superset of self.set1
        if other.set1.issubset(self.set2) and self.set1.issubset(other.set2):
            return True
        
        return False
    
    def is_subsplit(self, other: 'Split') -> bool:
        """
        Check if this split is a subsplit of another split.
        
        A split is a subsplit of another if one of its sides is a subset of one side
        of the other split, and the other side of this split is a subset of the other
        side of the other split. For example, 12|56 is a subsplit of 123|456.
        
        Parameters
        ----------
        other : Split
            The split to check against.
        
        Returns
        -------
        bool
            True if this split is a subsplit of the other, False otherwise.
        
        Raises
        ------
        ValueError
            If 'other' is not a Split instance.
        
        Examples
        --------
        >>> split1 = Split({1, 2, 6}, {3, 4, 5})
        >>> split2 = Split({1, 2}, {3, 4})
        >>> split2.is_subsplit(split1)
        True
        >>> split3 = Split({1, 2}, {3, 4, 5})
        >>> split3.is_subsplit(split1)
        True
        """
        if not isinstance(other, Split):
            raise ValueError("Not a Split instance")
        
        return (self.set1.issubset(other.set1) and self.set2.issubset(other.set2)) or \
               (self.set1.issubset(other.set2) and self.set2.issubset(other.set1))
    
    def induced_quartetsplits(self, include_trivial: bool = False) -> Set['Split']:
        """
        Return a set of all subsplits of size 4 of the split.
        
        Generates all quartet splits (2|2 splits) that can be induced from this
        split by selecting 2 elements from each side.
        
        Parameters
        ----------
        include_trivial : bool, optional
            If True, also include trivial quartet splits (1|3 splits).
            By default False.
        
        Returns
        -------
        Set[Split[T]]
            A set of quartet splits induced from this split.
        
        Examples
        --------
        >>> split = Split({1, 2, 3}, {4, 5, 6})
        >>> quartets = split.induced_quartetsplits()
        >>> len(quartets) > 0
        True
        
        Notes
        -----
        This method returns a Set[Split[T]] instead of QuartetSplitSet
        since QuartetSplitSet has been temporarily removed.
        """
        res: list[Split] = []
        
        # Generate 2|2 splits
        for s1 in itertools.combinations(self.set1, 2):
            for s2 in itertools.combinations(self.set2, 2):
                split = Split(set(s1), set(s2))
                res.append(split)
        
        # Optionally include trivial splits (1|3)
        if include_trivial:
            for s1 in itertools.combinations(self.set1, 1):
                for s2 in itertools.combinations(self.set2, 3):
                    split = Split(set(s1), set(s2))
                    res.append(split)
            
            for s1 in itertools.combinations(self.set1, 3):
                for s2 in itertools.combinations(self.set2, 1):
                    split = Split(set(s1), set(s2))
                    res.append(split)
        
        return set(res)
