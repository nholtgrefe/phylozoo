"""
Splits module.

This module provides classes for working with phylogenetic splits.
"""

import itertools
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


def induced_quartetsplits(split: Split, include_trivial: bool = False) -> set[Split]:
    """
    Return a set of all subsplits of size 4 of the split.
    
    Generates all quartet splits (2|2 splits) that can be induced from this
    split by selecting 2 elements from each side.
    
    Parameters
    ----------
    split : Split[T]
        The split to generate quartet splits from.
    include_trivial : bool, optional
        If True, also include trivial quartet splits (1|3 splits).
        By default False.
    
    Returns
    -------
    set[Split[T]]
        A set of quartet splits induced from this split.
    
    Examples
    --------
    >>> split = Split({1, 2, 3}, {4, 5, 6})
    >>> quartets = induced_quartetsplits(split)
    >>> len(quartets) > 0
    True
    
    Notes
    -----
    This function returns a set[Split[T]] instead of QuartetSplitSet
    since QuartetSplitSet has been temporarily removed.
    """
    res: list[Split] = []
    
    # Generate 2|2 splits
    for s1 in itertools.combinations(split.set1, 2):
        for s2 in itertools.combinations(split.set2, 2):
            quartet_split = Split(set(s1), set(s2))
            res.append(quartet_split)
    
    # Optionally include trivial splits (1|3)
    if include_trivial:
        for s1 in itertools.combinations(split.set1, 1):
            for s2 in itertools.combinations(split.set2, 3):
                quartet_split = Split(set(s1), set(s2))
                res.append(quartet_split)
        
        for s1 in itertools.combinations(split.set1, 3):
            for s2 in itertools.combinations(split.set2, 1):
                quartet_split = Split(set(s1), set(s2))
                res.append(quartet_split)
    
    return set(res)


def is_compatible(split1: Split, split2: Split) -> bool:
    """
    Check if two splits are compatible.
    
    Two splits are compatible if they have the same set of elements, and one
    of the sets of one split is a subset of one of the sets of the other split
    (and hence the other set is a superset).
    
    Parameters
    ----------
    split1 : Split[T]
        The first split to check.
    split2 : Split[T]
        The second split to check.
    
    Returns
    -------
    bool
        True if the splits are compatible, False otherwise.
    
    Raises
    ------
    PhyloZooValueError
        If either argument is not a Split instance.
    
    Examples
    --------
    >>> split1 = Split({1, 2}, {3, 4})
    >>> split2 = Split({1}, {2, 3, 4})
    >>> is_compatible(split1, split2)
    True
    >>> split3 = Split({1, 2, 3}, {4})
    >>> is_compatible(split1, split3)
    True
    >>> split4 = Split({1, 3}, {2, 4})
    >>> is_compatible(split1, split4)
    False
    """
    if not isinstance(split1, Split):
        raise PhyloZooValueError("First argument must be a Split instance")
    if not isinstance(split2, Split):
        raise PhyloZooValueError("Second argument must be a Split instance")
    
    # First check that they have the same elements
    if split1.elements != split2.elements:
        return False
    
    # Check all 4 compatibility cases with early returns
    # Case 1: split1.set1 is subset of split2.set1, so split1.set2 is superset of split2.set2
    if split1.set1.issubset(split2.set1) and split2.set2.issubset(split1.set2):
        return True
    
    # Case 2: split1.set1 is subset of split2.set2, so split1.set2 is superset of split2.set1
    if split1.set1.issubset(split2.set2) and split2.set1.issubset(split1.set2):
        return True
    
    # Case 3: split2.set1 is subset of split1.set1, so split2.set2 is superset of split1.set2
    if split2.set1.issubset(split1.set1) and split1.set2.issubset(split2.set2):
        return True
    
    # Case 4: split2.set1 is subset of split1.set2, so split2.set2 is superset of split1.set1
    if split2.set1.issubset(split1.set2) and split1.set1.issubset(split2.set2):
        return True
    
    return False


def is_subsplit(split1: Split, split2: Split) -> bool:
    """
    Check if one split is a subsplit of another split.
    
    A split is a subsplit of another if one of its sides is a subset of one side
    of the other split, and the other side of this split is a subset of the other
    side of the other split. For example, 12|56 is a subsplit of 123|456.
    
    Parameters
    ----------
    split1 : Split[T]
        The split to check if it is a subsplit.
    split2 : Split[T]
        The split to check against.
    
    Returns
    -------
    bool
        True if split1 is a subsplit of split2, False otherwise.
    
    Raises
    ------
    PhyloZooValueError
        If either argument is not a Split instance.
    
    Examples
    --------
    >>> split1 = Split({1, 2, 6}, {3, 4, 5})
    >>> split2 = Split({1, 2}, {3, 4})
    >>> is_subsplit(split2, split1)
    True
    >>> split3 = Split({1, 2}, {3, 4, 5})
    >>> is_subsplit(split3, split1)
    True
    """
    if not isinstance(split1, Split):
        raise PhyloZooValueError("First argument must be a Split instance")
    if not isinstance(split2, Split):
        raise PhyloZooValueError("Second argument must be a Split instance")
    
    return (split1.set1.issubset(split2.set1) and split1.set2.issubset(split2.set2)) or \
           (split1.set1.issubset(split2.set2) and split1.set2.issubset(split2.set1))
