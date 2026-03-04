"""
Split system classifications module.

This module provides functions for classifying split systems, such as checking
pairwise compatibility.
"""

import itertools
from typing import TYPE_CHECKING, TypeVar

from ...utils.exceptions import PhyloZooValueError
from .base import Split
from .splitsystem import SplitSystem

T = TypeVar('T')

if TYPE_CHECKING:
    pass


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
    >>> split3 = Split({1, 3}, {2, 4})
    >>> is_subsplit(split3, split1)
    False
    """
    if not isinstance(split1, Split):
        raise PhyloZooValueError("First argument must be a Split instance")
    if not isinstance(split2, Split):
        raise PhyloZooValueError("Second argument must be a Split instance")

    # Check subsplit condition
    # split1 is a subsplit of split2 if each side of split1 is a subset of
    # one of the sides of split2.
    return (
        (split1.set1.issubset(split2.set1) and split1.set2.issubset(split2.set2)) or
        (split1.set1.issubset(split2.set2) and split1.set2.issubset(split2.set1))
    )




def is_pairwise_compatible(system: SplitSystem) -> bool:
    """
    Check if all pairs of splits in the system are compatible.
    
    A split system is pairwise compatible if every pair of splits in the system
    is compatible with each other.
    
    Parameters
    ----------
    system : SplitSystem
        The split system to check.
    
    Returns
    -------
    bool
        True if all pairs of splits are compatible, False otherwise.
    
    Examples
    --------
    >>> from phylozoo.core.split import Split, SplitSystem
    >>> split1 = Split({1, 2}, {3, 4})
    >>> split2 = Split({1}, {2, 3, 4})
    >>> split3 = Split({1, 2, 3}, {4})
    >>> system = SplitSystem([split1, split2, split3])
    >>> is_pairwise_compatible(system)
    True
    >>> split4 = Split({1, 3}, {2, 4})
    >>> system2 = SplitSystem([split1, split4])
    >>> is_pairwise_compatible(system2)
    False
    """
    # Use itertools.combinations to check all pairs
    for split1, split2 in itertools.combinations(system, 2):
        if not is_compatible(split1, split2):
            return False
    
    return True


def has_all_trivial_splits(system: SplitSystem) -> bool:
    """
    Check if the split system contains all trivial splits.
    
    For a split system with n elements, there should be n trivial splits,
    where each trivial split has one element in one set and all other n-1
    elements in the other set.
    
    Parameters
    ----------
    system : SplitSystem
        The split system to check.
    
    Returns
    -------
    bool
        True if all trivial splits are present, False otherwise.
    
    Examples
    --------
    >>> from phylozoo.core.split import Split, SplitSystem
    >>> # System with 3 elements should have 3 trivial splits
    >>> split1 = Split({1}, {2, 3})
    >>> split2 = Split({2}, {1, 3})
    >>> split3 = Split({3}, {1, 2})
    >>> system = SplitSystem([split1, split2, split3])
    >>> has_all_trivial_splits(system)
    True
    >>> # Missing one trivial split
    >>> system2 = SplitSystem([split1, split2])
    >>> has_all_trivial_splits(system2)
    False
    """
    if len(system.elements) == 0:
        return True  # Empty system has all trivial splits (trivially)
    
    # Count trivial splits in the system
    trivial_count = sum(1 for split in system if split.is_trivial)
    
    # For n elements, there should be exactly n trivial splits
    return trivial_count == len(system.elements)


def is_tree_compatible(system: SplitSystem) -> bool:
    """
    Check if a split system is compatible with a tree.
    
    A split system is tree-compatible if:
    1. All pairs of splits are compatible (pairwise compatible)
    2. All trivial splits are present in the system
    
    Parameters
    ----------
    system : SplitSystem
        The split system to check.
    
    Returns
    -------
    bool
        True if the system is compatible with a tree, False otherwise.
    
    Examples
    --------
    >>> from phylozoo.core.split import Split, SplitSystem
    >>> # Tree-compatible system
    >>> split1 = Split({1}, {2, 3, 4})
    >>> split2 = Split({2}, {1, 3, 4})
    >>> split3 = Split({3}, {1, 2, 4})
    >>> split4 = Split({4}, {1, 2, 3})
    >>> split5 = Split({1, 2}, {3, 4})
    >>> system = SplitSystem([split1, split2, split3, split4, split5])
    >>> is_tree_compatible(system)
    True
    >>> # Incompatible system (splits conflict)
    >>> split6 = Split({1, 3}, {2, 4})
    >>> system2 = SplitSystem([split1, split2, split3, split4, split5, split6])
    >>> is_tree_compatible(system2)
    False
    """
    return is_pairwise_compatible(system) and has_all_trivial_splits(system)

