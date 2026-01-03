"""
Split system classifications module.

This module provides functions for classifying split systems (e.g., checking
if they are compatible with trees).
"""

import itertools
from typing import TYPE_CHECKING

from .base import Split, is_compatible
from .splitsystem import SplitSystem

if TYPE_CHECKING:
    pass


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
    trivial_count = sum(1 for split in system if split.is_trivial())
    
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

