"""
Split systems module.

This module provides classes for working with split systems.
"""

from typing import Iterator, TYPE_CHECKING

from ...utils.io import IOMixin
from .base import Split

class SplitSystem(IOMixin):
    """
    Class for a split system: set of full splits (complete partitions of elements).
    
    A split system is a collection of splits where each split covers the complete
    set of elements. This class validates that all splits cover the same element
    set and provides methods for working with split systems.
        
    Parameters
    ----------
    splits : set[Split] | list[Split], optional
        Set or list of splits. If a list is provided, it will be converted
        to a set to ensure uniqueness. By default None (empty set).
    
    Attributes
    ----------
    splits : frozenset[Split]
        Frozen set of splits (read-only after initialization).
    elements : frozenset
        Frozen set containing all elements appearing in the splits (read-only).
    
    Raises
    ------
    ValueError
        If not all splits cover the complete set of elements (i.e., not a set
        of full splits).
    
    Examples
    --------
    >>> split1 = Split({1, 2}, {3, 4})
    >>> split2 = Split({1, 3}, {2, 4})
    >>> system = SplitSystem([split1, split2])
    >>> len(system)
    2
    >>> system.elements == {1, 2, 3, 4}
    True
    >>> split1 in system
    True
    
    Notes
    -----
    The split system is immutable after initialization. Attempts to modify attributes
    will raise AttributeError.
    """
    
    __slots__ = ('_splits', '_elements', '_initialized')
    
    # I/O format configuration
    _default_format = 'nexus'
    _supported_formats = ['nexus']
    
    def __init__(self, splits: (set[Split] | list[Split]) | None = None) -> None:
        """
        Initialize a split system.
        
        Parameters
        ----------
        splits : set[Split] | list[Split], optional
            Set or list of splits. If a list is provided, it will be converted
            to a set to ensure uniqueness. By default None (empty set).
        
        Raises
        ------
        ValueError
            If not all splits cover the complete set of elements.
        """
        if splits is None:
            splits = set()
        elif not isinstance(splits, set):
            splits = set(splits)
        # If it's already a set, use it directly
        
        # Convert to frozenset for immutability
        splits_frozen = frozenset(splits)
        
        # Compute elements
        if len(splits_frozen) == 0:
            elements_set: set = set()
        else:
            # Get elements from first split
            first_split = next(iter(splits_frozen))
            elements_set = first_split.elements
        
        # Validate that all splits cover the same element set
        for split in splits_frozen:
            if split.elements != elements_set:
                raise ValueError("Not a set of full splits.")
        
        # Store as frozen sets for immutability
        self._splits: frozenset[Split] = splits_frozen
        self._elements: frozenset = frozenset(elements_set)
        self._initialized: bool = True
    
    def __setattr__(self, name: str, value: any) -> None:
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
        raise AttributeError(
            f"Cannot modify attribute '{name}'. SplitSystem is immutable."
        )
    
    @property
    def splits(self) -> frozenset[Split]:
        """
        Get the splits in the system (read-only).
        
        Returns
        -------
        frozenset[Split]
            Frozen set of splits.
        """
        return self._splits
    
    @property
    def elements(self) -> frozenset:
        """
        Get the elements in the system (read-only).
        
        Returns
        -------
        frozenset
            Frozen set of all elements appearing in the splits.
        """
        return self._elements
    
    def __repr__(self) -> str:
        """
        Return string representation of the split system.
        
        Returns
        -------
        str
            String representation.
        """
        return f"SplitSystem({list(self._splits)})"
    
    def __str__(self) -> str:
        """
        Return human-readable string representation of the split system.
        
        Displays the split system showing all splits, one per line.
        No truncation is applied.
        
        Returns
        -------
        str
            Human-readable string representation.
        
        Examples
        --------
        >>> split1 = Split({1, 2}, {3, 4})
        >>> split2 = Split({1, 3}, {2, 4})
        >>> system = SplitSystem([split1, split2])
        >>> str(system)
        'SplitSystem({\\n  Split(1 2 | 3 4),\\n  Split(1 3 | 2 4)\\n})'
        """
        n = len(self._splits)
        if n == 0:
            return "SplitSystem({})"
        
        # Sort splits for consistent display
        sorted_splits = sorted(self._splits, key=lambda s: (str(s.set1), str(s.set2)))
        
        # Show all splits, one per line
        splits_lines = [f"  {split}," for split in sorted_splits]
        # Remove trailing comma from last line
        if splits_lines:
            splits_lines[-1] = splits_lines[-1].rstrip(',')
        return f"SplitSystem({{\n" + "\n".join(splits_lines) + "\n})"
    
    def __iter__(self) -> Iterator[Split]:
        """
        Return an iterator over the splits.
        
        Returns
        -------
        Iterator[Split]
            Iterator over splits.
        """
        return iter(self._splits)
    
    def __contains__(self, split: Split) -> bool:
        """
        Check if a split is in the system.
        
        Parameters
        ----------
        split : Split
            Split to check.
        
        Returns
        -------
        bool
            True if the split is in the system, False otherwise.
        """
        return split in self._splits
    
    def __len__(self) -> int:
        """
        Return the number of splits in the system.
        
        Returns
        -------
        int
            Number of splits.
        """
        return len(self._splits)
    
    def induced_quartetsplits(self) -> set[Split]:
        """Stub for induced_quartetsplits function."""
        return set()
