"""
Weighted split systems module.

This module provides classes for working with weighted split systems.
"""

from typing import TYPE_CHECKING

from .base import Split
from .splitsystem import SplitSystem

class WeightedSplitSystem(SplitSystem):
    """
    Class for a weighted split system: set of full splits with positive weights.
    
    A weighted split system is a function that maps each possible split on a set
    of elements to a weight. This implementation inherits from SplitSystem and only
    stores splits with positive weights. Zero-weight splits are not allowed.
    
    Parameters
    ----------
    splits : set[Split] | list[Split], optional
        Set or list of splits. If a list is provided, it will be converted
        to a set to ensure uniqueness. By default None (empty set).
    weights : dict[Split, float], optional
        Dictionary mapping splits to their weights. All weights must be positive
        (zero weights are not allowed). By default None (empty dict, which results
        in an empty system).
    
    Attributes
    ----------
    splits : frozenset[Split]
        Frozen set of splits with positive weights (read-only after initialization).
    elements : frozenset
        Frozen set containing all elements appearing in the splits (read-only).
    weights : dict[Split, float]
        Dictionary mapping splits to their weights (read-only after initialization).
    total_weight : float
        Sum of all weights in the system (read-only).
    
    Raises
    ------
    ValueError
        If not all splits cover the complete set of elements, if any weight
        is not positive (zero or negative), or if a weight is provided for a split
        not in the splits.
    
    Examples
    --------
    >>> split1 = Split({1, 2}, {3, 4})
    >>> split2 = Split({1, 3}, {2, 4})
    >>> weights = {split1: 2.5, split2: 1.0}
    >>> system = WeightedSplitSystem([split1, split2], weights=weights)
    >>> len(system)
    2
    >>> system.get_weight(split1)
    2.5
    >>> system.get_weight(Split({1, 4}, {2, 3}))  # Split not in system
    0.0
    >>> system.total_weight
    3.5
    
    Notes
    -----
    The weighted split system is immutable after initialization. Attempts to modify
    attributes will raise AttributeError.
    """
    
    __slots__ = ('_splits', '_elements', '_initialized', '_weights', '_total_weight')
    
    # I/O format configuration (override parent)
    _default_format = 'nexus'
    _supported_formats = ['nexus']
    
    def __init__(
        self,
        splits: (set[Split] | list[Split]) | None = None,
        weights: dict[Split, float] | None = None,
    ) -> None:
        """
        Initialize a weighted split system.
        
        Parameters
        ----------
        splits : set[Split] | list[Split], optional
            Set or list of splits. If a list is provided, it will be converted
            to a set to ensure uniqueness. By default None (empty set).
        weights : dict[Split, float], optional
            Dictionary mapping splits to their weights. All weights must be positive.
            Only splits with positive weights are stored in the system. By default None
            (empty dict, which results in an empty system).
        
        Raises
        ------
        ValueError
            If not all splits cover the complete set of elements, if any weight
            is non-positive, or if a weight is provided for a split not in the splits.
        """
        if weights is None:
            weights = {}
        
        # Validate that all weights are positive (zero and negative are not allowed)
        for split, weight in weights.items():
            if weight <= 0:
                raise ValueError(f"Weight must be positive, got {weight} for split {split}")
        
        # Determine splits: if splits is provided, use it; otherwise derive from weights.keys()
        if splits is None:
            # If no splits provided, derive from weights
            splits = set(weights.keys())
        elif not isinstance(splits, set):
            splits = set(splits)
        
        # Validate that all weights correspond to splits
        for split in weights:
            if split not in splits:
                raise ValueError(f"Weight provided for split {split} not in splits")
        
        # Only include splits that have weights (all weights are positive at this point)
        splits_with_weights = {split for split in splits if split in weights}
        
        # Initialize parent with only splits that have weights
        super().__init__(splits_with_weights)
        
        # Store weights
        # Use object.__setattr__ to bypass immutability check during initialization
        object.__setattr__(self, '_weights', dict(weights))
        object.__setattr__(self, '_total_weight', sum(weights.values()))
    
    @property
    def weights(self) -> dict[Split, float]:
        """
        Get the weights dictionary (read-only).
        
        Returns
        -------
        dict[Split, float]
            Dictionary mapping splits to their weights.
        """
        return self._weights.copy()  # Return a copy for immutability
    
    @property
    def total_weight(self) -> float:
        """
        Get the sum of all weights in the system.
        
        Returns
        -------
        float
            Sum of all weights.
        """
        return self._total_weight
    
    def get_weight(self, split: Split) -> float:
        """
        Get the weight of a split.
        
        Returns 0.0 if the split is not in the system (i.e., has no weight assigned).
        
        Parameters
        ----------
        split : Split
            Split to get the weight for.
        
        Returns
        -------
        float
            Weight of the split, or 0.0 if the split is not in the system.
        
        Raises
        ------
        ValueError
            If the split does not cover the same elements as the split system.
        
        Examples
        --------
        >>> split1 = Split({1, 2}, {3, 4})
        >>> split2 = Split({1, 3}, {2, 4})
        >>> system = WeightedSplitSystem([split1, split2], weights={split1: 2.5, split2: 1.0})
        >>> system.get_weight(split1)
        2.5
        >>> system.get_weight(Split({1, 4}, {2, 3}))  # Split not in system
        0.0
        """
        # Validate that split covers the same elements as the system
        if split.elements != self._elements:
            raise ValueError(
                f"Split {split} does not cover the same elements as the split system. "
                f"Expected elements: {self._elements}, got: {split.elements}"
            )
        
        return self._weights.get(split, 0.0)
    
    def __repr__(self) -> str:
        """
        Return string representation of the weighted split system.
        
        Returns
        -------
        str
            String representation including weights.
        """
        weights_str = ", ".join(f"{split}: {weight}" for split, weight in self._weights.items())
        return f"WeightedSplitSystem(splits={list(self._splits)}, weights={{{weights_str}}})"


def to_weightedsplitsystem(
    system: SplitSystem,
    default_weight: float = 1.0,
) -> WeightedSplitSystem:
    """
    Convert a SplitSystem to a WeightedSplitSystem.
    
    Assigns the same weight (default_weight) to each split in the system.
    
    Parameters
    ----------
    system : SplitSystem
        The split system to convert.
    default_weight : float, optional
        The weight to assign to each split. Must be positive. By default 1.0.
    
    Returns
    -------
    WeightedSplitSystem
        A weighted split system with all splits having the specified weight.
    
    Raises
    ------
    ValueError
        If default_weight is not positive.
    
    Examples
    --------
    >>> split1 = Split({1, 2}, {3, 4})
    >>> split2 = Split({1, 3}, {2, 4})
    >>> system = SplitSystem([split1, split2])
    >>> weighted = to_weightedsplitsystem(system, default_weight=2.0)
    >>> isinstance(weighted, WeightedSplitSystem)
    True
    >>> weighted.get_weight(split1)
    2.0
    >>> weighted.get_weight(split2)
    2.0
    """
    if default_weight <= 0:
        raise ValueError(f"default_weight must be positive, got {default_weight}")
    
    weights = {split: default_weight for split in system.splits}
    return WeightedSplitSystem(list(system.splits), weights=weights)

