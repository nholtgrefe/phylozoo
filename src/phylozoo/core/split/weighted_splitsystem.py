"""
Weighted split systems module.

This module provides classes for working with weighted split systems. A weighted
split system assigns positive weights to each split in the system.
"""

from typing import TYPE_CHECKING

from ...utils.exceptions import PhyloZooValueError
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
    splits : set[Split] | list[Split] | dict[Split, float] | list[tuple[Split, float]] | None, optional
        Input splits with weights. Can be:
        - A set or list of splits (each assigned weight 1.0)
        - A dictionary mapping splits to their weights
        - A list of (split, weight) tuples
        By default None (empty system).
    
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
    PhyloZooValueError
        If not all splits cover the complete set of elements, if any weight
        is not positive (zero or negative), if duplicate splits are found, or if
        split elements don't match system elements.
    
    Examples
    --------
    From list of splits (weight 1.0 each):
    
    >>> split1 = Split({1, 2}, {3, 4})
    >>> split2 = Split({1, 3}, {2, 4})
    >>> system = WeightedSplitSystem([split1, split2])
    >>> system.get_weight(split1)
    1.0
    >>> system.get_weight(split2)
    1.0
    
    From dictionary with weights:
    
    >>> weights = {split1: 2.5, split2: 1.0}
    >>> system = WeightedSplitSystem(weights)
    >>> system.get_weight(split1)
    2.5
    >>> system.total_weight
    3.5
    
    From list of tuples:
    
    >>> system = WeightedSplitSystem([(split1, 0.8), (split2, 0.2)])
    >>> system.get_weight(split1)
    0.8
    
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
        splits: (
            set[Split]
            | list[Split]
            | dict[Split, float]
            | list[tuple[Split, float]]
            | None
        ) = None,
    ) -> None:
        """
        Initialize a weighted split system.
        
        Parameters
        ----------
        splits : set[Split] | list[Split] | dict[Split, float] | list[tuple[Split, float]] | None, optional
            Input splits with weights. Can be:
            - A set or list of splits (each assigned weight 1.0)
            - A dictionary mapping splits to their weights
            - A list of (split, weight) tuples
            By default None (empty system).
        
        Raises
        ------
        ValueError
            If not all splits cover the complete set of elements, or if any weight
            is non-positive.
        """
        if splits is None:
            splits = {}
        
        weights: dict[Split, float] = {}
        
        # Determine input type and extract splits and weights
        if isinstance(splits, dict):
            # Dictionary: splits are keys, values are weights
            weights = dict(splits)
        elif isinstance(splits, list):
            # List: check if it's list of tuples or list of splits
            if len(splits) == 0:
                weights = {}
            elif isinstance(splits[0], tuple):
                # List of (split, weight) tuples
                # Check for duplicate splits
                seen_splits: set[Split] = set()
                for split, weight in splits:
                    if split in seen_splits:
                        raise PhyloZooValueError(f"Duplicate split found: {split}")
                    seen_splits.add(split)
                weights = dict(splits)
            else:
                # List of splits: assign weight 1.0 to each
                # Check for duplicate splits
                seen_splits: set[Split] = set()
                for split in splits:
                    if split in seen_splits:
                        raise PhyloZooValueError(f"Duplicate split found: {split}")
                    seen_splits.add(split)
                weights = {split: 1.0 for split in splits}
        elif isinstance(splits, set):
            # Set of splits: assign weight 1.0 to each
            # Sets automatically handle uniqueness, so no need to check
            weights = {split: 1.0 for split in splits}
        else:
            raise PhyloZooValueError(
                f"Expected set[Split], list[Split], dict[Split, float], "
                f"list[tuple[Split, float]], or None, got {type(splits)}"
            )
        
        # Validate that all weights are positive (zero and negative are not allowed)
        for split, weight in weights.items():
            if weight <= 0:
                raise PhyloZooValueError(f"Weight must be positive, got {weight} for split {split}")
        
        # Only include splits that have weights (all weights are positive at this point)
        splits_with_weights = set(weights.keys())
        
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
        PhyloZooValueError
            If the split does not cover the same elements as the split system.
        
        Examples
        --------
        >>> split1 = Split({1, 2}, {3, 4})
        >>> split2 = Split({1, 3}, {2, 4})
        >>> system = WeightedSplitSystem({split1: 2.5, split2: 1.0})
        >>> system.get_weight(split1)
        2.5
        >>> system.get_weight(Split({1, 4}, {2, 3}))  # Split not in system
        0.0
        """
        # Validate that split covers the same elements as the system
        if split.elements != self._elements:
            raise PhyloZooValueError(
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
            String representation that can be used to recreate the object.
        """
        if len(self._weights) == 0:
            return "WeightedSplitSystem()"
        
        weights_str = ", ".join(f"{repr(split)}: {weight}" for split, weight in self._weights.items())
        return f"WeightedSplitSystem({{{weights_str}}})"
    
    def __str__(self) -> str:
        """
        Return human-readable string representation of the weighted split system.
        
        Displays the weighted split system showing all splits with their weights,
        one per line. No truncation is applied.
        
        Returns
        -------
        str
            Human-readable string representation.
        
        Examples
        --------
        >>> split1 = Split({1, 2}, {3, 4})
        >>> split2 = Split({1, 3}, {2, 4})
        >>> system = WeightedSplitSystem({split1: 2.5, split2: 1.0})
        >>> str(system)
        'WeightedSplitSystem({\\n  Split(1 2 | 3 4): 2.5,\\n  Split(1 3 | 2 4): 1.0\\n})'
        """
        n = len(self._splits)
        if n == 0:
            return "WeightedSplitSystem({})"
        
        # Sort splits for consistent display
        sorted_splits = sorted(self._splits, key=lambda s: (str(s.set1), str(s.set2)))
        
        # Show all splits with weights, one per line
        splits_lines = [
            f"  {split}: {self._weights[split]}," for split in sorted_splits
        ]
        # Remove trailing comma from last line
        if splits_lines:
            splits_lines[-1] = splits_lines[-1].rstrip(',')
        return f"WeightedSplitSystem({{\n" + "\n".join(splits_lines) + "\n})"


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
    PhyloZooValueError
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
        raise PhyloZooValueError(f"default_weight must be positive, got {default_weight}")
    
    weights = {split: default_weight for split in system.splits}
    return WeightedSplitSystem(weights)

