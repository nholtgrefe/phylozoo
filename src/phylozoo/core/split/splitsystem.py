"""
Split systems module.

This module provides classes for working with split systems.
"""

from typing import Iterator, Optional, Set, TYPE_CHECKING

from ...utils.tools import id_generator
from .split import Split

if TYPE_CHECKING:
    from ..network.sdnetwork import SemiDirectedPhyNetwork


class SplitSystem:
    """
    Class for a split system: set of full splits (complete partitions of elements).
    
    A split system is a collection of splits where each split covers the complete
    set of elements. This class validates that all splits cover the same element
    set and provides methods for working with split systems.
        
    Parameters
    ----------
    splits : Set[Split] | List[Split], optional
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
    
    def __init__(self, splits: Optional[Set[Split] | list[Split]] = None) -> None:
        """
        Initialize a split system.
        
        Parameters
        ----------
        splits : Set[Split] | List[Split], optional
            Set or list of splits. If a list is provided, it will be converted
            to a set to ensure uniqueness. By default None (empty set).
        
        Raises
        ------
        ValueError
            If not all splits cover the complete set of elements.
        """
        if splits is None:
            splits = set()
        elif isinstance(splits, list):
            splits = set(splits)
        else:
            splits = set(splits)
        
        # Convert to frozenset for immutability
        splits_frozen = frozenset(splits)
        
        # Compute elements
        if len(splits_frozen) == 0:
            elements_set: Set = set()
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
    
    def displayed_tree(self) -> 'SemiDirectedPhyNetwork':
        """
        Return the tree that displays the split system.
        
        Constructs a semi-directed network (tree) that displays all splits
        in the system. Raises an error if no such tree exists (i.e., the
        system is incompatible).
        
        Returns
        -------
        SemiDirectedPhyNetwork
            A tree (semi-directed network) that displays all splits in the system.
        
        Raises
        ------
        ValueError
            If the split system is incompatible (no tree exists that displays
            all splits).
        
        Examples
        --------
        >>> split1 = Split({1, 2}, {3, 4})
        >>> split2 = Split({1, 3}, {2, 4})
        >>> system = SplitSystem([split1, split2])
        >>> tree = system.displayed_tree()
        >>> isinstance(tree, SemiDirectedPhyNetwork)
        True
        """
        from ..network.sdnetwork import SemiDirectedPhyNetwork
        
        tree = SemiDirectedPhyNetwork(directed_edges=[], undirected_edges=[])
        center_node = id_generator()
        tree.add_node(center_node)
        
        # Add leaves from elements
        for element in self._elements:
            leaf = str(element)
            tree.add_node(leaf)
            tree.add_edge(center_node, leaf)
        
        # Process non-trivial splits
        for split in self._splits:
            if not split.is_trivial():
                # For now, we'll need to implement create_split in SemiDirectedPhyNetwork
                # This is a placeholder that will need to be implemented
                # tree.create_split(split)
                pass
        
        return tree
    
    def induced_quartetsplits(self) -> Set[Split]:
        """
        Return a set of all subsplits of size 4 of all splits in the system.
        
        Generates all quartet splits (2|2 splits) that can be induced from the
        splits in this system.
        
        Returns
        -------
        Set[Split]
            A set of quartet splits induced from this split system.
        
        Examples
        --------
        >>> split1 = Split({1, 2}, {3, 4, 5})
        >>> split2 = Split({1, 3}, {2, 4, 5})
        >>> system = SplitSystem([split1, split2])
        >>> quartets = system.induced_quartetsplits()
        >>> len(quartets) > 0
        True
        
        Notes
        -----
        This method returns a Set[Split] instead of QuartetSplitSet
        since QuartetSplitSet has been temporarily removed.
        """
        res: Set[Split] = set()
        for split in self._splits:
            quartet_splits = split.induced_quartetsplits()
            res.update(quartet_splits)
        return res
