"""
Circular ordering module.

This module provides classes for working with circular orderings.
"""

from typing import List, Optional, Set


class CircularOrdering:
    """
    A circular ordering of elements.

    This is a placeholder class for circular ordering functionality.
    """

    def __init__(self, elements: Optional[List[str]] = None) -> None:
        """
        Initialize a circular ordering.

        Parameters
        ----------
        elements : Optional[List[str]], optional
            List of elements in circular order, by default None
        """
        self.elements: List[str] = elements or []

    def __len__(self) -> int:
        """
        Return the number of elements in the ordering.

        Returns
        -------
        int
            Number of elements
        """
        return len(self.elements)

    def __repr__(self) -> str:
        """
        Return string representation of the circular ordering.

        Returns
        -------
        str
            String representation
        """
        return f"CircularOrdering(elements={len(self.elements)})"


class CircularSetOrdering:
    """
    A circular ordering of sets.

    This is a placeholder class for circular set ordering functionality.
    """

    def __init__(self, sets: Optional[List[Set[str]]] = None) -> None:
        """
        Initialize a circular set ordering.

        Parameters
        ----------
        sets : Optional[List[Set[str]]], optional
            List of sets in circular order, by default None
        """
        self.sets: List[Set[str]] = sets or []

    def __len__(self) -> int:
        """
        Return the number of sets in the ordering.

        Returns
        -------
        int
            Number of sets
        """
        return len(self.sets)

    def __repr__(self) -> str:
        """
        Return string representation of the circular set ordering.

        Returns
        -------
        str
            String representation
        """
        return f"CircularSetOrdering(sets={len(self.sets)})"
