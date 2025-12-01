"""
Trinet set module.

This module provides classes for working with sets of trinets.
"""

from typing import List, Optional

from phylozoo.trinet import Trinet


class TrinetSet:
    """
    A set of trinets.

    This is a placeholder class for trinet set functionality.
    """

    def __init__(self, trinets: Optional[List[Trinet]] = None) -> None:
        """
        Initialize a trinet set.

        Parameters
        ----------
        trinets : Optional[List[Trinet]], optional
            List of trinets, by default None
        """
        self.trinets: List[Trinet] = trinets or []

    def add(self, trinet: Trinet) -> None:
        """
        Add a trinet to the set.

        Parameters
        ----------
        trinet : Trinet
            Trinet to add
        """
        self.trinets.append(trinet)

    def __len__(self) -> int:
        """
        Return the number of trinets in the set.

        Returns
        -------
        int
            Number of trinets
        """
        return len(self.trinets)

    def __repr__(self) -> str:
        """
        Return string representation of the trinet set.

        Returns
        -------
        str
            String representation
        """
        return f"TrinetSet(trinets={len(self.trinets)})"


class DenseTrinetSet(TrinetSet):
    """
    A dense trinet set.

    This is a placeholder class for dense trinet set functionality.
    """

    pass
