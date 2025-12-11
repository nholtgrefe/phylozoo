"""
Distance matrix module.

This module provides classes for working with distance matrices.
"""

from typing import Dict, Tuple


class DistanceMatrix:
    """
    A distance matrix.

    This is a placeholder class for distance matrix functionality.
    """

    def __init__(self, distances: dict[tuple[str, str], float] | None = None) -> None:
        """
        Initialize a distance matrix.

        Parameters
        ----------
        distances : dict[tuple[str, str], float] | None, optional
            Dictionary mapping pairs of elements to distances, by default None
        """
        self.distances: dict[tuple[str, str], float] = distances or {}

    def get_distance(self, element1: str, element2: str) -> float | None:
        """
        Get the distance between two elements.

        Parameters
        ----------
        element1 : str
            First element
        element2 : str
            Second element

        Returns
        -------
        float | None
            Distance if found, None otherwise
        """
        return self.distances.get((element1, element2)) or self.distances.get((element2, element1))

    def set_distance(self, element1: str, element2: str, distance: float) -> None:
        """
        Set the distance between two elements.

        Parameters
        ----------
        element1 : str
            First element
        element2 : str
            Second element
        distance : float
            Distance value
        """
        self.distances[(element1, element2)] = distance

    def __repr__(self) -> str:
        """
        Return string representation of the distance matrix.

        Returns
        -------
        str
            String representation
        """
        return f"DistanceMatrix(distances={len(self.distances)})"
