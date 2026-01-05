"""
Quarnet module.

This module provides classes for working with quarnets (4-leaf networks).
"""

from typing import TYPE_CHECKING, Mapping

from ...core.quartet.qprofile import QuartetProfile
from ...core.primitives.circular_ordering import CircularOrdering

if TYPE_CHECKING:
    from ...core.quartet.base import Quartet


class FourCycle(QuartetProfile):
    """
    A four-cycle quarnet with a reticulation leaf.
    
    A FourCycle represents a quartet profile with exactly two resolved quartets
    that form a circular ordering. The quarnet has a reticulation (hybrid node)
    with one leaf below it (the reticulation leaf).
    
    Parameters
    ----------
    quartets : dict[Quartet, float] | Mapping[Quartet, float] | list[Quartet] | list[tuple[Quartet, float]]
        Input quartets. Must contain exactly 2 resolved quartets that form
        a circular ordering. Can be:
        - A dictionary mapping quartets to weights
        - A list of quartets (each assigned weight 1.0)
        - A list of (quartet, weight) tuples
    reticulation_leaf : str
        The leaf label that is below the reticulation in the four-cycle.
        Must be one of the 4 taxa in the quartets.
    
    Attributes
    ----------
    reticulation_leaf : str
        The leaf label below the reticulation.
    circular_ordering : CircularOrdering | None
        The circular ordering induced by the two quartets, or None if
        no such ordering exists.
    
    Raises
    ------
    ValueError
        If quartets is empty, if quartets have different taxa, if any
        weight is non-positive, if there are not exactly 2 resolved quartets,
        if the quartets do not form a circular ordering, or if the
        reticulation_leaf is not in the taxa.
    
    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> from phylozoo.core.quartet.base import Quartet
    >>> q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
    >>> q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
    >>> fourcycle = FourCycle([q1, q2], reticulation_leaf='A')
    >>> fourcycle.reticulation_leaf
    'A'
    >>> len(fourcycle)
    2
    """
    
    __slots__ = ('_reticulation_leaf',)
    
    def __init__(
        self,
        quartets: (
            dict['Quartet', float]
            | Mapping['Quartet', float]
            | list['Quartet']
            | list[tuple['Quartet', float]]
        ),
        reticulation_leaf: str,
    ) -> None:
        """
        Initialize a four-cycle quarnet.
        
        Parameters
        ----------
        quartets : dict[Quartet, float] | Mapping[Quartet, float] | list[Quartet] | list[tuple[Quartet, float]]
            Input quartets. Must contain exactly 2 resolved quartets that form
            a circular ordering. Can be:
            - A dictionary mapping quartets to weights
            - A list of quartets (each assigned weight 1.0)
            - A list of (quartet, weight) tuples
        reticulation_leaf : str
            The leaf label that is below the reticulation in the four-cycle.
            Must be one of the 4 taxa in the quartets.
        
        Raises
        ------
        ValueError
            If quartets is empty, if quartets have different taxa, if any
            weight is non-positive, if there are not exactly 2 resolved quartets,
            if the quartets do not form a circular ordering, or if the
            reticulation_leaf is not in the taxa.
        """
        # Initialize parent QuartetProfile
        super().__init__(quartets)
        
        # Validate that we have exactly 2 resolved quartets
        resolved_quartets = [q for q in self._quartets if q.is_resolved()]
        if len(resolved_quartets) != 2:
            raise ValueError(
                f"FourCycle must contain exactly 2 resolved quartets, "
                f"got {len(resolved_quartets)}"
            )
        
        # Validate that a circular ordering exists
        circular_orderings = self.circular_orderings
        if circular_orderings is None or len(circular_orderings) == 0:
            raise ValueError(
                "FourCycle quartets must form a circular ordering. "
                "The two resolved quartets do not have a common circular ordering."
            )
        
        # Validate that reticulation_leaf is in the taxa
        if reticulation_leaf not in self._taxa:
            raise ValueError(
                f"Reticulation leaf '{reticulation_leaf}' must be one of the "
                f"taxa: {set(self._taxa)}"
            )
        
        # Store reticulation leaf
        object.__setattr__(self, '_reticulation_leaf', reticulation_leaf)
    
    @property
    def reticulation_leaf(self) -> str:
        """
        Get the reticulation leaf.
        
        Returns
        -------
        str
            The leaf label below the reticulation.
        """
        return self._reticulation_leaf
    
    @property
    def circular_ordering(self) -> CircularOrdering | None:
        """
        Get the circular ordering induced by the two quartets.
        
        Returns
        -------
        CircularOrdering | None
            The circular ordering, or None if no ordering exists.
        """
        orderings = self.circular_orderings
        if orderings is None or len(orderings) == 0:
            return None
        # Return the first (and typically only) circular ordering
        return next(iter(orderings))
    
    def __repr__(self) -> str:
        """
        Return string representation of the four-cycle.
        
        Returns
        -------
        str
            String representation.
        """
        return (
            f"FourCycle(taxa={set(self._taxa)}, "
            f"reticulation_leaf='{self._reticulation_leaf}', "
            f"quartets={dict(self._quartets)})"
        )
