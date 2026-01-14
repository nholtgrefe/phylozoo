"""
Squirrel quartet profile module.

This module provides the SqQuartetProfile class for working with squirrel quartet profiles.
"""

from typing import TYPE_CHECKING, Mapping

from ...core.quartet.qprofile import QuartetProfile
from ...core.primitives.circular_ordering import CircularOrdering
from ...utils.exceptions import PhyloZooValueError

if TYPE_CHECKING:
    from ...core.quartet.base import Quartet


class SqQuartetProfile(QuartetProfile):
    """
    Squirrel quartet profile.
    
    A SqQuartetProfile is a quartet profile that contains either one or two
    resolved quartets (no unresolved or other options). If it contains two
    quartets (forming a cycle), it can optionally have a reticulation leaf.
    
    Parameters
    ----------
    quartets : dict[Quartet, float] | Mapping[Quartet, float] | list[Quartet] | list[tuple[Quartet, float]]
        Input quartets. Must contain exactly 1 or 2 resolved quartets. Can be:
        - A dictionary mapping quartets to weights
        - A list of quartets (each assigned weight 1.0)
        - A list of (quartet, weight) tuples
    reticulation_leaf : str | None, optional
        The leaf label that is below the reticulation in a four-cycle.
        Only valid when the profile contains exactly 2 quartets.
        Must be one of the 4 taxa in the quartets. By default None.
    
    Attributes
    ----------
    reticulation_leaf : str | None
        The leaf label below the reticulation (if profile has 2 quartets),
        or None if profile has 1 quartet.
    circular_ordering : CircularOrdering | None
        The circular ordering induced by the quartets (if 2 quartets),
        or None if 1 quartet or no ordering exists.
    
    Raises
    ------
    PhyloZooValueError
        If quartets is empty, if quartets have different taxa, if any
        weight is non-positive, if there are not exactly 1 or 2 resolved quartets,
        if there are unresolved quartets, if reticulation_leaf is provided
        but profile has only 1 quartet, if the quartets do not form a circular
        ordering when there are 2 quartets, or if the reticulation_leaf is
        not in the taxa.
    
    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> from phylozoo.core.quartet.base import Quartet
    >>> # Single quartet profile
    >>> q1 = Quartet(Split({'A', 'B'}, {'C', 'D'}))
    >>> sq_profile = SqQuartetProfile([q1])
    >>> len(sq_profile)
    1
    >>> sq_profile.reticulation_leaf is None
    True
    >>> # Two quartets (cycle) with reticulation leaf
    >>> q2 = Quartet(Split({'A', 'C'}, {'B', 'D'}))
    >>> sq_profile2 = SqQuartetProfile([q1, q2], reticulation_leaf='A')
    >>> sq_profile2.reticulation_leaf
    'A'
    >>> len(sq_profile2)
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
        reticulation_leaf: str | None = None,
    ) -> None:
        """
        Initialize a squirrel quartet profile.
        
        Parameters
        ----------
        quartets : dict[Quartet, float] | Mapping[Quartet, float] | list[Quartet] | list[tuple[Quartet, float]]
            Input quartets. Must contain exactly 1 or 2 resolved quartets. Can be:
            - A dictionary mapping quartets to weights
            - A list of quartets (each assigned weight 1.0)
            - A list of (quartet, weight) tuples
        reticulation_leaf : str | None, optional
            The leaf label that is below the reticulation in a four-cycle.
            Only valid when the profile contains exactly 2 quartets.
            Must be one of the 4 taxa in the quartets. By default None.
        
        Raises
        ------
        ValueError
            If quartets is empty, if quartets have different taxa, if any
            weight is non-positive, if there are not exactly 1 or 2 resolved quartets,
            if there are unresolved quartets, if reticulation_leaf is provided
            but profile has only 1 quartet, if the quartets do not form a circular
            ordering when there are 2 quartets, or if the reticulation_leaf is
            not in the taxa.
        """
        # Initialize parent QuartetProfile
        super().__init__(quartets)
        
        # Validate that we have exactly 1 or 2 resolved quartets
        resolved_quartets = [q for q in self._quartets if q.is_resolved()]
        unresolved_quartets = [q for q in self._quartets if not q.is_resolved()]
        
        if len(unresolved_quartets) > 0:
            raise PhyloZooValueError(
                f"SqQuartetProfile cannot contain unresolved quartets, "
                f"got {len(unresolved_quartets)} unresolved quartets"
            )
        
        if len(resolved_quartets) not in (1, 2):
            raise PhyloZooValueError(
                f"SqQuartetProfile must contain exactly 1 or 2 resolved quartets, "
                f"got {len(resolved_quartets)}"
            )
        
        # If 2 quartets, validate circular ordering
        if len(resolved_quartets) == 2:
            circular_orderings = self.circular_orderings
            if circular_orderings is None or len(circular_orderings) == 0:
                raise PhyloZooValueError(
                    "SqQuartetProfile with 2 quartets must form a circular ordering. "
                    "The two resolved quartets do not have a common circular ordering."
                )
            
            # Validate reticulation_leaf if provided
            if reticulation_leaf is not None:
                if reticulation_leaf not in self._taxa:
                    raise PhyloZooValueError(
                        f"Reticulation leaf '{reticulation_leaf}' must be one of the "
                        f"taxa: {set(self._taxa)}"
                    )
        else:
            # 1 quartet: reticulation_leaf should be None
            if reticulation_leaf is not None:
                raise PhyloZooValueError(
                    "reticulation_leaf can only be provided when profile has 2 quartets, "
                    f"got profile with {len(resolved_quartets)} quartet(s)"
                )
        
        # Store reticulation leaf
        object.__setattr__(self, '_reticulation_leaf', reticulation_leaf)
    
    @property
    def reticulation_leaf(self) -> str | None:
        """
        Get the reticulation leaf.
        
        Returns
        -------
        str | None
            The leaf label below the reticulation (if profile has 2 quartets),
            or None if profile has 1 quartet.
        """
        return self._reticulation_leaf
    
    @property
    def circular_ordering(self) -> CircularOrdering | None:
        """
        Get the circular ordering induced by the quartets.
        
        Returns
        -------
        CircularOrdering | None
            The circular ordering (if 2 quartets), or None if 1 quartet or
            no ordering exists.
        """
        if len(self._quartets) == 2:
            orderings = self.circular_orderings
            if orderings is None or len(orderings) == 0:
                return None
            # Return the first (and typically only) circular ordering
            return next(iter(orderings))
        return None
    
    def __repr__(self) -> str:
        """
        Return string representation of the squirrel quartet profile.
        
        Returns
        -------
        str
            String representation.
        """
        ret_leaf_str = f"reticulation_leaf='{self._reticulation_leaf}'" if self._reticulation_leaf else "reticulation_leaf=None"
        return (
            f"SqQuartetProfile(taxa={set(self._taxa)}, "
            f"{ret_leaf_str}, "
            f"quartets={dict(self._quartets)})"
        )

