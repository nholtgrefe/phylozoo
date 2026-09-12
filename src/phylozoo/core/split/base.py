"""
Splits module.

This module provides classes for working with phylogenetic splits. A split is a
2-partition {A, B} of a set of elements where A ∪ B equals the full set and
A ∩ B = ∅.
"""

from functools import cached_property
from typing import TypeVar

from ...utils.exceptions import PhyloZooValueError
from ..primitives.partition import Partition

T = TypeVar("T")


class Split(Partition[T]):
    """
    Class for 2-partitions of sets, child-class of the general Partition class.

    A split is a 2-partition of a set of elements. It takes as input two sets
    of elements that form the split.

    Parameters
    ----------
    set1 : set[T]
        First set of elements in the split.
    set2 : set[T]
        Second set of elements in the split.

    Raises
    ------
    PhyloZooValueError
        If the sets overlap (i.e., the split is invalid).

    Examples
    --------
    >>> split = Split({1, 2}, {3, 4})
    >>> split.is_trivial
    False
    >>> split.elements
    frozenset({1, 2, 3, 4})
    >>> split2 = Split({1}, {2, 3, 4})
    >>> split2.is_trivial
    True

    Attributes
    ----------
    set1, set2 : frozenset[T]
        The two sides of the split, in the canonical order used by
        :class:`~phylozoo.core.primitives.partition.Partition`. They are the very
        frozensets stored as the partition's parts (no copies), so a split costs about
        two frozensets of memory.
    elements : frozenset
        Set containing all elements from both sides of the split (inherited from Partition).
    """

    # Set in __init__ via object.__setattr__ (bypasses immutability of Partition);
    # declared here so static type checkers can see them.
    set1: frozenset[T]
    """First side of the split."""
    set2: frozenset[T]
    """Second side of the split."""

    def __init__(self, set1: set[T], set2: set[T]) -> None:
        """
        Initialize a split.

        Parameters
        ----------
        set1 : set[T]
            First set of elements in the split.
        set2 : set[T]
            Second set of elements in the split.

        Raises
        ------
        PhyloZooValueError
            If the sets overlap (i.e., the split is invalid) or if either set is empty.
        """
        # Validate that neither set is empty
        if len(set1) == 0 or len(set2) == 0:
            raise PhyloZooValueError("Split sets cannot be empty")

        # Initialize parent Partition with the two parts first
        super().__init__([set1, set2])

        # The sides are the partition's own canonical parts: sharing those frozensets
        # instead of copying them into sets is what keeps large split systems small.
        # Set via object.__setattr__ to bypass the immutability guard.
        object.__setattr__(self, "set1", self._parts[0])
        object.__setattr__(self, "set2", self._parts[1])

    @classmethod
    def _from_sides(
        cls, side1: frozenset[T], side2: frozenset[T], elements: frozenset[T]
    ) -> "Split[T]":
        """
        Build a split from two frozensets whose union is a known, shared ``elements`` set.

        For code that produces many splits over one taxon set (``induced_splits``,
        ``displayed_splits``), this avoids recomputing and re-storing the union of the
        sides for every split: all of them reference the single ``elements`` frozenset.
        The result is indistinguishable from ``Split(side1, side2)``.

        Parameters
        ----------
        side1, side2 : frozenset[T]
            The two sides. Must be non-empty, disjoint, and together equal to ``elements``.
        elements : frozenset[T]
            The full element set, shared between all splits built with it.

        Returns
        -------
        Split[T]
            The split.

        Raises
        ------
        PhyloZooValueError
            If a side is empty, the sides overlap, or they do not make up ``elements``.
        """
        if not side1 or not side2:
            raise PhyloZooValueError("Split sets cannot be empty")
        # Disjoint sides that both lie in ``elements`` and have its total size are
        # exactly a bipartition of it; none of these checks allocates a new set.
        if (
            len(side1) + len(side2) != len(elements)
            or not side1.isdisjoint(side2)
            or not side1 <= elements
            or not side2 <= elements
        ):
            raise PhyloZooValueError("Invalid partition: sets overlap")
        split = cls.__new__(cls)
        object.__setattr__(split, "_parts", Partition._canonical_form((side1, side2)))
        object.__setattr__(split, "_elements", elements)
        object.__setattr__(split, "_initialized", True)
        object.__setattr__(split, "set1", split._parts[0])
        object.__setattr__(split, "set2", split._parts[1])
        return split

    def __repr__(self) -> str:
        """
        Return string representation of the split.

        Returns
        -------
        str
            String representation.
        """
        # Shown as plain sets so the text stays ``Split({...}, {...})``, which recreates
        # the object, even though the sides are stored as frozensets.
        return f"Split({set(self.set1)}, {set(self.set2)})"

    @cached_property
    def is_trivial(self) -> bool:
        """
        Check if this is a trivial split.

        A trivial split is one where one of the sets has size 1.

        Returns
        -------
        bool
            True if the split is trivial, False otherwise.
        """
        return len(self.set1) == 1 or len(self.set2) == 1
