"""
Triplet base module.

This module provides the Triplet class for representing triplets (3-taxon rooted
trees). A triplet can be resolved (single trivial 1|2 split, where the singleton
side is a direct child of the root and the 2-element side forms a cherry below an
internal node) or unresolved (star tree, where the root has three direct leaf
children).
"""

from typing import TYPE_CHECKING

from ...utils.exceptions import PhyloZooValueError
from ..split.base import Split

if TYPE_CHECKING:
    from ..network.dnetwork import DirectedPhyNetwork


class Triplet:
    """
    Immutable triplet datatype representing a rooted tree on 3 taxa.

    A triplet can either have a single trivial 1|2 split representing a resolved
    rooted tree (the singleton side is a direct child of the root and the
    2-element side forms a cherry under an internal node), or be a star tree
    (represented by a set of 3 taxa) where the root has three direct leaf children.

    Parameters
    ----------
    split : Split | frozenset[str] | set[str]
        Either a trivial 1|2 split on 3 taxa, or a set/frozenset of exactly 3
        taxon labels for a star tree.

    Raises
    ------
    PhyloZooValueError
        If not exactly 3 taxa, or if split is not trivial (not 1|2).

    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> triplet = Triplet(Split({1}, {2, 3}))
    >>> triplet.taxa
    frozenset({1, 2, 3})
    >>> triplet.is_resolved()
    True
    >>> star_triplet = Triplet({1, 2, 3})
    >>> star_triplet.is_star()
    True
    """

    __slots__ = ('_taxa', '_split', '_initialized')

    def __init__(
        self,
        split: Split | frozenset[str] | set[str],
    ) -> None:
        """
        Initialize a triplet.

        Parameters
        ----------
        split : Split | frozenset[str] | set[str]
            Either a trivial 1|2 split on 3 taxa, or a set/frozenset of exactly
            3 taxon labels for a star tree.

        Raises
        ------
        PhyloZooValueError
            If not exactly 3 taxa, or if split is not trivial (not 1|2).
        """
        if isinstance(split, Split):
            # Extract taxa from split
            taxa_set = frozenset(split.elements)

            # Validate exactly 3 taxa
            if len(taxa_set) != 3:
                raise PhyloZooValueError(f"Split must have exactly 3 elements, got {len(taxa_set)}")

            # Validate split is 1|2 (trivial)
            if not split.is_trivial:
                raise PhyloZooValueError("Split must be a 1|2 split (trivial)")

            # Store split
            stored_split = split
        else:
            # Star tree: use provided taxa set
            taxa_set = frozenset(split)

            # Validate exactly 3 taxa
            if len(taxa_set) != 3:
                raise PhyloZooValueError(f"Taxa must have exactly 3 elements, got {len(taxa_set)}")

            # No split for star tree
            stored_split = None

        # Store as immutable
        object.__setattr__(self, '_taxa', taxa_set)
        object.__setattr__(self, '_split', stored_split)
        object.__setattr__(self, '_initialized', True)

    def __setattr__(self, name: str, value: any) -> None:
        """
        Prevent modification of attributes after initialization.

        Raises
        ------
        AttributeError
            If attempting to modify any attribute after initialization.
        """
        if not hasattr(self, '_initialized') or not self._initialized:
            super().__setattr__(name, value)
            return

        raise AttributeError(
            f"Cannot modify attribute '{name}'. Triplet is immutable."
        )

    @property
    def taxa(self) -> frozenset[str]:
        """
        Get the taxa of the triplet.

        Returns
        -------
        frozenset[str]
            The 3 taxon labels.
        """
        return self._taxa

    @property
    def split(self) -> Split | None:
        """
        Get the split of the triplet.

        Returns
        -------
        Split | None
            The split representing the topology, or None for star tree.
        """
        return self._split

    @property
    def outgroup(self) -> frozenset[str] | None:
        """
        Get the outgroup of the triplet (singleton side of the split).

        For a resolved triplet ``a|{b,c}``, returns ``frozenset({a})`` -- the
        single taxon that is a direct child of the root. For a star tree there
        is no outgroup and ``None`` is returned.

        Returns
        -------
        frozenset[str] | None
            The singleton side of the split, or None for star tree.
        """
        if self._split is None:
            return None
        # The trivial 1|2 split has one side of size 1 and one of size 2.
        if len(self._split.set1) == 1:
            return frozenset(self._split.set1)
        return frozenset(self._split.set2)

    @property
    def cherry(self) -> frozenset[str] | None:
        """
        Get the cherry of the triplet (2-element side of the split).

        For a resolved triplet ``a|{b,c}``, returns ``frozenset({b, c})`` -- the
        two taxa that share an internal parent below the root. For a star tree
        there is no cherry and ``None`` is returned.

        Returns
        -------
        frozenset[str] | None
            The 2-element side of the split, or None for star tree.
        """
        if self._split is None:
            return None
        if len(self._split.set1) == 2:
            return frozenset(self._split.set1)
        return frozenset(self._split.set2)

    def is_resolved(self) -> bool:
        """
        Check if the triplet is resolved (has a trivial split).

        Returns
        -------
        bool
            True if the triplet has a split (is resolved), False if it's a star tree.
        """
        return self._split is not None

    def is_star(self) -> bool:
        """
        Check if the triplet is a star tree (no split).

        Returns
        -------
        bool
            True if the triplet is a star tree, False if it's resolved.
        """
        return self._split is None

    def copy(self) -> 'Triplet':
        """
        Create a copy of the triplet.

        Returns
        -------
        Triplet
            A new Triplet instance with the same taxa and split.
        """
        if self._split is None:
            return Triplet(self._taxa)
        return Triplet(self._split)

    def to_network(self) -> 'DirectedPhyNetwork':
        """
        Convert the triplet to a DirectedPhyNetwork.

        For a resolved triplet ``a|{b,c}``, creates a rooted binary tree with
        a root node (with the outgroup ``a`` as one child and an internal node
        as the other child) and an internal node whose two children are the
        cherry leaves ``b`` and ``c``. For a star tree, creates a rooted tree
        with all three taxa as direct children of the root.

        Returns
        -------
        DirectedPhyNetwork
            A directed phylogenetic network representing the rooted triplet topology.
        """
        from ..network.dnetwork import DirectedPhyNetwork

        if self._split is None:
            # Star tree: 3 leaves directly under the root
            root_node = 4
            taxa_list = sorted(self._taxa)
            edges = [
                (root_node, taxa_list[0]),
                (root_node, taxa_list[1]),
                (root_node, taxa_list[2]),
            ]
            nodes = [
                (taxa_list[0], {'label': str(taxa_list[0])}),
                (taxa_list[1], {'label': str(taxa_list[1])}),
                (taxa_list[2], {'label': str(taxa_list[2])}),
            ]
        else:
            # Resolved triplet: root -> outgroup, root -> internal -> cherry
            outgroup_set = self.outgroup
            cherry_set = self.cherry
            outgroup_leaf = next(iter(outgroup_set))
            cherry_leaves = sorted(cherry_set)

            root_node = 4
            internal_node = 5

            edges = [
                (root_node, outgroup_leaf),
                (root_node, internal_node),
                (internal_node, cherry_leaves[0]),
                (internal_node, cherry_leaves[1]),
            ]
            nodes = [
                (outgroup_leaf, {'label': str(outgroup_leaf)}),
                (cherry_leaves[0], {'label': str(cherry_leaves[0])}),
                (cherry_leaves[1], {'label': str(cherry_leaves[1])}),
            ]

        return DirectedPhyNetwork(
            edges=edges,
            nodes=nodes,
        )

    def __hash__(self) -> int:
        """
        Return hash of the triplet.

        Returns
        -------
        int
            Hash value.
        """
        return hash((self._taxa, self._split))

    def __eq__(self, other: object) -> bool:
        """
        Check if two triplets are equal.

        Parameters
        ----------
        other : object
            Object to compare with.

        Returns
        -------
        bool
            True if triplets have same taxa and split, False otherwise.
        """
        if not isinstance(other, Triplet):
            return False
        return self._taxa == other._taxa and self._split == other._split

    def __repr__(self) -> str:
        """
        Return string representation of the triplet.

        Returns
        -------
        str
            String representation.
        """
        if self._split is None:
            return f"Triplet({set(self._taxa)})"
        return f"Triplet({self._split})"

    def __str__(self) -> str:
        """
        Return human-readable string representation of the triplet.

        For resolved triplets, displays as "Triplet(a | b c)" where ``a`` is
        the outgroup and ``b, c`` are the cherry leaves. For unresolved (star)
        triplets, displays as "Triplet(a b c)".

        Returns
        -------
        str
            Human-readable string representation.

        Examples
        --------
        >>> from phylozoo.core.split.base import Split
        >>> triplet = Triplet(Split({1}, {2, 3}))
        >>> str(triplet)
        'Triplet(1 | 2 3)'
        >>> star_triplet = Triplet({1, 2, 3})
        >>> str(star_triplet)
        'Triplet(1 2 3)'
        """
        if self._split is None:
            # Unresolved (star) triplet: show all taxa
            sorted_taxa = sorted(self._taxa, key=str)
            taxa_str = " ".join(str(taxon) for taxon in sorted_taxa)
            return f"Triplet({taxa_str})"
        else:
            # Resolved triplet: show outgroup | cherry
            outgroup_set = self.outgroup
            cherry_set = self.cherry
            sorted_outgroup = sorted(outgroup_set, key=str)
            sorted_cherry = sorted(cherry_set, key=str)
            outgroup_str = " ".join(str(taxon) for taxon in sorted_outgroup)
            cherry_str = " ".join(str(taxon) for taxon in sorted_cherry)
            return f"Triplet({outgroup_str} | {cherry_str})"
