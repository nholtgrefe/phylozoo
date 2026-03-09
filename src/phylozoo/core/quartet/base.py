"""
Quartet base module.

This module provides the Quartet class for representing quartets (4-taxon unrooted
trees). A quartet can be resolved (single non-trivial split) or unresolved (star tree).
"""

from typing import TYPE_CHECKING

from ...utils.exceptions import PhyloZooValueError
from ..primitives.circular_ordering import CircularOrdering
from ..split.base import Split

if TYPE_CHECKING:
    from ..network.sdnetwork import SemiDirectedPhyNetwork


class Quartet:
    """
    Immutable quartet datatype representing an unrooted tree on 4 taxa.
    
    A quartet can either have a single non-trivial split (2|2 split) representing
    a resolved tree, or be a star tree (represented by a set of 4 taxa) where all
    taxa are equivalent.
    
    Parameters
    ----------
    split : Split | frozenset[str] | set[str]
        Either a 2|2 split on 4 taxa, or a set/frozenset of exactly 4 taxon labels
        for a star tree.
    
    Raises
    ------
    PhyloZooValueError
        If not exactly 4 taxa, or if split is trivial (not 2|2).
    
    Examples
    --------
    >>> from phylozoo.core.split.base import Split
    >>> quartet = Quartet(Split({1, 2}, {3, 4}))
    >>> quartet.taxa
    frozenset({1, 2, 3, 4})
    >>> quartet.is_resolved()
    True
    >>> star_quartet = Quartet({1, 2, 3, 4})
    >>> star_quartet.is_star()
    True
    """
    
    __slots__ = ('_taxa', '_split', '_initialized', '_circular_orderings_cache')
    
    def __init__(
        self,
        split: Split | frozenset[str] | set[str],
    ) -> None:
        """
        Initialize a quartet.
        
        Parameters
        ----------
        split : Split | frozenset[str] | set[str]
            Either a 2|2 split on 4 taxa, or a set/frozenset of exactly 4 taxon
            labels for a star tree.
        
        Raises
        ------
        PhyloZooValueError
            If not exactly 4 taxa, or if split is trivial (not 2|2).
        """
        if isinstance(split, Split):
            # Extract taxa from split
            taxa_set = frozenset(split.elements)
            
            # Validate split is 2|2 (not trivial)
            if split.is_trivial:
                raise PhyloZooValueError("Split must be a 2|2 split (non-trivial)")
            
            # Validate exactly 4 taxa
            if len(taxa_set) != 4:
                raise PhyloZooValueError(f"Split must have exactly 4 elements, got {len(taxa_set)}")
            
            # Store split
            stored_split = split
        else:
            # Star tree: use provided taxa set
            taxa_set = frozenset(split)
            
            # Validate exactly 4 taxa
            if len(taxa_set) != 4:
                raise PhyloZooValueError(f"Taxa must have exactly 4 elements, got {len(taxa_set)}")
            
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
            f"Cannot modify attribute '{name}'. Quartet is immutable."
        )
    
    @property
    def taxa(self) -> frozenset[str]:
        """
        Get the taxa of the quartet.
        
        Returns
        -------
        frozenset[str]
            The 4 taxon labels.
        """
        return self._taxa
    
    @property
    def split(self) -> Split | None:
        """
        Get the split of the quartet.
        
        Returns
        -------
        Split | None
            The split representing the topology, or None for star tree.
        """
        return self._split
    
    def is_resolved(self) -> bool:
        """
        Check if the quartet is resolved (has a non-trivial split).
        
        Returns
        -------
        bool
            True if the quartet has a split (is resolved), False if it's a star tree.
        """
        return self._split is not None
    
    def is_star(self) -> bool:
        """
        Check if the quartet is a star tree (no split).
        
        Returns
        -------
        bool
            True if the quartet is a star tree, False if it's resolved.
        """
        return self._split is None
    
    def copy(self) -> 'Quartet':
        """
        Create a copy of the quartet.
        
        Returns
        -------
        Quartet
            A new Quartet instance with the same taxa and split.
        """
        if self._split is None:
            return Quartet(self._taxa)
        return Quartet(self._split)
    
    @property
    def circular_orderings(self) -> frozenset[CircularOrdering]:
        """
        Get circular orderings that are congruent with this quartet.
        
        For a resolved quartet with split {a,b}|{c,d}, returns the two
        circular orderings where a and b are neighbors and c and d are
        neighbors: abcd and abdc.
        
        For a star tree (unresolved quartet), returns all three circular
        orderings: abcd, acbd, and adbc.
        
        The result is cached after first computation.
        
        Returns
        -------
        frozenset[CircularOrdering]
            Set of circular orderings congruent with this quartet.
        """
        if hasattr(self, '_circular_orderings_cache'):
            return self._circular_orderings_cache
        
        taxa_list = sorted(self._taxa)
        a, b, c, d = taxa_list
        
        if self._split is None:
            # Star tree: return all three distinct orderings (abcd, acbd, abdc)
            result = frozenset([
                CircularOrdering([a, b, c, d]),
                CircularOrdering([a, c, b, d]),
                CircularOrdering([a, b, d, c]),
            ])
        else:
            # Resolved quartet: return the two orderings where split pairs are neighbors
            set1 = sorted(self._split.set1)
            set2 = sorted(self._split.set2)
            
            # The two taxa from set1 and set2
            x1, x2 = set1
            y1, y2 = set2
            
            # Return abcd and abdc where a,b are from one set and c,d from the other
            result = frozenset([
                CircularOrdering([x1, x2, y1, y2]),
                CircularOrdering([x1, x2, y2, y1]),
            ])
        
        # Cache the result
        object.__setattr__(self, '_circular_orderings_cache', result)
        return result
    
    def to_network(self) -> 'SemiDirectedPhyNetwork':
        """
        Convert the quartet to a SemiDirectedPhyNetwork.
        
        For a resolved quartet, creates a tree with two internal nodes:
        one connecting the two taxa on each side of the split, and these
        internal nodes are connected. For a star tree, creates a tree with
        all four taxa connected to a single internal node.
        
        Returns
        -------
        SemiDirectedPhyNetwork
            A semi-directed phylogenetic network representing the quartet topology.
        """
        from ..network.sdnetwork import SemiDirectedPhyNetwork
        
        taxa_list = sorted(self._taxa)
        
        if self._split is None:
            # Star tree: all 4 taxa connected to a single internal node
            internal_node = 5
            undirected_edges = [
                (internal_node, taxa_list[0]),
                (internal_node, taxa_list[1]),
                (internal_node, taxa_list[2]),
                (internal_node, taxa_list[3]),
            ]
            nodes = [
                (taxa_list[0], {'label': str(taxa_list[0])}),
                (taxa_list[1], {'label': str(taxa_list[1])}),
                (taxa_list[2], {'label': str(taxa_list[2])}),
                (taxa_list[3], {'label': str(taxa_list[3])}),
            ]
        else:
            # Resolved quartet: two internal nodes, one for each side of the split
            set1 = sorted(self._split.set1)
            set2 = sorted(self._split.set2)
            
            internal_node1 = 5
            internal_node2 = 6
            
            undirected_edges = [
                (internal_node1, set1[0]),
                (internal_node1, set1[1]),
                (internal_node2, set2[0]),
                (internal_node2, set2[1]),
                (internal_node1, internal_node2),
            ]
            nodes = [
                (set1[0], {'label': str(set1[0])}),
                (set1[1], {'label': str(set1[1])}),
                (set2[0], {'label': str(set2[0])}),
                (set2[1], {'label': str(set2[1])}),
            ]
        
        return SemiDirectedPhyNetwork(
            undirected_edges=undirected_edges,
            nodes=nodes,
        )
    
    def __hash__(self) -> int:
        """
        Return hash of the quartet.
        
        Returns
        -------
        int
            Hash value.
        """
        return hash((self._taxa, self._split))
    
    def __eq__(self, other: object) -> bool:
        """
        Check if two quartets are equal.
        
        Parameters
        ----------
        other : object
            Object to compare with.
        
        Returns
        -------
        bool
            True if quartets have same taxa and split, False otherwise.
        """
        if not isinstance(other, Quartet):
            return False
        return self._taxa == other._taxa and self._split == other._split
    
    def __repr__(self) -> str:
        """
        Return string representation of the quartet.
        
        Returns
        -------
        str
            String representation.
        """
        if self._split is None:
            return f"Quartet({set(self._taxa)})"
        return f"Quartet({self._split})"
    
    def __str__(self) -> str:
        """
        Return human-readable string representation of the quartet.
        
        For resolved quartets, displays as "Quartet(a b | c d)".
        For unresolved (star) quartets, displays as "Quartet(a b c d)".
        
        Returns
        -------
        str
            Human-readable string representation.
        
        Examples
        --------
        >>> from phylozoo.core.split.base import Split
        >>> quartet = Quartet(Split({1, 2}, {3, 4}))
        >>> str(quartet)
        'Quartet(1 2 | 3 4)'
        >>> star_quartet = Quartet({1, 2, 3, 4})
        >>> str(star_quartet)
        'Quartet(1 2 3 4)'
        """
        if self._split is None:
            # Unresolved (star) quartet: show all taxa
            sorted_taxa = sorted(self._taxa, key=str)
            taxa_str = " ".join(str(taxon) for taxon in sorted_taxa)
            return f"Quartet({taxa_str})"
        else:
            # Resolved quartet: show split format
            sorted_set1 = sorted(self._split.set1, key=str)
            sorted_set2 = sorted(self._split.set2, key=str)
            set1_str = " ".join(str(taxon) for taxon in sorted_set1)
            set2_str = " ".join(str(taxon) for taxon in sorted_set2)
            return f"Quartet({set1_str} | {set2_str})"

