"""
Semi-directed network module.

This module provides classes and functions for working with semi-directed phylogenetic networks.
"""

import warnings
from functools import cached_property
from typing import Any, Dict, Iterator, List, Optional, Set, Tuple, TypeVar, Union

from .m_phynetwork import MixedPhyNetwork

T = TypeVar('T')


class SemiDirectedPhyNetwork(MixedPhyNetwork):
    """
    A semi-directed phylogenetic network.
    
    A SemiDirectedPhyNetwork is a mixed phylogenetic network where:
    - Exactly the non-hybrid edges are undirected
    - All hybrid edges remain directed
    
    This is a subclass of MixedPhyNetwork with additional constraints on edge directions.
    
    Parameters
    ----------
    directed_edges : List[Union[Tuple[T, T], Tuple[T, T, int], Dict[str, Any]]], optional
        List of directed edges. These should be hybrid edges only.
        Can be:
        - (u, v) tuples (key auto-generated)
        - (u, v, key) tuples (explicit key)
        - Dict with 'u', 'v' keys and optional 'key' and edge attributes
        
        Edge attributes can be set via dict format. Suggested attributes:
        - branch_length (float): Branch length for the edge
        - bootstrap (float): Bootstrap support value (must be in [0.0, 1.0])
        - gamma (float): For hybrid edges only - inheritance probability (must be in [0.0, 1.0]).
          **Gamma values can only be set on edges that point into hybrid nodes** (hybrid edges).
          Attempting to set gamma on a non-hybrid edge will raise a ValueError during validation.
          If ANY gamma is specified for edges entering a hybrid node,
          then ALL edges entering that hybrid node must have gamma values summing to 1.0.
        
        Must be provided (can be empty list for empty network, but a warning will be raised).
        By default None.
    undirected_edges : List[Union[Tuple[T, T], Tuple[T, T, int], Dict[str, Any]]], optional
        List of undirected edges. These should be non-hybrid edges only.
        Can be:
        - (u, v) tuples (key auto-generated)
        - (u, v, key) tuples (explicit key)
        - Dict with 'u', 'v' keys and optional 'key' and edge attributes
        
        Edge attributes can be set via dict format. Same suggested attributes as directed edges.
        By default None.
    taxa : Optional[Dict[T, str] | List[Tuple[T, str]]], optional
        Taxon labels for leaves. Can be:
        - Dict mapping leaf node IDs to taxon labels: {leaf_id: "taxon"}
        - List of (leaf_id, taxon_label) tuples
        The IDs in this mapping must be leaves (no outgoing directed edges).
        Not all leaves need to be covered - uncovered leaves get auto-generated labels.
        By default None.
    internal_node_labels : Optional[Dict[T, str] | List[Tuple[T, str]]], optional
        Labels for internal nodes (optional). Can be:
        - Dict mapping node IDs to labels: {node_id: "label"}
        - List of (node_id, label) tuples
        Only internal nodes (non-leaves) can be labeled via this parameter.
        Leaves must be labeled via the 'taxa' parameter. By default None.
    
    Notes
    -----
    This class is immutable after initialization. To create a network,
    build it using MixedMultiGraph and then create a SemiDirectedPhyNetwork from it,
    or initialize it with edges and taxa.
    
    The `validate()` method ensures that:
    - All hybrid edges are directed
    - All non-hybrid edges are undirected
    
    Examples
    --------
    >>> # Initialize with taxa mapping
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2)],
    ...     taxa={1: "A", 2: "B"}
    ... )
    >>> net.taxa
    {'A', 'B'}
    >>> # Network with hybrid node
    >>> net2 = SemiDirectedPhyNetwork(
    ...     directed_edges=[(5, 4), (6, 4)],  # Hybrid edges (must be directed)
    ...     undirected_edges=[(4, 1)],  # Tree edge (must be undirected)
    ...     taxa={1: "A"}
    ... )
    >>> net2.hybrid_nodes
    [4]
    """
    
    def __init__(
        self,
        directed_edges: Optional[List[Union[Tuple[T, T], Tuple[T, T, int], Dict[str, Any]]]] = None,
        undirected_edges: Optional[List[Union[Tuple[T, T], Tuple[T, T, int], Dict[str, Any]]]] = None,
        taxa: Optional[Dict[T, str] | List[Tuple[T, str]]] = None,
        internal_node_labels: Optional[Dict[T, str] | List[Tuple[T, str]]] = None,
    ) -> None:
        """
        Initialize a semi-directed phylogenetic network.
        
        Parameters
        ----------
        directed_edges : Optional[List[Union[Tuple[T, T], Tuple[T, T, int], Dict[str, Any]]]], optional
            List of directed edges. Must be provided (can be empty list for empty network,
            but a warning will be raised). By default None.
        undirected_edges : Optional[List[Union[Tuple[T, T], Tuple[T, T, int], Dict[str, Any]]]], optional
            List of undirected edges. Must be provided (can be empty list for empty network).
            By default None.
        taxa : Optional[Dict[T, str] | List[Tuple[T, str]]], optional
            Taxon labels for leaves. Can be:
            - Dict mapping leaf IDs to taxon labels: {leaf_id: "taxon"}
            - List of (leaf_id, taxon_label) tuples
            The IDs must be leaves (no outgoing directed edges). Not all leaves need to be
            covered - uncovered leaves get auto-generated labels. By default None.
        internal_node_labels : Optional[Dict[T, str] | List[Tuple[T, str]]], optional
            Labels for internal nodes (optional). Can be:
            - Dict mapping node IDs to labels: {node_id: "label"}
            - List of (node_id, label) tuples
            Only internal nodes (non-leaves) can be labeled via this parameter.
            Leaves must be labeled via 'taxa'. By default None.
        
        Examples
        --------
        >>> net = SemiDirectedPhyNetwork(
        ...     undirected_edges=[(3, 1), (3, 2)],
        ...     taxa={1: "A", 2: "B"}
        ... )
        >>> net = SemiDirectedPhyNetwork(
        ...     directed_edges=[(5, 4)],
        ...     undirected_edges=[(4, 1)],
        ...     taxa={1: "A"}
        ... )
        """
        # Call parent constructor
        super().__init__(
            directed_edges=directed_edges,
            undirected_edges=undirected_edges,
            taxa=taxa,
            internal_node_labels=internal_node_labels
        )
    
    def validate(self) -> bool:
        """
        Validate the network structure and edge attributes.
        
        Checks:
        1. Network is connected
        2. All hybrid edges are directed
        3. All non-hybrid edges are undirected
        
        Returns
        -------
        bool
            True if the network is valid, False otherwise.
        
        Raises
        ------
        ValueError
            If any structural or edge attribute constraints are violated.
        
        Notes
        -----
        This method performs comprehensive validation of the network structure
        and edge attributes. Empty networks (no nodes) are considered valid.
        """
        # First call parent validation (checks connectivity)
        if not super().validate():
            return False
        
        # Empty networks are valid
        if self.number_of_nodes() == 0:
            return True
        
        # Identify potential hybrid nodes (nodes with total in-degree >= 2)
        # Total in-degree includes both directed and undirected edges
        potential_hybrid_nodes = set()
        for v in self._graph.nodes:
            # Count incoming edges (directed + undirected)
            directed_in = self._graph.indegree(v)
            undirected_in = sum(1 for _ in self._graph.incident_undirected_edges(v))
            total_in = directed_in + undirected_in
            if total_in >= 2:
                potential_hybrid_nodes.add(v)
        
        # For each potential hybrid node, all incoming edges must be directed
        for v in potential_hybrid_nodes:
            # Check all undirected edges incident to this node
            for u, w, key in self._graph._undirected.edges(keys=True):
                if w == v:  # Incoming undirected edge
                    raise ValueError(
                        f"Node {v} has in-degree >= 2, so it is a hybrid node. "
                        f"Undirected edge ({u}, {v}, key={key}) entering hybrid node {v} "
                        f"must be directed. In a semi-directed network, all hybrid edges must be directed."
                    )
        
        # Get hybrid edges (directed edges pointing into hybrid nodes)
        hybrid_edges_set = set(self.hybrid_edges)
        
        # Check that all directed edges are hybrid edges
        for u, v, key in self._graph._directed.edges(keys=True):
            if (u, v) not in hybrid_edges_set:
                raise ValueError(
                    f"Directed edge ({u}, {v}, key={key}) is not a hybrid edge. "
                    f"In a semi-directed network, only hybrid edges can be directed."
                )
        
        # Check that all undirected edges are not hybrid edges
        for u, v, key in self._graph._undirected.edges(keys=True):
            if (u, v) in hybrid_edges_set:
                raise ValueError(
                    f"Undirected edge ({u}, {v}, key={key}) is a hybrid edge. "
                    f"In a semi-directed network, hybrid edges must be directed."
                )
        
        return True
    
    @cached_property
    def hybrid_nodes(self) -> List[T]:
        """
        Return a list of all hybrid nodes.
        
        In a semi-directed network, a hybrid node is a node with in-degree >= 2
        from directed edges (regardless of out-degree, since outgoing edges may be undirected).
        
        Returns
        -------
        List[T]
            List of hybrid node identifiers.
        
        Examples
        --------
        >>> net = SemiDirectedPhyNetwork(
        ...     directed_edges=[(5, 4), (6, 4)],
        ...     undirected_edges=[(4, 1)],
        ...     taxa={1: "A"}
        ... )
        >>> net.hybrid_nodes
        [4]
        """
        return [
            v for v in self._graph.nodes
            if self._graph.indegree(v) >= 2
        ]
    
    @cached_property
    def hybrid_edges(self) -> List[Tuple[T, T]]:
        """
        Return a list of all hybrid edges (directed edges pointing into hybrid nodes).
        
        Returns
        -------
        List[Tuple[T, T]]
            List of (source, target) tuples for hybrid edges.
        
        Examples
        --------
        >>> net = SemiDirectedPhyNetwork(
        ...     directed_edges=[(3, 2), (4, 2)],
        ...     taxa={2: "A"}
        ... )
        >>> net.hybrid_edges
        [(3, 2), (4, 2)]
        """
        res = []
        for v in self.hybrid_nodes:
            for p in self._graph._directed.predecessors(v):
                res.append((p, v))
        return res
    
    @cached_property
    def tree_edges(self) -> List[Tuple[T, T]]:
        """
        Return a list of all tree edges (undirected edges).
        
        In a semi-directed network, tree edges are all undirected edges.
        
        Returns
        -------
        List[Tuple[T, T]]
            List of (source, target) tuples for tree edges.
        
        Examples
        --------
        >>> net = SemiDirectedPhyNetwork(
        ...     undirected_edges=[(3, 1), (3, 2)],
        ...     taxa={1: "A", 2: "B"}
        ... )
        >>> len(net.tree_edges)
        2
        """
        # In semi-directed networks, tree edges are all undirected edges
        return list(self._graph._undirected.edges())
    
    def __repr__(self) -> str:
        """
        Return string representation of the network.
        
        Returns
        -------
        str
            String representation showing nodes, edges, and taxa count.
        """
        sorted_taxa = sorted(self.taxa)
        n_taxa = len(sorted_taxa)
        
        # Truncate taxon list at 10, add dots if longer
        if n_taxa <= 10:
            taxa_list_str = ", ".join(sorted_taxa)
        else:
            taxa_list_str = ", ".join(sorted_taxa[:10]) + ", ..."
        
        return (
            f"SemiDirectedPhyNetwork(nodes={self.number_of_nodes()}, "
            f"directed_edges={self._graph._directed.number_of_edges()}, "
            f"undirected_edges={self._graph._undirected.number_of_edges()}, "
            f"taxa={n_taxa}, "
            f"taxa_list=[{taxa_list_str}])"
        )


def random_semi_directed_network(n_leaves: int, seed: Optional[int] = None) -> SemiDirectedPhyNetwork:
    """
    Generate a random semi-directed network.
    
    Parameters
    ----------
    n_leaves : int
        Number of leaves in the network
    seed : Optional[int], optional
        Random seed, by default None
    
    Returns
    -------
    SemiDirectedPhyNetwork
        A random semi-directed network (placeholder - returns empty network)
    
    Notes
    -----
    This is a placeholder function. Implement actual random generation logic here.
    """
    # Placeholder implementation
    return SemiDirectedPhyNetwork(directed_edges=[], undirected_edges=[])
