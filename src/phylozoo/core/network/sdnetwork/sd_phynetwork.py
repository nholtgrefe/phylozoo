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
    
    A SemiDirectedPhyNetwork is a mixed multi-graph that can be obtained from
    a directed phylogenetic LSA (Least Stable Ancestor) network by:
    1. Suppressing the outdegree-2 root node (if it exists)
    2. Undirecting all non-hybrid edges
    
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
        
        Checks all constraints from MixedPhyNetwork plus additional semi-directed network constraints:
        - All hybrid edges are directed
        - All non-hybrid edges are undirected
        
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
        This method calls the parent validation (which checks connectivity, internal node degrees,
        indegree constraints, bootstrap, and gamma) and then performs additional semi-directed
        network specific checks.
        """
        # Call parent validation (checks connectivity, internal node degrees, indegree constraints,
        # bootstrap, and gamma constraints)
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            if not super().validate():
                return False
        
        # Empty networks are valid
        if self.number_of_nodes() == 0:
            return True
        
        # Additional semi-directed network constraints
        # TODO: Implement additional validation logic for semi-directed networks
        # - All hybrid edges are directed
        # - All non-hybrid edges are undirected
        
        return True
    
    def copy(self) -> 'SemiDirectedPhyNetwork':
        """
        Create a copy of the network.
        
        Returns a shallow copy of the network. Cached properties are not
        copied but will be recomputed on first access.
        
        Returns
        -------
        SemiDirectedPhyNetwork
            A copy of the network.
        
        Examples
        --------
        >>> net = SemiDirectedPhyNetwork(undirected_edges=[(3, 1)], taxa={1: "A"})
        >>> net2 = net.copy()
        >>> net.number_of_nodes()
        2
        >>> net2.number_of_nodes()
        2
        >>> isinstance(net2, SemiDirectedPhyNetwork)
        True
        """
        # Create new network by copying internal structures
        new_net = SemiDirectedPhyNetwork.__new__(SemiDirectedPhyNetwork)
        new_net._graph = self._graph.copy()
        new_net._node_to_label = self._node_to_label.copy()
        new_net._label_to_node = self._label_to_node.copy()
        return new_net
    
    @cached_property
    def level(self) -> int:
        """
        Return the level of the network.
        
        Placeholder: To be implemented.
        
        Returns
        -------
        int
            Level of the network (placeholder, returns 0 for now).
        """
        # TODO: Implement level calculation
        return 0
    
    def __repr__(self) -> str:
        """
        Return string representation of the network.
        
        Returns
        -------
        str
            String representation showing nodes, edges, level, taxa count, and taxon list.
        
        Examples
        --------
        >>> net = SemiDirectedPhyNetwork(undirected_edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"})
        >>> repr(net)
        'SemiDirectedPhyNetwork(nodes=3, edges=2, level=0, taxa=2, taxa_list=[A, B])'
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
            f"edges={self.number_of_edges()}, "
            f"level={self.level}, "
            f"taxa={n_taxa}, "
            f"taxa_list=[{taxa_list_str}])"
        )


