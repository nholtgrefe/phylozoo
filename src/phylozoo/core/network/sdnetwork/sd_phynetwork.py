"""
Semi-directed network module.

This module provides classes and functions for working with semi-directed phylogenetic networks.
"""

import warnings
from functools import cached_property
from typing import Any, Dict, Iterator, List, Optional, Set, Tuple, TypeVar, Union

from ...primitives.m_multigraph.mm_operations import is_connected
from .m_phynetwork import MixedPhyNetwork

T = TypeVar('T')


class SemiDirectedPhyNetwork(MixedPhyNetwork):
    """
    A semi-directed phylogenetic network.
    
    A SemiDirectedPhyNetwork is a mixed multigraph (with both directed and undirected edges)
    representing a phylogenetic network structure. It is a subclass of MixedPhyNetwork with
    additional constraints: exactly the non-hybrid edges are undirected, and all hybrid edges
    remain directed. This ensures the network does not have undirected cycles.
    
    A SemiDirectedPhyNetwork can be obtained from a directed phylogenetic LSA (Least Stable
    Ancestor) network by:
    
    1. Suppressing the outdegree-2 root node (if it exists)
    2. Undirecting all non-hybrid edges
    
    The network consists of:
    
    - **Leaf nodes**: Nodes with no outgoing directed edges, each with a taxon label
    - **Tree nodes**: Internal nodes with in-degree 0 and total degree >= 3
    - **Hybrid nodes**: Internal nodes with in-degree >= 2 and total degree = in-degree + 1
    
    The `validate()` method checks whether the network adheres to structural constraints
    and semi-directed network constraints upon initialization.
    
    This class uses composition with MixedMultiGraph for graph structure and adds
    phylogenetic-specific features like taxa labels, node labels, and network topology methods.
    
    Leaves refer to node IDs, and taxa refer to the labels of leaves. Leaves must always
    have labels (taxa). Internal nodes may or may not be labeled. Node IDs are separate
    from labels, allowing flexible node identification.
    
    Each label must be unique across all nodes. Attempting to use a duplicate label
    will raise a ValueError.
    
    Parameters
    ----------
    directed_edges : Optional[List[Union[Tuple[T, T], Tuple[T, T, int], Dict[str, Any]]]], optional
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
          
          Note: If you need more forgiving behavior (e.g., partial gamma values that
          don't sum to 1.0), consider using a different edge attribute name (e.g., 'gamma2')
          which will not be validated.
        
        These attributes are validated during initialization.
        
        Other edge attributes are also supported (any key-value pairs can be added),
        but the above three are suggested for phylogenetic networks.
        
        By default None.
    undirected_edges : Optional[List[Union[Tuple[T, T], Tuple[T, T, int], Dict[str, Any]]]], optional
        List of undirected edges. These should be non-hybrid edges only.
        Can be:
        - (u, v) tuples (key auto-generated)
        - (u, v, key) tuples (explicit key)
        - Dict with 'u', 'v' keys and optional 'key' and edge attributes
        
        Edge attributes can be set via dict format. Suggested attributes:
        - branch_length (float): Branch length for the edge
        - bootstrap (float): Bootstrap support value (must be in [0.0, 1.0])
        
        **Note: Undirected edges cannot have gamma values.**
        
        These attributes are validated during initialization.
        
        Other edge attributes are also supported (any key-value pairs can be added),
        but the above two are suggested for phylogenetic networks.
        
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
    
    Attributes
    ----------
    _graph : MixedMultiGraph[T]
        Internal graph structure using MixedMultiGraph.
        **Warning:** Do not modify directly. Use class methods instead.
    _node_to_label : Dict[T, str]
        Mapping from node IDs to labels. Only nodes with explicit labels are included.
        Leaves always have labels (taxa), but internal nodes may be unlabeled.
    _label_to_node : Dict[str, T]
        Reverse mapping from labels to node IDs (for quick lookup).
    
    Notes
    -----
    This class is immutable after initialization. To create a network,
    build it using MixedMultiGraph and then create a SemiDirectedPhyNetwork from it,
    or initialize it with edges and taxa.
    
    Examples
    --------
    >>> # Initialize with taxa mapping
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 4)],
    ...     taxa={1: "A", 2: "B", 4: "C"}
    ... )
    >>> net.taxa
    {'A', 'B', 'C'}
    >>> # Partial taxa mapping - uncovered leaves get auto-generated labels
    >>> net2 = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 4), (3, 5)],
    ...     taxa={1: "A"}
    ... )
    >>> net2.taxa  # 2, 4, and 5 are auto-labeled
    {'A', '2', '4', '5'}
    >>> # Network with branch lengths and bootstrap support
    >>> net3 = SemiDirectedPhyNetwork(
    ...     undirected_edges=[
    ...         {'u': 3, 'v': 1, 'branch_length': 0.5, 'bootstrap': 0.95},
    ...         {'u': 3, 'v': 2, 'branch_length': 0.3, 'bootstrap': 0.87},
    ...         (3, 4)
    ...     ],
    ...     taxa={1: "A", 2: "B", 4: "C"}
    ... )
    >>> net3.get_branch_length(3, 1)
    0.5
    >>> net3.get_bootstrap(3, 1)
    0.95
    >>> # Network with hybrid node and gamma values
    >>> net4 = SemiDirectedPhyNetwork(
    ...     directed_edges=[
    ...         {'u': 5, 'v': 4, 'gamma': 0.6},  # Hybrid edge
    ...         {'u': 6, 'v': 4, 'gamma': 0.4}  # Hybrid edge (Sum = 1.0)
    ...     ],
    ...     undirected_edges=[(4, 1), (4, 2), (4, 3)],  # Tree edges
    ...     taxa={1: "A", 2: "B", 3: "C"}
    ... )
    >>> net4.get_gamma(5, 4)
    0.6
    >>> net4.get_gamma(6, 4)
    0.4
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
    
    def _validate_semidir_constraint(self) -> bool:
        """
        Validate semi-directed network constraint.
        
        This method should be implemented to validate semi-directed network specific constraints.
        Currently always returns True as a placeholder.
        
        Returns
        -------
        bool
            Always returns True (placeholder).
        
        Notes
        -----
        This method should validate that:
        - All hybrid edges are directed
        - All non-hybrid edges are undirected
        """
        # TODO: Implement semi-directed network constraint validation
        return True
    
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
        This method performs the same validation checks as MixedPhyNetwork, but replaces
        the mixed network constraint check with a semi-directed network constraint check.
        """
        # Empty networks are valid
        if self.number_of_nodes() == 0:
            return True
        
        # 1. Check that network is connected (weakly connected)
        if not is_connected(self._graph):
            raise ValueError(
                "Network is not connected. All nodes must be in a single connected component."
            )
        
        # 2. Validate degree constraints
        self._validate_degree_constraints()
        
        # 3. Validate semi-directed network constraint (replaces mixed network constraint)
        self._validate_semidir_constraint()
        
        # 4. Validate bootstrap constraints
        self._validate_bootstrap_constraints()
        
        # 5. Validate gamma constraints
        self._validate_gamma_constraints()
        
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


