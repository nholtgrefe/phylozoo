"""
Semi-directed network module.

This module provides classes and functions for working with semi-directed phylogenetic networks.
"""

import threading
import warnings
from functools import cached_property
from typing import Any, Dict, Iterator, List, Optional, Set, Tuple, TypeVar, Union

from ...primitives.m_multigraph.mm_operations import (
    is_connected,
    orient_away_from_vertex,
    source_components,
)
from .m_phynetwork import MixedPhyNetwork

T = TypeVar('T')

# Thread-local flag to prevent infinite recursion in _validate_semidir_constraint
_validating_semidir = threading.local()


class SemiDirectedPhyNetwork(MixedPhyNetwork):
    """
    A semi-directed phylogenetic network.
    
    A SemiDirectedPhyNetwork is a mixed multigraph (with both directed and undirected edges)
    representing a phylogenetic network structure. It is a subclass of MixedPhyNetwork with
    additional constraints: exactly the non-hybrid edges are undirected, and all hybrid edges
    remain directed. This ensures the network does not have undirected cycles.
    
    A SemiDirectedPhyNetwork can be obtained from a directed phylogenetic LSA (Least Stable
    Ancestor) network by:
    
    1. Undirecting all non-hybrid edges
    2. Suppressing degree-2 nodes
    
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
        
        The validation proceeds via round-trip conversion:
        
        1. Compute source components of the mixed graph; there must be exactly one.
        2. Pick one non-leaf node from that source component as a root.
        3. Orient all edges away from the root (using ``orient_away_from_vertex``).
        4. Build a ``DirectedPhyNetwork`` from the oriented graph.
        5. Check if the chosen root is the LSA node. If not, the SD network is invalid.
        6. Convert it back to a semi-directed network with ``to_sd_network``.
        7. Assert that the resulting semi-directed network matches the original
           (same nodes, same directed edges, same undirected edges).
           
        Note: Any non-leaf node in the source component should yield the same result
        (becomes the LSA when oriented). We only need to check one.
        
        Returns
        -------
        bool
            True if the semi-directed constraints are satisfied.
        
        Raises
        ------
        ValueError
            If the network does not satisfy semi-directed network constraints.
        
        Notes
        -----
        This method validates that:
        - All hybrid edges are directed
        - All non-hybrid edges are undirected
        - The network can be oriented and converted back to the same structure
        """
        # Skip validation if we're already inside _validate_semidir_constraint (prevents infinite recursion)
        if getattr(_validating_semidir, 'active', False):
            return True
        
        # Handle trivial cases: 0 or 1 node networks are always valid
        n_nodes = self.number_of_nodes()
        if n_nodes <= 1:
            return True
        
        # Handle two-node case: only valid if it's an undirected edge with two leaves
        if n_nodes == 2:
            nodes_list = list(self._graph.nodes())
            u, v = nodes_list[0], nodes_list[1]
            
            # Check if both are leaves (outdegree 0)
            if self._graph.outdegree(u) != 0 or self._graph.outdegree(v) != 0:
                raise ValueError(
                    "Two-node semi-directed network must have two leaves connected by an undirected edge"
                )
            
            # Check if there's exactly one undirected edge between them
            undirected_count = len(list(self._graph._undirected.edges(u, v, keys=True)))
            if undirected_count != 1:
                raise ValueError(
                    "Two-node semi-directed network must have exactly one undirected edge"
                )
            
            # Check no directed edges
            if len(list(self._graph._directed.edges(u, v, keys=True))) > 0:
                raise ValueError(
                    "Two-node semi-directed network cannot have directed edges"
                )
            
            return True
        
        # Set flag to prevent infinite recursion
        _validating_semidir.active = True
        try:
            # Step 1: find source components (do this before imports to fail fast)
            components = source_components(self._graph)
            if len(components) != 1:
                raise ValueError(
                    f"Semi-directed network must have exactly one source component; found {len(components)}"
                )
            nodes_in_component, _, _ = components[0]
            if not nodes_in_component:
                raise ValueError("Source component is empty")
            
            # Local imports to avoid circular dependencies
            from ...network.dnetwork.operations import to_sd_network
            from ...network.dnetwork.d_phynetwork import DirectedPhyNetwork
            from ...network.dnetwork.classifications import is_LSA_network
            
            # Step 2: Pick a non-leaf vertex r from the source component
            # We cannot pick a leaf node as root, as it won't become the LSA when oriented
            # Pick any internal node (non-leaf based on the network's definition)
            internal_in_component = [node for node in nodes_in_component if node in self.internal_nodes]
            if not internal_in_component:
                raise ValueError(
                    "Semi-directed network constraint validation failed: "
                    "source component contains only leaf nodes."
                )
            root = internal_in_component[0]
            
            # Step 3: Orient away from r
            # If not fully directed or not valid, raise error immediately
            oriented_dm = orient_away_from_vertex(self._graph, root)
            
            # Step 4: Build a DirectedPhyNetwork from the oriented graph
            directed_edges: List[Dict[str, Any]] = []
            for u, v, key, data in oriented_dm.edges(keys=True, data=True):
                edge_dict: Dict[str, Any] = {"u": u, "v": v}
                if key != 0:
                    edge_dict["key"] = key
                if data:
                    edge_dict.update(data)
                directed_edges.append(edge_dict)
            
            # Build taxa mapping based on the oriented graph (leaf = outdegree 0)
            oriented_leaves: Set[Any] = {
                node for node in oriented_dm.nodes() if oriented_dm.outdegree(node) == 0
            }
            taxa_mapping: Dict[Any, str] = {
                leaf: self.get_label(leaf) or str(leaf)
                for leaf in oriented_leaves
            }
            
            # Build DirectedPhyNetwork - if not valid, raise error immediately
            try:
                d_network = DirectedPhyNetwork(
                    edges=directed_edges,
                    taxa=taxa_mapping if taxa_mapping else None,
                    internal_node_labels=None,
                )
            except ValueError as e:
                raise ValueError(
                    f"Semi-directed network constraint validation failed: "
                    f"oriented network is not a valid DirectedPhyNetwork: {e}"
                )
            
            # Step 5: Check if r is the LSA of the DirectedPhyNetwork
            # If not, raise error immediately
            if not is_LSA_network(d_network):
                raise ValueError(
                    "Semi-directed network constraint validation failed: "
                    "the chosen root is not the LSA node of the oriented network."
                )
            
            # Step 6: Convert back to semi-directed using to_sd_network
            # Since the network is already an LSA network, to_sd_network won't need to convert it
            # to an LSA network again.
            try:
                sd_roundtrip = to_sd_network(d_network)
            except ValueError as e:
                # If conversion fails, the SD network is invalid
                raise ValueError(
                    f"Semi-directed network constraint validation failed: "
                    f"conversion back to semi-directed failed: {e}"
                )
            
            # Step 7: Compare structure only (nodes, directed edges, undirected edges)
            # Compare nodes
            if set(self._graph.nodes()) != set(sd_roundtrip._graph.nodes()):
                raise ValueError(
                    "Semi-directed network constraint validation failed: "
                    "nodes do not match after round-trip conversion"
                )
            
            # Compare directed edges (only structure, not attributes)
            orig_directed = {
                (u, v, key if key is not None else 0)
                for u, v, key in self._graph.directed_edges_iter(keys=True)
            }
            roundtrip_directed = {
                (u, v, key if key is not None else 0)
                for u, v, key in sd_roundtrip._graph.directed_edges_iter(keys=True)
            }
            if orig_directed != roundtrip_directed:
                raise ValueError(
                    "Semi-directed network constraint validation failed: "
                    "directed edges do not match after round-trip conversion"
                )
            
            # Compare undirected edges (only structure, not attributes)
            orig_undirected = {
                (min(u, v), max(u, v), key if key is not None else 0)
                for u, v, key in self._graph.undirected_edges_iter(keys=True)
            }
            roundtrip_undirected = {
                (min(u, v), max(u, v), key if key is not None else 0)
                for u, v, key in sd_roundtrip._graph.undirected_edges_iter(keys=True)
            }
            if orig_undirected != roundtrip_undirected:
                raise ValueError(
                    "Semi-directed network constraint validation failed: "
                    "undirected edges do not match after round-trip conversion"
                )
            
            # All checks passed
            return True
            
        finally:
            # Clear flag
            _validating_semidir.active = False
    
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


