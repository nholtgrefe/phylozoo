"""
Mixed network module.

This module provides classes and functions for working with mixed phylogenetic networks.
"""

import math
import warnings
from functools import cached_property
from typing import Any, Dict, Iterator, List, Optional, Set, Tuple, TypeVar, Union

from ...primitives.m_multigraph import MixedMultiGraph
from ...primitives.m_multigraph.operations import is_connected

T = TypeVar('T')


class MixedPhyNetwork:
    """
    A mixed phylogenetic network.
    
    A MixedPhyNetwork is a mixed multigraph (with both directed and undirected edges)
    representing a phylogenetic network structure. This is an abstract network type which
    may have undirected cycles, used for canonical forms and to address unidentifiability
    issues. For semi-directed phylogenetic networks without undirected cycles, use the
    SemiDirectedPhyNetwork subclass.
    
    A MixedPhyNetwork can be obtained from a directed phylogenetic LSA (Least Stable
    Ancestor) network by:
    
    1. Undirecting all non-hybrid edges
    2. Optionally undirecting all hybrid edges for selected hybrid nodes
       (i.e., if one hybrid edge is undirected, all partner hybrid edges are undirected)
    3. Suppressing degree-2 nodes

    The network consists of:
    
    - **Leaf nodes**: Nodes with no outgoing directed edges, each with a taxon label
    - **Tree nodes**: Internal nodes with in-degree 0 and total degree >= 3
    - **Hybrid nodes**: Internal nodes with in-degree >= 2 and total degree = in-degree + 1
    
    The `validate()` method checks whether the network adheres to structural constraints
    upon initialization.
    
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
        List of directed edges. Can be:
        - (u, v) tuples (key auto-generated)
        - (u, v, key) tuples (explicit key)
        - Dict with 'u', 'v' keys and optional 'key' and edge attributes
        
        Edge attributes can be set via dict format. Suggested attributes:
        - branch_length (float): Branch length for the edge
        - bootstrap (float): Bootstrap support value (must be in [0.0, 1.0])
        - gamma (float): For hybrid edges only - inheritance probability (must be in [0.0, 1.0]).
          **Gamma values can only be set on directed edges that point into hybrid nodes**
          (hybrid edges). Undirected edges cannot have gamma values.
          Attempting to set gamma on an undirected or non-hybrid edge will raise a
          ValueError during validation. If ANY gamma is specified for edges entering a hybrid node,
          then ALL edges entering that hybrid node must have gamma values summing to 1.0.
          
          Note: If you need more forgiving behavior (e.g., partial gamma values that
          don't sum to 1.0), consider using a different edge attribute name (e.g., 'gamma2')
          which will not be validated.
        
        These attributes are validated during initialization.
        
        Other edge attributes are also supported (any key-value pairs can be added),
        but the above three are suggested for phylogenetic networks.
        
        By default None.
    undirected_edges : Optional[List[Union[Tuple[T, T], Tuple[T, T, int], Dict[str, Any]]]], optional
        List of undirected edges. Can be:
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
    build it using MixedMultiGraph and then create a MixedPhyNetwork from it,
    or initialize it with edges and taxa.
    
    Examples
    --------
    >>> # Initialize with taxa mapping
    >>> net = MixedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 4)],
    ...     taxa={1: "A", 2: "B", 4: "C"}
    ... )
    >>> net.taxa
    {'A', 'B', 'C'}
    >>> # Partial taxa mapping - uncovered leaves get auto-generated labels
    >>> net2 = MixedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 4), (3, 5)],
    ...     taxa={1: "A"}
    ... )
    >>> net2.taxa  # 2, 4, and 5 are auto-labeled
    {'A', '2', '4', '5'}
    >>> # Network with branch lengths and bootstrap support
    >>> net3 = MixedPhyNetwork(
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
    >>> net4 = MixedPhyNetwork(
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
        Initialize a mixed phylogenetic network.
        
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
        >>> net = MixedPhyNetwork(
        ...     undirected_edges=[(3, 1), (3, 2)],
        ...     taxa={1: "A", 2: "B"}
        ... )
        >>> net = MixedPhyNetwork(
        ...     directed_edges=[(5, 4)],
        ...     undirected_edges=[(4, 1)],
        ...     taxa={1: "A"}
        ... )
        """
        # Warn if both edge lists are empty
        directed_list = directed_edges if directed_edges is not None else []
        undirected_list = undirected_edges if undirected_edges is not None else []
        
        if not directed_list and not undirected_list:
            warnings.warn(
                "Initializing MixedPhyNetwork with empty edges lists. "
                "This creates an empty network.",
                UserWarning,
                stacklevel=2
            )
        
        self._graph: MixedMultiGraph[T] = MixedMultiGraph(
            directed_edges=directed_list,
            undirected_edges=undirected_list
        )
        self._node_to_label: Dict[T, str] = {}
        self._label_to_node: Dict[str, T] = {}
        
        self._initialize_labels(taxa, internal_node_labels)
        
        # Validate the network
        if not self.validate():
            raise ValueError("Network validation failed")
    
    def _initialize_labels(
        self,
        taxa: Optional[Dict[T, str] | List[Tuple[T, str]]],
        internal_node_labels: Optional[Dict[T, str] | List[Tuple[T, str]]],
    ) -> None:
        """
        Initialize node labels from taxa and internal node labels.
        
        This helper method processes taxa labels for leaves, auto-labels uncovered leaves,
        and processes internal node labels.
        
        Parameters
        ----------
        taxa : Optional[Dict[T, str] | List[Tuple[T, str]]]
            Taxon labels for leaves.
        internal_node_labels : Optional[Dict[T, str] | List[Tuple[T, str]]]
            Labels for internal nodes.
        
        Raises
        ------
        ValueError
            If taxa or internal_node_labels have invalid format, or if nodes in taxa
            are not leaves, or if nodes in internal_node_labels are leaves.
        """
        # Identify all leaves (nodes with no outgoing directed edges)
        all_leaves: Set[T] = self.leaves
        
        # Process taxa (taxon labels for leaves)
        covered_leaves: Set[T] = set()
        if taxa:
            if isinstance(taxa, dict):
                # Dict: leaf_id -> taxon_label
                for leaf_id, taxon_label in taxa.items():
                    # Verify this is indeed a leaf
                    if leaf_id not in all_leaves:
                        raise ValueError(
                            f"Node {leaf_id} in taxa mapping is not a leaf "
                            f"(has outgoing directed edges)"
                        )
                    self._set_label(leaf_id, taxon_label)
                    covered_leaves.add(leaf_id)
            elif isinstance(taxa, list):
                # List of (leaf_id, taxon_label) tuples
                for leaf_id, taxon_label in taxa:
                    # Verify this is indeed a leaf
                    if leaf_id not in all_leaves:
                        raise ValueError(
                            f"Node {leaf_id} in taxa mapping is not a leaf "
                            f"(has outgoing directed edges)"
                        )
                    self._set_label(leaf_id, taxon_label)
                    covered_leaves.add(leaf_id)
            else:
                raise ValueError(
                    "taxa must be a dict (leaf_id -> taxon_label) "
                    "or a list of (leaf_id, taxon_label) tuples"
                )
        
        # Auto-label uncovered leaves
        uncovered_leaves = all_leaves - covered_leaves
        
        for leaf_id in uncovered_leaves:
            # Auto-generate label based on node_id
            label = str(leaf_id)
            # If label already exists (e.g., from taxa mapping), make it unique
            if label in self._label_to_node:
                counter = 1
                unique_label = f"{label}_{counter}"
                while unique_label in self._label_to_node:
                    counter += 1
                    unique_label = f"{label}_{counter}"
                label = unique_label
            self._set_label(leaf_id, label)
        
        # Process internal node labels (optional)
        # Only set labels for internal nodes that are explicitly provided
        if internal_node_labels:
            if isinstance(internal_node_labels, dict):
                # Dict: node_id -> label
                for node_id, label in internal_node_labels.items():
                    if node_id in all_leaves:
                        raise ValueError(
                            f"Node {node_id} in internal_node_labels is a leaf. "
                            f"Use the 'taxa' parameter to set leaf labels (taxa)."
                        )
                    # Only set labels for internal nodes
                    self._set_label(node_id, label)
            elif isinstance(internal_node_labels, list):
                # List of (node_id, label) tuples
                for node_id, label in internal_node_labels:
                    if node_id in all_leaves:
                        raise ValueError(
                            f"Node {node_id} in internal_node_labels is a leaf. "
                            f"Use the 'taxa' parameter to set leaf labels (taxa)."
                        )
                    # Only set labels for internal nodes
                    self._set_label(node_id, label)
            else:
                raise ValueError(
                    "internal_node_labels must be a dict (node_id -> label) "
                    "or a list of (node_id, label) tuples"
                )
        
        # Note: Internal nodes without explicit labels remain unlabeled
        # Only leaves are guaranteed to have labels (taxa)
    
    def _validate_bootstrap_constraints(self) -> None:
        """
        Validate bootstrap values are within [0.0, 1.0].
        
        Raises
        ------
        ValueError
            If any edge has a bootstrap value outside [0.0, 1.0].
        """
        # Quick check: if no bootstrap values are set, skip validation
        has_bootstrap = False
        for u, v, key, data in self._graph._directed.edges(keys=True, data=True):
            if 'bootstrap' in data:
                has_bootstrap = True
                break
        if not has_bootstrap:
            for u, v, key, data in self._graph._undirected.edges(keys=True, data=True):
                if 'bootstrap' in data:
                    has_bootstrap = True
                    break
        if not has_bootstrap:
            return
        
        # Check directed edges
        for u, v, key, data in self._graph._directed.edges(keys=True, data=True):
            if 'bootstrap' in data:
                bootstrap = data['bootstrap']
                if not isinstance(bootstrap, (int, float)):
                    raise ValueError(
                        f"Bootstrap value on directed edge ({u}, {v}, key={key}) must be numeric, "
                        f"got {type(bootstrap).__name__}"
                    )
                if math.isnan(bootstrap) or bootstrap < 0.0 or bootstrap > 1.0:
                    raise ValueError(
                        f"Bootstrap value on directed edge ({u}, {v}, key={key}) is {bootstrap}, "
                        f"but must be in [0.0, 1.0]"
                    )
        
        # Check undirected edges
        for u, v, key, data in self._graph._undirected.edges(keys=True, data=True):
            if 'bootstrap' in data:
                bootstrap = data['bootstrap']
                if not isinstance(bootstrap, (int, float)):
                    raise ValueError(
                        f"Bootstrap value on undirected edge ({u}, {v}, key={key}) must be numeric, "
                        f"got {type(bootstrap).__name__}"
                    )
                if math.isnan(bootstrap) or bootstrap < 0.0 or bootstrap > 1.0:
                    raise ValueError(
                        f"Bootstrap value on undirected edge ({u}, {v}, key={key}) is {bootstrap}, "
                        f"but must be in [0.0, 1.0]"
                    )
    
    def _validate_gamma_constraints(self) -> None:
        """
        Validate gamma constraints for hybrid nodes.
        
        Gamma values can only be set on hybrid edges (directed edges pointing into hybrid nodes).
        For each hybrid node, if ANY incoming edge has a gamma value, then
        ALL incoming edges (including parallel edges) must have gamma values,
        and they must sum to exactly 1.0. If no gamma values are specified,
        no validation is performed.
        
        Raises
        ------
        ValueError
            If gamma constraints are violated, including gamma set on non-hybrid edges
            or undirected edges.
        """
        # Quick check: if no gamma values are set, skip validation
        has_gamma = False
        for u, v, key, data in self._graph._directed.edges(keys=True, data=True):
            if 'gamma' in data:
                has_gamma = True
                break
        if not has_gamma:
            for u, v, key, data in self._graph._undirected.edges(keys=True, data=True):
                if 'gamma' in data:
                    has_gamma = True
                    break
        if not has_gamma:
            return
        
        # First, check that gamma is only set on directed hybrid edges
        hybrid_edges_set = set(self.hybrid_edges)
        
        # Check directed edges
        for u, v, key, data in self._graph._directed.edges(keys=True, data=True):
            if 'gamma' in data:
                # Check if this edge is a hybrid edge
                if (u, v) not in hybrid_edges_set:
                    raise ValueError(
                        f"Gamma value can only be set on hybrid edges (edges pointing into "
                        f"hybrid nodes). Directed edge ({u}, {v}, key={key}) is not a hybrid edge."
                    )
        
        # Check that gamma is not set on undirected edges
        for u, v, key, data in self._graph._undirected.edges(keys=True, data=True):
            if 'gamma' in data:
                raise ValueError(
                    f"Gamma values cannot be set on undirected edges. "
                    f"Undirected edge ({u}, {v}, key={key}) cannot have gamma values."
                )
        
        # Then validate gamma constraints for hybrid nodes
        for hybrid_node in self.hybrid_nodes:
            gamma_values: List[float] = []
            incoming_edges: List[Tuple[T, T, int]] = []
            
            # Use incident_parent_edges to get all incoming directed edges (including parallel edges)
            for edge in self.incident_parent_edges(hybrid_node, keys=True, data=True):
                if len(edge) == 4:  # (u, v, key, data)
                    u, v, key, edge_data = edge
                    incoming_edges.append((u, v, key))
                    gamma = edge_data.get('gamma')
                    if gamma is not None:
                        # Validate gamma is numeric
                        if not isinstance(gamma, (int, float)):
                            raise ValueError(
                                f"Gamma value on edge ({u}, {v}, key={key}) entering hybrid node "
                                f"{hybrid_node} must be numeric, got {type(gamma).__name__}"
                            )
                        # Validate gamma is in [0.0, 1.0]
                        if gamma < 0.0 or gamma > 1.0:
                            raise ValueError(
                                f"Gamma value on edge ({u}, {v}, key={key}) entering hybrid node "
                                f"{hybrid_node} is {gamma}, but must be in [0.0, 1.0]"
                            )
                        gamma_values.append(gamma)
            
            # If any gamma values are set, ALL edges must have gamma values
            if len(gamma_values) > 0:
                # Check if all incoming edges have gamma values
                missing_edges: List[str] = []
                for u, v, key in incoming_edges:
                    gamma = self.get_edge_attribute(u, v, key, 'gamma')
                    if gamma is None:
                        missing_edges.append(f"({u}, {v}, key={key})")
                
                if missing_edges:
                    raise ValueError(
                        f"Hybrid node {hybrid_node} has some edges with gamma values "
                        f"but others without. If ANY gamma is specified, ALL incoming edges "
                        f"must have gamma values. Missing gamma on edges: {', '.join(missing_edges)}"
                    )
                
                # All gammas are present, check they sum to 1.0
                gamma_sum = sum(gamma_values)
                if abs(gamma_sum - 1.0) > 1e-10:
                    raise ValueError(
                        f"Hybrid node {hybrid_node} has gamma values that sum to {gamma_sum}, "
                        f"but must sum to exactly 1.0"
                    )
    
    def _validate_mixednetwork_constraint(self) -> bool:
        """
        Validate mixed network constraint.
        
        This is a placeholder for future mixed network specific validation.
        Currently always returns True but issues a warning that additional
        validation checks may be added later.
        
        Returns
        -------
        bool
            Always returns True.
        
        Warns
        -----
        UserWarning
            Always issued to indicate that additional validation checks may be added later.
        """
        warnings.warn(
            "Validity is not fully checked for MixedPhyNetworks. "
            "Additional validation checks may be added later. "
            "For validated networks, use SemiDirectedPhyNetwork.",
            UserWarning,
            stacklevel=3
        )
        return True
    
    def _validate_degree_constraints(self) -> None:
        """
        Validate degree constraints for nodes.
        
        Checks:
        1. All internal nodes have degree >= 3
        2. Each node has indegree either 0 or total_degree-1
        
        Raises
        ------
        ValueError
            If any degree constraints are violated.
        """
        # 1. Check that all internal nodes have degree >= 3
        internal_nodes = self.internal_nodes
        for node in internal_nodes:
            degree = self._graph.degree(node)
            if degree < 3:
                raise ValueError(
                    f"Internal node {node} has degree {degree}, but all internal nodes "
                    f"must have degree >= 3."
                )
        
        # 2. Check that each node has indegree either 0 or total_degree-1
        # This constraint applies to all nodes (including leaves)
        for node in self._graph.nodes:
            indegree = self._graph.indegree(node)
            total_degree = self._graph.degree(node)
            if indegree != 0 and indegree != total_degree - 1:
                raise ValueError(
                    f"Node {node} has indegree {indegree} and total degree {total_degree}. "
                    f"Each node must have indegree either 0 or total_degree-1."
                )
    
    def validate(self) -> bool:
        """
        Validate the network structure and edge attributes.
        
        Checks:
        1. Network is connected (weakly connected)
        2. All internal nodes have degree >= 3
        3. Each node has indegree either 0 or undirected_degree-1
        4. Mixed network constraint (issues warning that additional checks may be added)
        5. Bootstrap values are in [0.0, 1.0]
        6. Gamma constraints (gamma only on hybrid edges, sum to 1.0 if specified)
        
        Returns
        -------
        bool
            True if the network is valid, False otherwise.
        
        Raises
        ------
        ValueError
            If any structural or edge attribute constraints are violated.
        
        Warns
        -----
        UserWarning
            Always issued to indicate that additional validation checks may be added later.
        
        Notes
        -----
        This method performs validation checks but issues a warning that additional
        checks may be added later. For fully validated networks with additional
        constraints, use SemiDirectedPhyNetwork.
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
        
        # 3. Validate mixed network constraint (issues warning)
        self._validate_mixednetwork_constraint()
        
        # 4. Validate bootstrap constraints
        self._validate_bootstrap_constraints()
        
        # 5. Validate gamma constraints
        self._validate_gamma_constraints()
        
        return True
    
    # ========== Label Operations ==========
    
    def get_label(self, node_id: T) -> Optional[str]:
        """
        Get the label for a node.
        
        Parameters
        ----------
        node_id : T
            Node identifier.
        
        Returns
        -------
        Optional[str]
            Label for the node. Returns None if node has no label.
            Leaves always have labels (taxa), but internal nodes may be unlabeled.
        
        Examples
        --------
        >>> net = MixedPhyNetwork(undirected_edges=[(3, 1)], taxa={1: "A"})
        >>> net.get_label(1)
        'A'
        >>> net.get_label(3) is None
        True
        """
        return self._node_to_label.get(node_id)
    
    def get_node_id(self, label: str) -> Optional[T]:
        """
        Get the node ID for a label.
        
        Parameters
        ----------
        label : str
            Node label.
        
        Returns
        -------
        Optional[T]
            Node ID if found, None otherwise.
        
        Examples
        --------
        >>> net = MixedPhyNetwork(undirected_edges=[(3, 1)], taxa={1: "A"})
        >>> net.get_node_id("A")
        1
        """
        return self._label_to_node.get(label)
    
    def _set_label(self, node_id: T, label: str) -> None:
        """
        Set or update the label for a node (private method for initialization only).
        
        Each label must be unique across all nodes. If the label is already
        used by a different node, a ValueError is raised.
        
        Parameters
        ----------
        node_id : T
            Node identifier.
        label : str
            Label to set. Must be unique across all nodes.
        
        Raises
        ------
        ValueError
            If node_id is not in the network, or if the label is already
            used by a different node.
        """
        if node_id not in self._graph:
            raise ValueError(f"Node {node_id} is not in the network")
        
        old_label = self._node_to_label.get(node_id)
        
        # Check if label is already used by a different node
        # Allow same node to update its own label
        if label in self._label_to_node:
            existing_node = self._label_to_node[label]
            if existing_node != node_id:
                raise ValueError(
                    f"Label '{label}' is already used by node {existing_node}. "
                    f"Each label must be unique."
                )
            # Same node, same label - nothing to do
            return
        
        # Update mappings
        self._node_to_label[node_id] = label
        if old_label and old_label in self._label_to_node and old_label != label:
            del self._label_to_node[old_label]
        self._label_to_node[label] = node_id
    
    # ========== Edge Attribute Access (Read-Only) ==========
    
    def get_edge_attribute(
        self,
        u: T,
        v: T,
        key: Optional[int] = None,
        attr: str = 'branch_length'
    ) -> Optional[Any]:
        """
        Get an edge attribute value (read-only).
        
        Parameters
        ----------
        u, v : T
            Edge endpoints.
        key : Optional[int], optional
            Edge key for parallel edges. If None and multiple parallel edges exist,
            raises ValueError. Must specify key when parallel edges exist.
        attr : str, default 'branch_length'
            Attribute name. Suggested attributes: 'branch_length', 'bootstrap', 'gamma'.
            Any edge attribute can be accessed.
        
        Returns
        -------
        Optional[Any]
            Attribute value, or None if not set.
        
        Raises
        ------
        ValueError
            If key is None and multiple parallel edges exist between u and v.
        
        Notes
        -----
        Networks will not have undirected and directed edges with the same endpoints,
        so this method searches both directed and undirected edges (directed first).
        """
        # Try directed edges first
        if self._graph._directed.has_edge(u, v, key=key):
            # Get all directed edges between u and v
            edges_data = self._graph._directed[u].get(v, {})
            if not edges_data:
                return None
            
            num_edges = len(edges_data)
            if num_edges == 0:
                return None
            elif num_edges == 1:
                first_key = next(iter(edges_data))
                return edges_data[first_key].get(attr)
            else:
                if key is None:
                    raise ValueError(
                        f"Multiple parallel directed edges exist between {u} and {v}. "
                        f"Must specify 'key' parameter to get attribute from a specific edge."
                    )
                if key not in edges_data:
                    return None
                return edges_data[key].get(attr)
        
        # Try undirected edges
        if self._graph._undirected.has_edge(u, v, key=key):
            # Get all undirected edges between u and v
            edges_data = self._graph._undirected[u].get(v, {})
            if not edges_data:
                return None
            
            num_edges = len(edges_data)
            if num_edges == 0:
                return None
            elif num_edges == 1:
                first_key = next(iter(edges_data))
                return edges_data[first_key].get(attr)
            else:
                if key is None:
                    raise ValueError(
                        f"Multiple parallel undirected edges exist between {u} and {v}. "
                        f"Must specify 'key' parameter to get attribute from a specific edge."
                    )
                if key not in edges_data:
                    return None
                return edges_data[key].get(attr)
        
        return None
    
    def get_branch_length(
        self,
        u: T,
        v: T,
        key: Optional[int] = None
    ) -> Optional[float]:
        """
        Get branch length for an edge.
        
        Parameters
        ----------
        u, v : T
            Edge endpoints.
        key : Optional[int], optional
            Edge key for parallel edges. Required if multiple parallel edges exist.
        
        Returns
        -------
        Optional[float]
            Branch length, or None if not set.
        """
        return self.get_edge_attribute(u, v, key, 'branch_length')
    
    def get_bootstrap(
        self,
        u: T,
        v: T,
        key: Optional[int] = None
    ) -> Optional[float]:
        """
        Get bootstrap support for an edge.
        
        Bootstrap values are typically in the range 0.0 to 1.0.
        
        Parameters
        ----------
        u, v : T
            Edge endpoints.
        key : Optional[int], optional
            Edge key for parallel edges. Required if multiple parallel edges exist.
        
        Returns
        -------
        Optional[float]
            Bootstrap support value (typically 0.0 to 1.0), or None if not set.
        """
        return self.get_edge_attribute(u, v, key, 'bootstrap')
    
    def get_gamma(
        self,
        u: T,
        v: T,
        key: Optional[int] = None
    ) -> Optional[float]:
        """
        Get gamma value for a hybrid edge.
        
        Gamma values can only be set on directed hybrid edges (directed edges pointing
        into hybrid nodes). Undirected edges cannot have gamma values.
        If ANY gamma value is specified for edges entering a hybrid node, then
        ALL edges entering that hybrid node must have gamma values, and they must
        sum to exactly 1.0.
        
        Parameters
        ----------
        u, v : T
            Edge endpoints (v must be a hybrid node, edge must be directed).
        key : Optional[int], optional
            Edge key for parallel edges. Required if multiple parallel edges exist.
        
        Returns
        -------
        Optional[float]
            Gamma value, or None if not set.
        
        Raises
        ------
        ValueError
            If the edge is undirected (gamma can only be on directed edges).
        """
        # Gamma is only on directed hybrid edges - check if edge is undirected
        if self._graph._undirected.has_edge(u, v, key=key):
            raise ValueError(
                f"Gamma values cannot be set on undirected edges. "
                f"Edge ({u}, {v}, key={key}) is undirected."
            )
        # Gamma is only on directed hybrid edges
        if not self._graph._directed.has_edge(u, v, key=key):
            return None
        return self.get_edge_attribute(u, v, key, 'gamma')
    
    # ========== Graph Query Operations (Delegated) ==========
    
    def number_of_nodes(self) -> int:
        """
        Return the number of nodes.
        
        Returns
        -------
        int
            Number of nodes.
        """
        return self._graph.number_of_nodes()
    
    def number_of_edges(self) -> int:
        """
        Return the number of edges.
        
        Returns
        -------
        int
            Number of edges (directed + undirected).
        """
        return self._graph.number_of_edges()
    
    def has_edge(
        self,
        u: T,
        v: T,
        key: Optional[int] = None,
        directed: Optional[bool] = None
    ) -> bool:
        """
        Check if edge exists.
        
        Parameters
        ----------
        u : T
            Source node.
        v : T
            Target node.
        key : Optional[int], optional
            Edge key. By default None.
        directed : Optional[bool], optional
            If True, only check directed edges. If False, only check undirected edges.
            If None, check both. By default None.
        
        Returns
        -------
        bool
            True if edge exists, False otherwise.
        """
        if directed is None:
            return self._graph._directed.has_edge(u, v, key=key) or \
                   self._graph._undirected.has_edge(u, v, key=key)
        elif directed:
            return self._graph._directed.has_edge(u, v, key=key)
        else:
            return self._graph._undirected.has_edge(u, v, key=key)
    
    def degree(self, v: T) -> int:
        """
        Return the total degree of node v.
        
        Parameters
        ----------
        v : T
            Node identifier.
        
        Returns
        -------
        int
            Total degree (undirected + in-degree + out-degree).
        """
        return self._graph.degree(v)
    
    def indegree(self, v: T) -> int:
        """
        Return the in-degree of node v (directed edges only).
        
        Parameters
        ----------
        v : T
            Node identifier.
        
        Returns
        -------
        int
            In-degree of v.
        """
        return self._graph.indegree(v)
    
    def outdegree(self, v: T) -> int:
        """
        Return the out-degree of node v (directed edges only).
        
        Parameters
        ----------
        v : T
            Node identifier.
        
        Returns
        -------
        int
            Out-degree of v.
        """
        return self._graph.outdegree(v)
    
    def undirected_degree(self, v: T) -> int:
        """
        Return the undirected degree of node v.
        
        Parameters
        ----------
        v : T
            Node identifier.
        
        Returns
        -------
        int
            Undirected degree of v.
        """
        return self._graph.undirected_degree(v)
    
    def incident_parent_edges(
        self,
        v: T,
        keys: bool = False,
        data: bool = False
    ) -> Iterator[Tuple[T, T] | Tuple[T, T, int] | Tuple[T, T, int, Dict[str, Any]]]:
        """
        Return an iterator over directed edges entering node v (from parent nodes).
        
        Parameters
        ----------
        v : T
            Node identifier.
        keys : bool, optional
            If True, return edge keys. By default False.
        data : bool, optional
            If True, return edge data. By default False.
        
        Returns
        -------
        Iterator
            Iterator over incoming directed edges as (u, v) or (u, v, key) or (u, v, key, data).
        """
        return self._graph.incident_parent_edges(v, keys=keys, data=data)
    
    def incident_child_edges(
        self,
        v: T,
        keys: bool = False,
        data: bool = False
    ) -> Iterator[Tuple[T, T] | Tuple[T, T, int] | Tuple[T, T, int, Dict[str, Any]]]:
        """
        Return an iterator over directed edges leaving node v (to child nodes).
        
        Parameters
        ----------
        v : T
            Node identifier.
        keys : bool, optional
            If True, return edge keys. By default False.
        data : bool, optional
            If True, return edge data. By default False.
        
        Returns
        -------
        Iterator
            Iterator over outgoing directed edges as (v, u) or (v, u, key) or (v, u, key, data).
        """
        return self._graph.incident_child_edges(v, keys=keys, data=data)
    
    def incident_undirected_edges(
        self,
        v: T,
        keys: bool = False,
        data: bool = False
    ) -> Iterator[Tuple[T, T] | Tuple[T, T, int] | Tuple[T, T, int, Dict[str, Any]]]:
        """
        Return an iterator over undirected edges incident to node v.
        
        Parameters
        ----------
        v : T
            Node identifier.
        keys : bool, optional
            If True, return edge keys. By default False.
        data : bool, optional
            If True, return edge data. By default False.
        
        Returns
        -------
        Iterator
            Iterator over undirected edges as (u, v) or (u, v, key) or (u, v, key, data).
        """
        return self._graph.incident_undirected_edges(v, keys=keys, data=data)
    
    def neighbors(self, v: T) -> Iterator[T]:
        """
        Return an iterator over neighbors of node v.
        
        Neighbors include nodes connected by both directed and undirected edges.

        Parameters
        ----------
        v : T
            Node identifier.
        
        Returns
        -------
        Iterator[T]
            Iterator over neighbors.
        
        Examples
        --------
        >>> net = MixedPhyNetwork(
        ...     directed_edges=[(1, 2)],
        ...     undirected_edges=[(2, 3)],
        ...     taxa={3: "A"}
        ... )
        >>> sorted(net.neighbors(2))
        [1, 3]
        """
        return self._graph.neighbors(v)
    
    # ========== Phylogenetic-Specific Methods ==========
    
    @cached_property
    def leaves(self) -> Set[T]:
        """
        Get the set of leaf node IDs (nodes with no degree-1).
        
        Returns
        -------
        Set[T]
            Set of leaf node identifiers. Returns a new set (which is mutable).
        """
        return {node for node in self._graph.nodes if self._graph.degree(node) == 1}
    
    @cached_property
    def taxa(self) -> Set[str]:
        """
        Get the set of taxon labels (labels of leaves).
        
        Returns
        -------
        Set[str]
            Set of taxon labels.
        
        Examples
        --------
        >>> net = MixedPhyNetwork(undirected_edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"})
        >>> net.taxa
        {'A', 'B'}
        """
        return {self._node_to_label[leaf] for leaf in self.leaves}
    
    @cached_property
    def hybrid_nodes(self) -> List[T]:
        """
        Return a list of all hybrid nodes.
        
        A hybrid node is a node with in-degree >= 2 and total degree = in-degree + 1.
        
        Returns
        -------
        List[T]
            List of hybrid node identifiers.
        
        Examples
        --------
        >>> net = MixedPhyNetwork(
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
            and self._graph.degree(v) == self._graph.indegree(v) + 1
        ]
    
    @cached_property
    def hybrid_edges(self) -> List[Tuple[T, T]]:
        """
        Return a list of all hybrid edges.
        
        Hybrid edges are all directed edges.
        
        Returns
        -------
        List[Tuple[T, T]]
            List of (source, target) tuples for hybrid edges.
        
        Examples
        --------
        >>> net = MixedPhyNetwork(
        ...     directed_edges=[(3, 2), (4, 2)],
        ...     undirected_edges=[(2, 1)],
        ...     taxa={1: "A"}
        ... )
        >>> net.hybrid_edges
        [(3, 2), (4, 2)]
        """
        return list(self._graph._directed.edges())
    
    @cached_property
    def tree_nodes(self) -> List[T]:
        """
        Return a list of all tree nodes.
        
        A tree node is a node with in-degree 0 and total degree >= 3.
        
        Returns
        -------
        List[T]
            List of tree node identifiers.
        
        Examples
        --------
        >>> net = MixedPhyNetwork(
        ...     undirected_edges=[(1, 2), (1, 3), (1, 4)],
        ...     taxa={2: "A", 3: "B", 4: "C"}
        ... )
        >>> net.tree_nodes
        [1]
        """
        return [
            v for v in self._graph.nodes
            if self._graph.indegree(v) == 0
            and self._graph.degree(v) >= 3
        ]
    
    @cached_property
    def internal_nodes(self) -> List[T]:
        """
        Return a list of all internal nodes.
        
        Internal nodes are all nodes that are not leaves.
        
        Returns
        -------
        List[T]
            List of internal node identifiers.
        
        Examples
        --------
        >>> net = MixedPhyNetwork(
        ...     directed_edges=[(5, 4), (6, 4)],
        ...     undirected_edges=[(4, 1)],
        ...     taxa={1: "A"}
        ... )
        >>> sorted(net.internal_nodes)
        [4, 5, 6]
        """
        leaves_set = set(self.leaves)
        return [
            v for v in self._graph.nodes
            if v not in leaves_set
        ]
    
    @cached_property
    def tree_edges(self) -> List[Tuple[T, T]]:
        """
        Return a list of all tree edges.
        
        Tree edges are simply all undirected edges.
        
        Returns
        -------
        List[Tuple[T, T]]
            List of (source, target) tuples for tree edges.
        
        Examples
        --------
        >>> net = MixedPhyNetwork(
        ...     directed_edges=[(5, 4), (6, 4)],
        ...     undirected_edges=[(4, 1)],
        ...     taxa={1: "A"}
        ... )
        >>> sorted(net.tree_edges)
        [(4, 1)]
        """
        return list(self._graph._undirected.edges())

    def __repr__(self) -> str:
        """
        Return string representation of the network.

        Returns
        -------
        str
            String representation showing nodes, edges, level, taxa count, and taxon list.
        
        Examples
        --------
        >>> net = MixedPhyNetwork(undirected_edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"})
        >>> repr(net)
        'MixedPhyNetwork(nodes=3, edges=2, taxa=2, taxa_list=[A, B])'
        """
        sorted_taxa = sorted(self.taxa)
        n_taxa = len(sorted_taxa)
        
        # Truncate taxon list at 10, add dots if longer
        if n_taxa <= 10:
            taxa_list_str = ", ".join(sorted_taxa)
        else:
            taxa_list_str = ", ".join(sorted_taxa[:10]) + ", ..."
        
        return (
            f"MixedPhyNetwork(nodes={self.number_of_nodes()}, "
            f"edges={self.number_of_edges()}, "
            f"taxa={n_taxa}, "
            f"taxa_list=[{taxa_list_str}])"
        )
    
    def __contains__(self, node_id: T) -> bool:
        """
        Check if node is in the network.
        
        Parameters
        ----------
        node_id : T
            Node identifier to check.
        
        Returns
        -------
        bool
            True if node is in the network, False otherwise.
        """
        return node_id in self._graph
    
    def __iter__(self) -> Iterator[T]:
        """
        Iterate over nodes.
        
        Returns
        -------
        Iterator[T]
            Iterator over node identifiers.
        """
        return iter(self._graph.nodes)
    
    def __len__(self) -> int:
        """
        Return the number of nodes.
        
        Returns
        -------
        int
            Number of nodes.
        """
        return self.number_of_nodes()
    
    # ========== Graph Operations ==========
    
    def copy(self) -> 'MixedPhyNetwork':
        """
        Create a copy of the network.
        
        Returns a shallow copy of the network. Cached properties are not
        copied but will be recomputed on first access.
        
        Returns
        -------
        MixedPhyNetwork
            A copy of the network.
        
        Examples
        --------
        >>> net = MixedPhyNetwork(undirected_edges=[(3, 1)], taxa={1: "A"})
        >>> net2 = net.copy()
        >>> net.number_of_nodes()
        2
        >>> net2.number_of_nodes()
        2
        """
        # Create new network by copying internal structures
        new_net = MixedPhyNetwork.__new__(MixedPhyNetwork)
        new_net._graph = self._graph.copy()
        new_net._node_to_label = self._node_to_label.copy()
        new_net._label_to_node = self._label_to_node.copy()
        return new_net
