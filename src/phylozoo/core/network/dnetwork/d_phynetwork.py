"""
Directed network module.

This module provides classes and functions for working with directed phylogenetic networks.
"""

import warnings
from functools import cached_property
from typing import Any, Dict, Iterator, List, Optional, Set, Tuple, TypeVar, Union

from ...primitives.d_multigraph.dm_graph import DirectedMultiGraph

T = TypeVar('T')


class DirectedPhyNetwork:
    """
    A directed phylogenetic network.
    
    This class uses composition with DirectedMultiGraph for graph structure and adds
    phylogenetic-specific features like leaves, node labels, and network topology methods.
    
    Leaves refer to node IDs, and taxa refer to the labels of leaves. Leaves must always
    have labels (taxa). Internal nodes may or may not be labeled. Node IDs are separate
    from labels, allowing flexible node identification.
    
    Each label must be unique across all nodes. Attempting to use a duplicate label
    will raise a ValueError.
    
    Parameters
    ----------
    edges : List[Union[Tuple[T, T], Tuple[T, T, int], Dict[str, Any]]]
        List of directed edges. Can be:
        - (u, v) tuples (key auto-generated)
        - (u, v, key) tuples (explicit key)
        - Dict with 'u', 'v' keys and optional 'key' and edge attributes
        
        Edge attributes can be set via dict format. Suggested attributes:
        - branch_length (float): Branch length for the edge
        - bootstrap (float): Bootstrap support value (typically 0.0 to 1.0)
        - gamma (float): For hybrid edges only - inheritance probability.
          The sum of gamma values for all edges entering a hybrid node must equal 1.0.
          This is validated during initialization.
        
        Other edge attributes are also supported (any key-value pairs can be added),
        but only the above three are suggested for phylogenetic networks.
        
        Must be provided (can be empty list for empty network, but a warning will be raised).
    taxa : Optional[Dict[T, str] | List[Tuple[T, str]]], optional
        Taxon labels for leaves. Can be:
        - Dict mapping leaf node IDs to taxon labels: {leaf_id: "taxon"}
        - List of (leaf_id, taxon_label) tuples
        The IDs in this mapping must be leaves (no outgoing edges).
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
    _graph : DirectedMultiGraph[T]
        Internal graph structure using DirectedMultiGraph.
        **Warning:** Do not modify directly. Use class methods instead.
    _node_to_label : Dict[T, str]
        Mapping from node IDs to labels. Only nodes with explicit labels are included.
        Leaves always have labels (taxa), but internal nodes may be unlabeled.
    _label_to_node : Dict[str, T]
        Reverse mapping from labels to node IDs (for quick lookup).
    
    Notes
    -----
    This class is immutable after initialization. To create a network,
    build it using DirectedMultiGraph and then create a DirectedPhyNetwork from it,
    or initialize it with edges and taxa.
    
    Examples
    --------
    >>> # Initialize with taxa mapping
    >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"})
    >>> net.taxa
    {'A', 'B'}
    >>> net.root_node
    3
    >>> net.is_tree()
    True
    >>> # Partial taxa mapping - uncovered leaves get auto-generated labels
    >>> net2 = DirectedPhyNetwork(edges=[(3, 1), (3, 2), (3, 4)], taxa={1: "A"})
    >>> net2.taxa  # 2 and 4 are auto-labeled
    {'A', '2', '4'}
    >>> # Network with branch lengths and bootstrap support
    >>> net3 = DirectedPhyNetwork(
    ...     edges=[
    ...         {'u': 3, 'v': 1, 'branch_length': 0.5, 'bootstrap': 0.95},
    ...         {'u': 3, 'v': 2, 'branch_length': 0.3, 'bootstrap': 0.87}
    ...     ],
    ...     taxa={1: "A", 2: "B"}
    ... )
    >>> net3.get_branch_length(3, 1)
    0.5
    >>> net3.get_bootstrap(3, 1)
    0.95
    >>> # Network with hybrid node and gamma values
    >>> net4 = DirectedPhyNetwork(
    ...     edges=[
    ...         {'u': 5, 'v': 4, 'gamma': 0.6},
    ...         {'u': 6, 'v': 4, 'gamma': 0.4},  # Sum = 1.0
    ...         {'u': 4, 'v': 1, 'branch_length': 0.2}
    ...     ],
    ...     taxa={1: "A"}
    ... )
    >>> net4.get_gamma(5, 4)
    0.6
    >>> net4.get_gamma(6, 4)
    0.4
    """
    
    def __init__(
        self,
        edges: List[Union[Tuple[T, T], Tuple[T, T, int], Dict[str, Any]]],
        taxa: Optional[Dict[T, str] | List[Tuple[T, str]]] = None,
        internal_node_labels: Optional[Dict[T, str] | List[Tuple[T, str]]] = None,
    ) -> None:
        """
        Initialize a directed phylogenetic network.
        
        Parameters
        ----------
        edges : List[Union[Tuple[T, T], Tuple[T, T, int], Dict[str, Any]]]
            List of directed edges. Must be provided (can be empty list for empty network,
            but a warning will be raised).
        taxa : Optional[Dict[T, str] | List[Tuple[T, str]]], optional
            Taxon labels for leaves. Can be:
            - Dict mapping leaf IDs to taxon labels: {leaf_id: "taxon"}
            - List of (leaf_id, taxon_label) tuples
            The IDs must be leaves (no outgoing edges). Not all leaves need to be
            covered - uncovered leaves get auto-generated labels. By default None.
        internal_node_labels : Optional[Dict[T, str] | List[Tuple[T, str]]], optional
            Labels for internal nodes (optional). Can be:
            - Dict mapping node IDs to labels: {node_id: "label"}
            - List of (node_id, label) tuples
            Only internal nodes (non-leaves) can be labeled via this parameter.
            Leaves must be labeled via 'taxa'. By default None.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"})  # Dict: leaf_id -> taxon
        >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa=[(1, "A"), (2, "B")])  # List of tuples
        >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"}, internal_node_labels={3: "root"})
        >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"}, internal_node_labels=[(3, "root")])
        >>> # With edge attributes
        >>> net = DirectedPhyNetwork(
        ...     edges=[
        ...         {'u': 3, 'v': 1, 'branch_length': 0.5, 'bootstrap': 0.95},
        ...         {'u': 3, 'v': 2, 'branch_length': 0.3}
        ...     ],
        ...     taxa={1: "A", 2: "B"}
        ... )
        """
        # Warn if edges is empty
        if not edges:
            warnings.warn(
                "Initializing DirectedPhyNetwork with empty edges list. "
                "This creates an empty network.",
                UserWarning,
                stacklevel=2
            )
        
        self._graph: DirectedMultiGraph[T] = DirectedMultiGraph(edges=edges)
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
        # Identify all leaves (nodes with no outgoing edges)
        all_leaves: Set[T] = {
            node for node in self._graph.nodes
            if self._graph.outdegree(node) == 0
        }
        
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
                            f"(has outgoing edges)"
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
                            f"(has outgoing edges)"
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
    
    def validate(self) -> bool:
        """
        Validate the network structure and edge attributes.
        
        Checks:
        - Network structure is valid
        - For each hybrid node, if any incoming edge has a gamma value,
          the sum of all gamma values on ALL parallel incoming edges must equal 1.0
        
        Returns
        -------
        bool
            True if the network is valid, False otherwise.
        
        Raises
        ------
        ValueError
            If gamma constraints are violated.
        
        Notes
        -----
        This method checks that the network structure is valid and that
        gamma values for hybrid edges sum to 1.0 per hybrid node.
        """
        # Check gamma constraints for hybrid nodes
        for hybrid_node in self.hybrid_nodes:
            gamma_sum = 0.0
            # Use incident_parent_edges to get all incoming edges (including parallel edges)
            for edge in self.incident_parent_edges(hybrid_node, keys=True, data=True):
                if len(edge) == 4:  # (u, v, key, data)
                    _, _, _, edge_data = edge
                    gamma = edge_data.get('gamma')
                    if gamma is not None:
                        gamma_sum += gamma
            
            # If any gamma values are set, they must sum to 1.0
            if gamma_sum > 0 and abs(gamma_sum - 1.0) > 1e-10:
                raise ValueError(
                    f"Hybrid node {hybrid_node} has gamma values that sum to {gamma_sum}, "
                    f"but must sum to 1.0"
                )
        
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
        >>> net = DirectedPhyNetwork(edges=[(3, 1)], taxa={1: "A"})
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
        >>> net = DirectedPhyNetwork(edges=[(3, 1)], taxa={1: "A"})
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
    
    def get_edge_attribute(self, u: T, v: T, key: Optional[int] = None, attr: str = 'branch_length') -> Optional[Any]:
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
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(
        ...     edges=[
        ...         {'u': 3, 'v': 1, 'key': 0, 'branch_length': 0.5},
        ...         {'u': 3, 'v': 1, 'key': 1, 'branch_length': 0.7}  # Parallel edge
        ...     ],
        ...     taxa={1: "A"}
        ... )
        >>> net.get_edge_attribute(3, 1, key=0, attr='branch_length')
        0.5
        >>> net.get_edge_attribute(3, 1, key=1, attr='branch_length')
        0.7
        >>> net.get_edge_attribute(3, 1, attr='branch_length')  # Raises ValueError
        """
        if not self.has_edge(u, v):
            return None
        
        # Check if parallel edges exist
        edges_data = self._graph._graph[u][v]
        num_edges = len(edges_data)
        
        if num_edges == 0:
            return None
        elif num_edges == 1:
            # Single edge - key not needed
            first_key = next(iter(edges_data))
            return edges_data[first_key].get(attr)
        else:
            # Multiple parallel edges - key is required
            if key is None:
                raise ValueError(
                    f"Multiple parallel edges exist between {u} and {v}. "
                    f"Must specify 'key' parameter to get attribute from a specific edge."
                )
            if key not in edges_data:
                return None
            return edges_data[key].get(attr)
    
    def get_edge_attributes(self, u: T, v: T, attr: str = 'branch_length') -> Dict[int, Any]:
        """
        Get attribute values for all parallel edges between u and v.
        
        Parameters
        ----------
        u, v : T
            Edge endpoints.
        attr : str, default 'branch_length'
            Attribute name. Suggested attributes: 'branch_length', 'bootstrap', 'gamma'.
            Any edge attribute can be accessed.
        
        Returns
        -------
        Dict[int, Any]
            Dictionary mapping edge keys to attribute values.
            Keys without the attribute are omitted.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(
        ...     edges=[
        ...         {'u': 3, 'v': 1, 'key': 0, 'branch_length': 0.5},
        ...         {'u': 3, 'v': 1, 'key': 1, 'branch_length': 0.7}
        ...     ],
        ...     taxa={1: "A"}
        ... )
        >>> net.get_edge_attributes(3, 1, 'branch_length')
        {0: 0.5, 1: 0.7}
        """
        if not self.has_edge(u, v):
            return {}
        
        result: Dict[int, Any] = {}
        for edge_key, edge_data in self._graph._graph[u][v].items():
            if attr in edge_data:
                result[edge_key] = edge_data[attr]
        return result
    
    def get_branch_length(self, u: T, v: T, key: Optional[int] = None) -> Optional[float]:
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
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(
        ...     edges=[{'u': 3, 'v': 1, 'branch_length': 0.5}],
        ...     taxa={1: "A"}
        ... )
        >>> net.get_branch_length(3, 1)
        0.5
        """
        return self.get_edge_attribute(u, v, key, 'branch_length')
    
    def get_bootstrap(self, u: T, v: T, key: Optional[int] = None) -> Optional[float]:
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
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(
        ...     edges=[{'u': 3, 'v': 1, 'bootstrap': 0.95}],
        ...     taxa={1: "A"}
        ... )
        >>> net.get_bootstrap(3, 1)
        0.95
        """
        return self.get_edge_attribute(u, v, key, 'bootstrap')
    
    def get_gamma(self, u: T, v: T, key: Optional[int] = None) -> Optional[float]:
        """
        Get gamma value for a hybrid edge.
        
        Gamma values represent inheritance probabilities for hybrid nodes.
        The sum of gamma values for all edges entering a hybrid node must equal 1.0.
        
        Parameters
        ----------
        u, v : T
            Edge endpoints (v must be a hybrid node).
        key : Optional[int], optional
            Edge key for parallel edges. Required if multiple parallel edges exist.
        
        Returns
        -------
        Optional[float]
            Gamma value, or None if not set.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(
        ...     edges=[
        ...         {'u': 5, 'v': 4, 'gamma': 0.6},
        ...         {'u': 6, 'v': 4, 'gamma': 0.4},  # Sum = 1.0
        ...         {'u': 4, 'v': 1}
        ...     ],
        ...     taxa={1: "A"}
        ... )
        >>> net.get_gamma(5, 4)
        0.6
        >>> net.get_gamma(6, 4)
        0.4
        """
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
            Number of edges.
        """
        return self._graph.number_of_edges()
    
    def has_edge(self, u: T, v: T, key: Optional[int] = None) -> bool:
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
        
        Returns
        -------
        bool
            True if edge exists, False otherwise.
        """
        return self._graph.has_edge(u, v, key)
    
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
            Total degree (in-degree + out-degree).
        """
        return self._graph.degree(v)
    
    def indegree(self, v: T) -> int:
        """
        Return the in-degree of node v.
        
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
        Return the out-degree of node v.
        
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
    
    def parents(self, v: T) -> Iterator[T]:
        """
        Return an iterator over parent nodes of v.
        
        Parameters
        ----------
        v : T
            Node identifier.
        
        Returns
        -------
        Iterator[T]
            Iterator over parent nodes.
        """
        return self._graph.predecessors(v)
    
    def children(self, v: T) -> Iterator[T]:
        """
        Return an iterator over child nodes of v.
        
        Parameters
        ----------
        v : T
            Node identifier.
        
        Returns
        -------
        Iterator[T]
            Iterator over child nodes.
        """
        return self._graph.successors(v)
    
    def neighbors(self, v: T) -> Iterator[T]:
        """
        Return an iterator over neighbors of node v (parents and children).
        
        Parameters
        ----------
        v : T
            Node identifier.
        
        Returns
        -------
        Iterator[T]
            Iterator over neighbors.
        """
        return self._graph.neighbors(v)
    
    def incident_parent_edges(self, v: T, keys: bool = False, data: bool = False) -> Iterator[Tuple[T, T] | Tuple[T, T, int] | Tuple[T, T, int, Dict[str, Any]]]:
        """
        Return an iterator over edges entering node v (from parent nodes).
        
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
            Iterator over incoming edges as (u, v) or (u, v, key) or (u, v, key, data).
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(
        ...     edges=[
        ...         {'u': 1, 'v': 2, 'branch_length': 0.5},
        ...         {'u': 3, 'v': 2, 'branch_length': 0.3}
        ...     ],
        ...     taxa={2: "A"}
        ... )
        >>> list(net.incident_parent_edges(2))
        [(1, 2), (3, 2)]
        >>> list(net.incident_parent_edges(2, data=True))
        [(1, 2, {'branch_length': 0.5}), (3, 2, {'branch_length': 0.3})]
        """
        return self._graph.incident_parent_edges(v, keys=keys, data=data)
    
    def incident_child_edges(self, v: T, keys: bool = False, data: bool = False) -> Iterator[Tuple[T, T] | Tuple[T, T, int] | Tuple[T, T, int, Dict[str, Any]]]:
        """
        Return an iterator over edges leaving node v (to child nodes).
        
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
            Iterator over outgoing edges as (v, u) or (v, u, key) or (v, u, key, data).
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(
        ...     edges=[
        ...         {'u': 1, 'v': 2, 'branch_length': 0.5},
        ...         {'u': 1, 'v': 3, 'branch_length': 0.3}
        ...     ],
        ...     taxa={2: "A", 3: "B"}
        ... )
        >>> list(net.incident_child_edges(1))
        [(1, 2), (1, 3)]
        >>> list(net.incident_child_edges(1, data=True))
        [(1, 2, {'branch_length': 0.5}), (1, 3, {'branch_length': 0.3})]
        """
        return self._graph.incident_child_edges(v, keys=keys, data=data)
    
    # ========== Phylogenetic-Specific Methods ==========
    
    @cached_property
    def leaves(self) -> Set[T]:
        """
        Get the set of leaf node IDs (nodes with no outgoing edges).
        
        Returns
        -------
        Set[T]
            Set of leaf node identifiers. Returns a new set (which is mutable).
        """
        return {node for node in self._graph.nodes if self._graph.outdegree(node) == 0}
    
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
        >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"})
        >>> net.taxa
        {'A', 'B'}
        """
        return {self._node_to_label[leaf] for leaf in self.leaves}
    
    @cached_property
    def internal_nodes(self) -> List[T]:
        """
        Return a list of all internal (non-leaf) nodes.
        
        Returns
        -------
        List[T]
            List of internal node identifiers.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"})
        >>> net.internal_nodes
        [3]
        """
        return [v for v in self._graph.nodes if v not in self.leaves]
    
    @cached_property
    def root_node(self) -> T:
        """
        Return the root node of the network.
        
        The root is the node with in-degree 0.
        
        Returns
        -------
        T
            Root node identifier.
        
        Raises
        ------
        ValueError
            If there is no root node or multiple root nodes.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"})
        >>> net.root_node
        3
        """
        roots = [v for v in self._graph.nodes if self._graph.indegree(v) == 0]
        if len(roots) == 0:
            raise ValueError("Network has no root node")
        if len(roots) > 1:
            raise ValueError(f"Network has multiple root nodes: {roots}")
        return roots[0]
    
    @cached_property
    def hybrid_nodes(self) -> List[T]:
        """
        Return a list of all hybrid nodes.
        
        A hybrid node is a node with in-degree 2 and total degree 3.
        
        Returns
        -------
        List[T]
            List of hybrid node identifiers.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(edges=[(5, 4), (6, 4), (4, 1)], taxa={1: "A"})
        >>> net.hybrid_nodes
        [4]
        """
        return [
            v for v in self._graph.nodes
            if self._graph.indegree(v) == 2 and self._graph.degree(v) == 3
        ]
    
    @cached_property
    def hybrid_edges(self) -> List[Tuple[T, T]]:
        """
        Return a list of all hybrid edges.
        
        Hybrid edges are edges that point into hybrid nodes.
        
        Returns
        -------
        List[Tuple[T, T]]
            List of (source, target) tuples for hybrid edges.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(edges=[(3, 2), (4, 2)], taxa={2: "A"})
        >>> net.hybrid_edges
        [(3, 2), (4, 2)]
        """
        res = []
        for v in self.hybrid_nodes:
            for p in self._graph.predecessors(v):
                res.append((p, v))
        return res
    
    @cached_property
    def non_hybrid_edges(self) -> List[Tuple[T, T]]:
        """
        Return a list of all non-hybrid edges.
        
        Returns
        -------
        List[Tuple[T, T]]
            List of (source, target) tuples for non-hybrid edges.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"})
        >>> len(net.non_hybrid_edges)
        2
        """
        hybrid_edges = set(self.hybrid_edges)
        # edges() returns EdgeView which iterates as (u, v) tuples
        return [(u, v) for (u, v) in self._graph.edges() if (u, v) not in hybrid_edges]
    
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
    
    def is_tree(self) -> bool:
        """
        Check if the network is a tree.
        
        A tree has no hybrid nodes.
        
        Returns
        -------
        bool
            True if the network is a tree, False otherwise.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"})
        >>> net.is_tree()
        True
        """
        return len(self.hybrid_nodes) == 0
    
    # ========== Graph Operations ==========
    
    def copy(self) -> 'DirectedPhyNetwork':
        """
        Create a copy of the network.
        
        Returns a shallow copy of the network. Cached properties are not
        copied but will be recomputed on first access.
        
        Returns
        -------
        DirectedPhyNetwork
            A copy of the network.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(edges=[(3, 1)], taxa={1: "A"})
        >>> net2 = net.copy()
        >>> net.number_of_nodes()
        2
        >>> net2.number_of_nodes()
        2
        """
        # Create new network by copying internal structures
        new_net = DirectedPhyNetwork.__new__(DirectedPhyNetwork)
        new_net._graph = self._graph.copy()
        new_net._node_to_label = self._node_to_label.copy()
        new_net._label_to_node = self._label_to_node.copy()
        return new_net
    
    def __repr__(self) -> str:
        """
        Return string representation of the network.
        
        Returns
        -------
        str
            String representation showing nodes, edges, level, taxa count, and taxon list.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], taxa={1: "A", 2: "B"})
        >>> repr(net)
        'DirectedPhyNetwork(nodes=3, edges=2, level=0, taxa=2, taxa_list=[A, B])'
        """
        sorted_taxa = sorted(self.taxa)
        n_taxa = len(sorted_taxa)
        
        # Truncate taxon list at 10, add dots if longer
        if n_taxa <= 10:
            taxa_list_str = ", ".join(sorted_taxa)
        else:
            taxa_list_str = ", ".join(sorted_taxa[:10]) + ", ..."
        
        return (
            f"DirectedPhyNetwork(nodes={self.number_of_nodes()}, "
            f"edges={self.number_of_edges()}, "
            f"level={self.level}, "
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

