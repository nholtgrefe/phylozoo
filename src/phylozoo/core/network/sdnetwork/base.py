"""
Mixed network module.

This module provides classes and functions for working with mixed phylogenetic networks.
"""

import math
import warnings
from functools import cached_property
from typing import Any, Dict, Iterator, List, Optional, Set, Tuple, TypeVar, Union

from ...primitives.m_multigraph import MixedMultiGraph
from ...primitives.m_multigraph.operations import is_connected, has_self_loops

T = TypeVar('T')


class MixedPhyNetwork:
    """
    A mixed phylogenetic network.
    
    A MixedPhyNetwork is a weakly connected, mixed multigraph (with both directed and 
    undirected edges) representing a phylogenetic network structure. This is an abstract 
    network type which may have undirected cycles, used for canonical forms and to address 
    unidentifiability issues. For semi-directed phylogenetic networks without undirected 
    cycles, use the SemiDirectedPhyNetwork subclass.
    
    It consists of:
    - **Leaf nodes**: Nodes with no outgoing directed edges, each with a taxon label
    - **Tree nodes**: Internal nodes with in-degree 0 and total degree >= 3
    - **Hybrid nodes**: Internal nodes with in-degree >= 2 and total degree = in-degree + 1
    
    A MixedPhyNetwork is obtained from a directed phylogenetic LSA (Least Stable
    Ancestor) network by undirecting all non-hybrid edges, optionally undirecting all 
    hybrid edges for selected hybrid nodes (if one hybrid edge is undirected, all partner 
    hybrid edges are undirected), and suppressing degree-2 nodes.
    
    Notes
    -----
    The class uses composition with ``MixedMultiGraph`` and is immutable after initialization; 
    construct via ``nodes``/``directed_edges``/``undirected_edges``, from a prebuilt 
    ``MixedMultiGraph``, or load from a file/eNewick string.

    Parameters
    ----------
    directed_edges : Optional[List[Union[Tuple[T, T], Tuple[T, T, int], Dict[str, Any]]]], optional
        List of directed edges. Formats:
        - (u, v) tuples (key auto-generated)
        - (u, v, key) tuples (explicit key)
        - Dict with 'u', 'v' and optional 'key' (for parallel edges) plus edge attributes
        
        Edge attributes (validated):
        - branch_length (float)
        - bootstrap (float in [0.0, 1.0])
        - gamma (float in [0.0, 1.0], hybrid edges only; for each hybrid node, all 
        incoming gammas must sum to 1.0) 
        Use a different attribute name (e.g., 'gamma2') for non-validated and/or additional
        attributes.
        
        Can be empty or None for empty/single-node networks. By default None.
    undirected_edges : Optional[List[Union[Tuple[T, T], Tuple[T, T, int], Dict[str, Any]]]], optional
        List of undirected edges. Formats:
        - (u, v) tuples (key auto-generated)
        - (u, v, key) tuples (explicit key)
        - Dict with 'u', 'v' and optional 'key' (for parallel edges) plus edge attributes
        
        Edge attributes (validated):
        - branch_length (float)
        - bootstrap (float in [0.0, 1.0])
        
        Note: Undirected edges cannot have gamma values.
        
        Can be empty or None for empty/single-node networks. By default None.
    nodes : Optional[List[Union[T, Tuple[T, Dict[str, Any]]]]], optional
        List of nodes. Formats:
        - Simple node IDs: `1`, `"node1"`, etc.
        - Tuples: `(node_id, {'label': '...','attr': ...})`
        
        Node attributes (validated):
        - label: string, unique across all nodes; use another key for non-string data.
        
        Leaves without labels are auto-labeled. Leaf-labels are referred to as `taxa`.
        Use a different attribute name (e.g., 'label2') for non-validated and/or additional
        attributes.

        Can be empty or None. By default None.
    
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
    
    Examples
    --------
    >>> # Initialize with nodes and labels
    >>> net = MixedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 4)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
    ... )
    >>> net.taxa
    {'A', 'B', 'C'}
    >>> # Partial labels - uncovered leaves get auto-generated labels
    >>> net2 = MixedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 4), (3, 5)],
    ...     nodes=[(1, {'label': 'A'})]
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
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (4, {'label': 'C'})]
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
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (3, {'label': 'C'})]
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
        nodes: Optional[List[Union[T, Tuple[T, Dict[str, Any]]]]] = None,
    ) -> None:
        """
        Initialize a mixed phylogenetic network.
        
        Parameters
        ----------
        directed_edges : Optional[List[Union[Tuple[T, T], Tuple[T, T, int], Dict[str, Any]]]], optional
            List of directed edges. Formats:
            - (u, v) tuples (key auto-generated)
            - (u, v, key) tuples (explicit key)
            - Dict with 'u', 'v' and optional 'key' plus edge attributes
            
            Edge attributes (validated):
            - branch_length (float)
            - bootstrap (float in [0.0, 1.0])
            - gamma (float in [0.0, 1.0], hybrid edges only; all incoming gammas must sum to 1.0)
              Use a different attribute name (e.g., 'gamma2') for non-validated values.
            
            Can be empty list or None for empty/single-node networks. By default None.
        undirected_edges : Optional[List[Union[Tuple[T, T], Tuple[T, T, int], Dict[str, Any]]]], optional
            List of undirected edges. Formats:
            - (u, v) tuples (key auto-generated)
            - (u, v, key) tuples (explicit key)
            - Dict with 'u', 'v' and optional 'key' plus edge attributes
            
            Edge attributes (validated):
            - branch_length (float)
            - bootstrap (float in [0.0, 1.0])
            
            Can be empty list or None for empty/single-node networks. By default None.
        nodes : Optional[List[Union[T, Tuple[T, Dict[str, Any]]]]], optional
            List of nodes. Formats:
            - Simple node IDs: `1`, `"node1"`, etc.
            - Tuples: `(node_id, {'label': '...','attr': ...})` (NetworkX-style)
            
            Node attributes (validated):
            - label: string, unique across all nodes; use another key for non-string data.
            
            Leaves without labels will get auto-generated labels. Can be empty list or None.
            By default None.
        
        Examples
        --------
        >>> # Simple network with labels
        >>> net = MixedPhyNetwork(
        ...     undirected_edges=[(3, 1), (3, 2)],
        ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        ... )
        >>> # Mix of labeled and unlabeled nodes
        >>> net = MixedPhyNetwork(
        ...     undirected_edges=[(3, 1), (3, 2)],
        ...     nodes=[3, (1, {'label': 'A'})]  # 3 has no label, 1 has label 'A'
        ... )
        >>> # Single-node network
        >>> net = MixedPhyNetwork(nodes=[(1, {'label': 'A'})])
        >>> # Empty network
        >>> net = MixedPhyNetwork(directed_edges=[], undirected_edges=[])
        >>> # With edge attributes
        >>> net = MixedPhyNetwork(
        ...     undirected_edges=[
        ...         {'u': 3, 'v': 1, 'branch_length': 0.5, 'bootstrap': 0.95},
        ...         {'u': 3, 'v': 2, 'branch_length': 0.3}
        ...     ],
        ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        ... )
        """
        # Default to empty lists
        if directed_edges is None:
            directed_edges = []
        if undirected_edges is None:
            undirected_edges = []
        
        # Warn if both edge lists are empty (but allow it)
        if not directed_edges and not undirected_edges:
            warnings.warn(
                "Initializing MixedPhyNetwork with empty edges lists. "
                "This creates an empty network or a network with isolated nodes.",
                UserWarning,
                stacklevel=2
            )
        
        self._graph: MixedMultiGraph[T] = MixedMultiGraph(
            directed_edges=directed_edges,
            undirected_edges=undirected_edges
        )
        self._node_to_label: Dict[T, str] = {}
        self._label_to_node: Dict[str, T] = {}
        
        # Step 1: Add all nodes to the graph (including attributes like label)
        # Labels are validated and added to dictionaries during this step
        self._add_nodes_to_graph(nodes)
        
        # Step 2: Auto-label any uncovered leaves
        # All leaves are guaranteed to have labels after this step
        self._auto_label_unlabeled_leaves()
        
        # Step 3: Validate the network structure
        if not self.validate():
            raise ValueError("Network validation failed")
    
    def _add_label_to_dicts(self, node_id: T, label: str) -> None:
        """
        Add a node-label mapping to both label dictionaries with validation.
        
        This helper method validates the label (string type, uniqueness) and adds
        it to both _node_to_label and _label_to_node dictionaries. This is the
        primary validation point for label constraints.
        
        Parameters
        ----------
        node_id : T
            Node identifier.
        label : str
            Label to add. Must be a string and unique.
        
        Raises
        ------
        ValueError
            If label is not a string, or if the label is already used by a different node.
        """
        # Validate string type
        if not isinstance(label, str):
            raise ValueError(
                f"Node {node_id} has non-string label '{label}' (type: {type(label).__name__}). "
                f"Labels must be strings. For non-string metadata, store it under a "
                f"different node attribute instead of 'label'."
            )
        
        # Check for duplicate labels
        if label in self._label_to_node and self._label_to_node[label] != node_id:
            existing_node = self._label_to_node[label]
            raise ValueError(
                f"Label '{label}' is already used by node {existing_node}. "
                f"Each label must be unique."
            )
        
        # Add mapping
        self._node_to_label[node_id] = label
        self._label_to_node[label] = node_id
    
    def _add_nodes_to_graph(
        self,
        nodes: Optional[List[Union[T, Tuple[T, Dict[str, Any]]]]],
    ) -> None:
        """
        Add all nodes to the underlying graph with their attributes.
        
        If a node has a 'label' attribute, it is immediately validated and added to
        the label dictionaries using _add_label_to_dicts. Label validation (string
        type and uniqueness) is performed during this step.
        
        Parameters
        ----------
        nodes : Optional[List[Union[T, Tuple[T, Dict[str, Any]]]]]
            Node specifications as simple IDs or (node_id, attr_dict) tuples.
        
        Raises
        ------
        ValueError
            If a node tuple does not provide a dict of attributes, or if duplicate labels are found.
        """
        if not nodes:
            return
        
        for node_spec in nodes:
            if isinstance(node_spec, tuple) and len(node_spec) == 2:
                node_id, attrs = node_spec
                if not isinstance(attrs, dict):
                    raise ValueError(
                        f"Node tuple must be (node_id, dict_of_attributes), got {node_spec}"
                    )
                self._graph.add_node(node_id, **attrs)
                # If label is present, validate and add it to dictionaries immediately
                if 'label' in attrs:
                    self._add_label_to_dicts(node_id, attrs['label'])
            else:
                node_id = node_spec  # type: ignore[assignment]
                self._graph.add_node(node_id)
    
    def _auto_label_unlabeled_leaves(self) -> None:
        """
        Auto-label any leaves that do not yet have labels.
        
        Generates labels based on node IDs, ensuring no duplicate labels. This method
        guarantees that all leaves have labels after execution. Labels are validated
        (string type and uniqueness) via _add_label_to_dicts.
        
        The graph node attribute is also updated to store the label.
        """
        all_leaves: Set[T] = self.leaves
        labeled_leaves = {leaf for leaf in all_leaves if leaf in self._node_to_label}
        uncovered_leaves = all_leaves - labeled_leaves
        
        for leaf_id in uncovered_leaves:
            label = str(leaf_id)
            if label in self._label_to_node:
                counter = 1
                unique_label = f"{label}_{counter}"
                while unique_label in self._label_to_node:
                    counter += 1
                    unique_label = f"{label}_{counter}"
                label = unique_label
            # Add to dictionaries using helper (validates string type and uniqueness),
            # then update graph attribute
            self._add_label_to_dicts(leaf_id, label)
            # Update graph node attribute (need to update both directed and undirected graphs)
            if leaf_id in self._graph._directed:
                self._graph._directed.nodes[leaf_id]['label'] = label
            if leaf_id in self._graph._undirected:
                self._graph._undirected.nodes[leaf_id]['label'] = label
            if leaf_id in self._graph._combined:
                self._graph._combined.nodes[leaf_id]['label'] = label
    
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
        2. No self-loops
        3. All internal nodes have degree >= 3
        4. Each node has indegree either 0 or total_degree-1
        5. Mixed network constraint (issues warning that additional checks may be added)
        6. Bootstrap values are in [0.0, 1.0]
        7. Gamma constraints (gamma only on hybrid edges, sum to 1.0 if specified)
        
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
            Also issued for empty and single-node networks.
        
        Notes
        -----
        This method performs validation checks but issues a warning that additional
        checks may be added later. For fully validated networks with additional
        constraints, use SemiDirectedPhyNetwork.
        Empty networks (no nodes) are considered valid but will raise a warning.
        Single-node networks are also valid but will raise a warning.
        """
        # Empty networks are valid
        if self.number_of_nodes() == 0:
            warnings.warn(
                "Empty network (no nodes) detected. While valid, this may not be useful for phylogenetic analysis.",
                UserWarning,
                stacklevel=2
            )
            return True

        # Single-node networks are valid only if they have no self-loops
        if self.number_of_nodes() == 1:
            if has_self_loops(self._graph):
                raise ValueError("Self-loops are not allowed in MixedPhyNetwork.")
            warnings.warn(
                "Single-node network detected. While valid, this may not be useful for phylogenetic analysis.",
                UserWarning,
                stacklevel=2
            )
            return True
        
        # 1. Check that network is connected (weakly connected)
        if not is_connected(self._graph):
            raise ValueError(
                "Network is not connected. All nodes must be in a single connected component."
            )
        
        # 2. Disallow self-loops
        if has_self_loops(self._graph):
            raise ValueError("Self-loops are not allowed in MixedPhyNetwork.")
        
        # 3. Validate degree constraints
        self._validate_degree_constraints()
        
        # 4. Validate mixed network constraint (issues warning)
        self._validate_mixednetwork_constraint()
        
        # 5. Validate bootstrap constraints
        self._validate_bootstrap_constraints()
        
        # 6. Validate gamma constraints
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
        >>> net = MixedPhyNetwork(undirected_edges=[(3, 1)], nodes=[(1, {'label': 'A'})])
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
        >>> net = MixedPhyNetwork(undirected_edges=[(3, 1)], nodes=[(1, {'label': 'A'})])
        >>> net.get_node_id("A")
        1
        """
        return self._label_to_node.get(label)
    
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
        ...     nodes=[(3, {'label': 'A'})]
        ... )
        >>> sorted(net.neighbors(2))
        [1, 3]
        """
        return self._graph.neighbors(v)
    
    # ========== Phylogenetic-Specific Methods ==========
    
    @cached_property
    def leaves(self) -> Set[T]:
        """
        Get the set of leaf node IDs (nodes with degree 1, or degree 0 for single-node networks).
        
        Returns
        -------
        Set[T]
            Set of leaf node identifiers. Returns a new set (which is mutable).
        """
        # Single-node network: the node is a leaf even though it has degree 0
        if self.number_of_nodes() == 1:
            return set(self._graph.nodes)
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
        >>> net = MixedPhyNetwork(undirected_edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        >>> net.taxa
        {'A', 'B'}
        """
        return {self._node_to_label[leaf] for leaf in self.leaves}
    
    @cached_property
    def hybrid_nodes(self) -> Set[T]:
        """
        Get the set of all hybrid nodes.
        
        A hybrid node is a node with in-degree >= 2 and total degree = in-degree + 1.
        
        Returns
        -------
        Set[T]
            Set of hybrid node identifiers. Returns a new set (which is mutable).
        
        Examples
        --------
        >>> net = MixedPhyNetwork(
        ...     directed_edges=[(5, 4), (6, 4)],
        ...     undirected_edges=[(4, 1)],
        ...     nodes=[(1, {'label': 'A'})]
        ... )
        >>> net.hybrid_nodes
        {4}
        """
        return {
            v for v in self._graph.nodes
            if self._graph.indegree(v) >= 2
            and self._graph.degree(v) == self._graph.indegree(v) + 1
        }
    
    @cached_property
    def hybrid_edges(self) -> Set[Tuple[T, T]]:
        """
        Get the set of all hybrid edges.
        
        Hybrid edges are all directed edges.
        
        Returns
        -------
        Set[Tuple[T, T]]
            Set of (source, target) tuples for hybrid edges. Returns a new set (which is mutable).
        
        Examples
        --------
        >>> net = MixedPhyNetwork(
        ...     directed_edges=[(3, 2), (4, 2)],
        ...     undirected_edges=[(2, 1)],
        ...     nodes=[(1, {'label': 'A'})]
        ... )
        >>> net.hybrid_edges
        {(3, 2), (4, 2)}
        """
        return set(self._graph._directed.edges())
    
    @cached_property
    def tree_nodes(self) -> Set[T]:
        """
        Get the set of all tree nodes.
        
        A tree node is a node with in-degree 0 and total degree >= 3.
        
        Returns
        -------
        Set[T]
            Set of tree node identifiers. Returns a new set (which is mutable).
        
        Examples
        --------
        >>> net = MixedPhyNetwork(
        ...     undirected_edges=[(1, 2), (1, 3), (1, 4)],
        ...     nodes=[(2, {'label': 'A'}), (3, {'label': 'B'}), (4, {'label': 'C'})]
        ... )
        >>> net.tree_nodes
        {1}
        """
        return {
            v for v in self._graph.nodes
            if self._graph.indegree(v) == 0
            and self._graph.degree(v) >= 3
        }
    
    @cached_property
    def internal_nodes(self) -> Set[T]:
        """
        Get the set of all internal nodes.
        
        Internal nodes are all nodes that are not leaves.
        
        Returns
        -------
        Set[T]
            Set of internal node identifiers. Returns a new set (which is mutable).
        
        Examples
        --------
        >>> net = MixedPhyNetwork(
        ...     directed_edges=[(5, 4), (6, 4)],
        ...     undirected_edges=[(4, 1)],
        ...     nodes=[(1, {'label': 'A'})]
        ... )
        >>> sorted(net.internal_nodes)
        [4, 5, 6]
        """
        leaves = self.leaves
        return {
            v for v in self._graph.nodes
            if v not in leaves
        }
    
    @cached_property
    def tree_edges(self) -> Set[Tuple[T, T]]:
        """
        Get the set of all tree edges.
        
        Tree edges are simply all undirected edges.
        
        Returns
        -------
        Set[Tuple[T, T]]
            Set of (source, target) tuples for tree edges. Returns a new set (which is mutable).
        
        Examples
        --------
        >>> net = MixedPhyNetwork(
        ...     directed_edges=[(5, 4), (6, 4)],
        ...     undirected_edges=[(4, 1)],
        ...     nodes=[(1, {'label': 'A'})]
        ... )
        >>> sorted(net.tree_edges)
        [(4, 1)]
        """
        return set(self._graph._undirected.edges())

    def __repr__(self) -> str:
        """
        Return string representation of the network.

        Returns
        -------
        str
            String representation showing nodes, edges, level, taxa count, and taxon list.
        
        Examples
        --------
        >>> net = MixedPhyNetwork(undirected_edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
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
        >>> net = MixedPhyNetwork(undirected_edges=[(3, 1)], nodes=[(1, {'label': 'A'})])
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
