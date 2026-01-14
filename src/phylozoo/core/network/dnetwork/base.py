"""
Directed network module.

This module provides classes and functions for working with directed phylogenetic networks.
"""

import math
import warnings
from functools import cached_property
from typing import Any, Iterator, TypeVar

import networkx as nx

from ....utils.exceptions import (
    PhyloZooNetworkDegreeError, 
    PhyloZooNetworkStructureError, 
    PhyloZooNetworkAttributeError, 
    PhyloZooValueError,
    PhyloZooTypeError,
    PhyloZooEmptyNetworkWarning,
    PhyloZooSingleNodeNetworkWarning,
)

from ...primitives.d_multigraph import DirectedMultiGraph
from ...primitives.d_multigraph.features import is_connected, has_self_loops
from ....utils.validation import validation_aware
from ....utils.io import IOMixin

T = TypeVar('T')


@validation_aware(allowed=["validate", "_validate_*"], default=["validate"])
class DirectedPhyNetwork(IOMixin):
    """
    A directed phylogenetic network.
    
    A DirectedPhyNetwork is a weakly connected, directed acyclic multi-graph (DAG) representing 
    a phylogenetic network structure. It consists of:
    - **Root node**: Exactly one node with in-degree 0
    - **Leaf nodes**: Nodes with out-degree 0, each with in-degree 1 and a taxon label
    - **Tree nodes**: Internal nodes (non-root, non-leaf) with in-degree 1 and out-degree >= 2
    - **Hybrid nodes**: Internal nodes with in-degree >= 2 and out-degree 1
    For technical reasons, an empty network or single-node network (where root and leaf are the same node) 
    is also valid. Internal nodes may be unlabeled.
    
    Notes
    -----
    The class uses composition with ``DirectedMultiGraph`` and is immutable after initialization; 
    construct via ``nodes``/``edges``, from a prebuilt ``DirectedMultiGraph``, or load from a file/eNewick 
    string.

    Parameters
    ----------
    edges : list[tuple[T, T] | tuple[T, T, int] | dict[str, Any]] | None, optional
        List of directed edges. Formats:
        - (u, v) tuples (key auto-generated)
        - (u, v, key) tuples (explicit key)
        - Dict with 'u', 'v' and optional 'key' (for parallel edges) plus edge attributes
        
        Edge attributes (validated):
        - branch_length (float; for set of parallel edges, all must have equal branch_length)
        - bootstrap (float in [0.0, 1.0])
        - gamma (float in [0.0, 1.0], hybrid edges only; for each hybrid node, all 
        incoming gammas must sum to 1.0) 
        Use a different attribute name (e.g., 'gamma2') for non-validated and/or additional
        attritbutes.
        
        Can be empty or None for empty/single-node networks. By default None.
    nodes : list[T | tuple[T, dict[str, Any]]] | None, optional
        List of nodes. Formats:
        - Simple node IDs: `1`, `"node1"`, etc.
        - Tuples: `(node_id, {'label': '...','attr': ...})`
        
        Node attributes (validated):
        - label: string, unique across all nodes; use another key for non-string data.
        
        Leaves without labels are auto-labeled. Leaf-labels are referred to as `taxa`.
        Use a different attribute name (e.g., 'label2') for non-validated and/or additional
        attritbutes.

        Can be empty or None. By default None.
    attributes : dict[str, Any] | None, optional
        Optional dictionary of graph-level attributes to store with the network.
        These attributes are stored in the underlying graph's `.graph` attribute
        and are preserved through copy operations. Can be used to store metadata
        like provenance, source file, creation date, etc.
        By default None.
    
    Attributes
    ----------
    _graph : DirectedMultiGraph[T]
        Internal graph structure using DirectedMultiGraph.
        **Warning:** Do not modify directly. Use class methods instead.
    _node_to_label : dict[T, str]
        Mapping from node IDs to labels. Only nodes with explicit labels are included.
        Leaves always have labels (taxa), but internal nodes may be unlabeled.
    _label_to_node : dict[str, T]
        Reverse mapping from labels to node IDs (for quick lookup).
    
    Examples
    --------
    >>> # Initialize with nodes and labels
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> net.taxa
    {'A', 'B'}
    >>> net.root_node
    3
    >>> net.is_tree()
    True
    >>> # Partial labels - uncovered leaves get auto-generated labels
    >>> net2 = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2), (3, 4)],
    ...     nodes=[(1, {'label': 'A'})]
    ... )
    >>> net2.taxa  # 2 and 4 are auto-labeled
    {'A', '2', '4'}
    >>> # Network with branch lengths and bootstrap support
    >>> net3 = DirectedPhyNetwork(
    ...     edges=[
    ...         {'u': 3, 'v': 1, 'branch_length': 0.5, 'bootstrap': 0.95},
    ...         {'u': 3, 'v': 2, 'branch_length': 0.3, 'bootstrap': 0.87}
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> net3.get_branch_length(3, 1)
    0.5
    >>> net3.get_bootstrap(3, 1)
    0.95
    >>> # Network with hybrid node and gamma values
    >>> net4 = DirectedPhyNetwork(
    ...     edges=[
    ...         (7, 5), (7, 6),  # Root to tree nodes
    ...         {'u': 5, 'v': 4, 'gamma': 0.6},  # Hybrid edge
    ...         {'u': 6, 'v': 4, 'gamma': 0.4},  # Hybrid edge (Sum = 1.0)
    ...         (5, 8), (6, 9),  # Tree nodes also have other children
    ...         (4, 1)  # Hybrid to leaf
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
    ... )
    >>> net4.get_gamma(5, 4)
    0.6
    >>> net4.get_gamma(6, 4)
    0.4
    """
    
    # I/O format configuration
    _default_format = 'enewick'
    _supported_formats = ['enewick', 'dot']
    
    def __init__(
        self,
        edges: list[tuple[T, T] | tuple[T, T, int] | dict[str, Any]] | None = None,
        nodes: list[T | tuple[T, dict[str, Any]]] | None = None,
        attributes: dict[str, Any] | None = None,
    ) -> None:
        """
        Initialize a directed phylogenetic network.
        
        Parameters
        ----------
        edges : list[tuple[T, T] | tuple[T, T, int] | dict[str, Any]] | None, optional
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
        nodes : list[T | tuple[T, dict[str, Any]]] | None, optional
            List of nodes. Formats:
            - Simple node IDs: `1`, `"node1"`, etc.
            - Tuples: `(node_id, {'label': '...','attr': ...})` (NetworkX-style)
            
            Node attributes (validated):
            - label: string, unique across all nodes; use another key for non-string data.
            
            Leaves without labels will get auto-generated labels. Can be empty list or None.
            By default None.
        attributes : dict[str, Any] | None, optional
            Optional dictionary of graph-level attributes to store with the network.
            These attributes are stored in the underlying graph's `.graph` attribute
            and are preserved through copy operations. Can be used to store metadata
            like provenance, source file, creation date, etc.
            By default None.
        
        Examples
        --------
        >>> # Simple network with labels
        >>> net = DirectedPhyNetwork(
        ...     edges=[(3, 1), (3, 2)],
        ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        ... )
        >>> # Mix of labeled and unlabeled nodes
        >>> net = DirectedPhyNetwork(
        ...     edges=[(3, 1), (3, 2)],
        ...     nodes=[3, (1, {'label': 'A'})]  # 3 has no label, 1 has label 'A'
        ... )
        >>> # Single-node network
        >>> net = DirectedPhyNetwork(nodes=[(1, {'label': 'A'})])
        >>> # Empty network
        >>> net = DirectedPhyNetwork(edges=[])
        >>> # With edge attributes
        >>> net = DirectedPhyNetwork(
        ...     edges=[
        ...         {'u': 3, 'v': 1, 'branch_length': 0.5, 'bootstrap': 0.95},
        ...         {'u': 3, 'v': 2, 'branch_length': 0.3}
        ...     ],
        ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
        ... )
        """
        # Default to empty list
        if edges is None:
            edges = []
        
        self._graph: DirectedMultiGraph[T] = DirectedMultiGraph(edges=edges, attributes=attributes)
        self._node_to_label: dict[T, str] = {}
        self._label_to_node: dict[str, T] = {}
        
        # Step 1: Add all nodes to the graph (including attributes like label)
        # Labels are validated and added to dictionaries during this step
        self._add_nodes_to_graph(nodes)
        
        # Step 2: Auto-label any uncovered leaves
        # All leaves are guaranteed to have labels after this step
        self._auto_label_unlabeled_leaves()
        
        # Step 3: Validate the network structure
        self.validate()
    
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
        PhyloZooValueError
            If the label is already used by a different node.
        PhyloZooTypeError
            If the label is not a string.
        """
        # Validate string type
        if not isinstance(label, str):
            raise PhyloZooTypeError(
                f"Node {node_id} has non-string label '{label}' (type: {type(label).__name__}). "
                f"Labels must be strings. For non-string metadata, store it under a "
                f"different node attribute instead of 'label'."
            )
        
        # Check for duplicate labels
        if label in self._label_to_node and self._label_to_node[label] != node_id:
            existing_node = self._label_to_node[label]
            raise PhyloZooValueError(
                f"Label '{label}' is already used by node {existing_node}. "
                f"Each label must be unique."
            )
        
        # Add mapping
        self._node_to_label[node_id] = label
        self._label_to_node[label] = node_id
    
    def _add_nodes_to_graph(
        self,
        nodes: list[T | tuple[T, dict[str, Any]]] | None,
    ) -> None:
        """
        Add all nodes to the underlying graph with their attributes.
        
        If a node has a 'label' attribute, it is immediately validated and added to
        the label dictionaries using _add_label_to_dicts. Label validation (string
        type and uniqueness) is performed during this step.
        
        Parameters
        ----------
        nodes : list[T | tuple[T, dict[str, Any]]] | None
            Node specifications as simple IDs or (node_id, attr_dict) tuples.
        
        Raises
        ------
        PhyloZooTypeError
            If a node tuple does not provide a dict of attributes, or if duplicate labels are found.
        """
        if not nodes:
            return
        
        for node_spec in nodes:
            if isinstance(node_spec, tuple) and len(node_spec) == 2:
                node_id, attrs = node_spec
                if not isinstance(attrs, dict):
                    raise PhyloZooTypeError(
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
        all_leaves: set[T] = self.leaves
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
            self._graph._graph.nodes[leaf_id]['label'] = label
    
    def _validate_bootstrap_constraints(self) -> None:
        """
        Validate bootstrap values are within [0.0, 1.0].
        
        Raises
        ------
        PhyloZooNetworkAttributeError
            If any edge has a bootstrap value outside [0.0, 1.0].
        """
        # Quick check: if no bootstrap values are set, skip validation
        has_bootstrap = False
        for u, v, key, data in self._graph.edges(keys=True, data=True):
            if 'bootstrap' in data:
                has_bootstrap = True
                break
        if not has_bootstrap:
            return
        
        for u, v, key, data in self._graph.edges(keys=True, data=True):
            if 'bootstrap' in data:
                bootstrap = data['bootstrap']
                if not isinstance(bootstrap, (int, float)):
                    raise PhyloZooNetworkAttributeError(
                        f"Bootstrap value on edge ({u}, {v}, key={key}) must be numeric, "
                        f"got {type(bootstrap).__name__}"
                    )
                if math.isnan(bootstrap) or bootstrap < 0.0 or bootstrap > 1.0:
                    raise PhyloZooNetworkAttributeError(
                        f"Bootstrap value on edge ({u}, {v}, key={key}) is {bootstrap}, "
                        f"but must be in [0.0, 1.0]"
                    )
    
    def _validate_gamma_constraints(self) -> None:
        """
        Validate gamma constraints for hybrid nodes.
        
        Gamma values can only be set on hybrid edges (edges pointing into hybrid nodes).
        For each hybrid node, if ANY incoming edge has a gamma value, then
        ALL incoming edges (including parallel edges) must have gamma values,
        and they must sum to exactly 1.0. If no gamma values are specified,
        no validation is performed.
        
        Raises
        ------
        PhyloZooNetworkAttributeError
            If gamma constraints are violated, including gamma set on non-hybrid edges.
        """
        # Quick check: if no gamma values are set, skip validation
        has_gamma = False
        for u, v, key, data in self._graph.edges(keys=True, data=True):
            if 'gamma' in data:
                has_gamma = True
                break
        if not has_gamma:
            return
        
        # First, check that gamma is only set on hybrid edges
        hybrid_edges_set = self.hybrid_edges  # Now contains (u, v, key) tuples
        
        for u, v, key, data in self._graph.edges(keys=True, data=True):
            if 'gamma' in data:
                # Check if this edge is a hybrid edge
                if (u, v, key) not in hybrid_edges_set:
                    raise PhyloZooNetworkAttributeError(
                        f"Gamma value can only be set on hybrid edges (edges pointing into "
                        f"hybrid nodes). Edge ({u}, {v}, key={key}) is not a hybrid edge."
                    )
        
        # Then validate gamma constraints for hybrid nodes
        for hybrid_node in self.hybrid_nodes:
            gamma_values: list[float] = []
            incoming_edges: list[tuple[T, T, int]] = []
            
            # Use incident_parent_edges to get all incoming edges (including parallel edges)
            for edge in self.incident_parent_edges(hybrid_node, keys=True, data=True):
                if len(edge) == 4:  # (u, v, key, data)
                    u, v, key, edge_data = edge
                    incoming_edges.append((u, v, key))
                    gamma = edge_data.get('gamma')
                    if gamma is not None:
                        # Validate gamma is numeric
                        if not isinstance(gamma, (int, float)):
                            raise PhyloZooNetworkAttributeError(
                                f"Gamma value on edge ({u}, {v}, key={key}) entering hybrid node "
                                f"{hybrid_node} must be numeric, got {type(gamma).__name__}"
                            )
                        # Validate gamma is in [0.0, 1.0]
                        if gamma < 0.0 or gamma > 1.0:
                            raise PhyloZooNetworkAttributeError(
                                f"Gamma value on edge ({u}, {v}, key={key}) entering hybrid node "
                                f"{hybrid_node} is {gamma}, but must be in [0.0, 1.0]"
                            )
                        gamma_values.append(gamma)
            
            # If any gamma values are set, ALL edges must have gamma values
            if len(gamma_values) > 0:
                # Check if all incoming edges have gamma values
                missing_edges: list[str] = []
                for u, v, key in incoming_edges:
                    gamma = self.get_edge_attribute(u, v, key, 'gamma')
                    if gamma is None:
                        missing_edges.append(f"({u}, {v}, key={key})")
                
                if missing_edges:
                    raise PhyloZooNetworkAttributeError(
                        f"Hybrid node {hybrid_node} has some edges with gamma values "
                        f"but others without. If ANY gamma is specified, ALL incoming edges "
                        f"must have gamma values. Missing gamma on edges: {', '.join(missing_edges)}"
                    )
                
                # All gammas are present, check they sum to 1.0
                gamma_sum = sum(gamma_values)
                if abs(gamma_sum - 1.0) > 1e-10:
                    raise PhyloZooNetworkAttributeError(
                        f"Hybrid node {hybrid_node} has gamma values that sum to {gamma_sum}, "
                        f"but must sum to exactly 1.0"
                    )
    
    def _validate_degree_constraints(self) -> None:
        """
        Validate degree constraints for nodes.
        
        Checks:
        1. Single root node (exactly one node with in-degree 0)
        2. Leaf nodes have in-degree 1 and out-degree 0
        3. Internal nodes have either (in-degree 1 and out-degree >= 2) or
           (in-degree >= 2 and out-degree 1)
        
        Raises
        ------
        PhyloZooNetworkDegreeError
            If any degree constraints are violated.
        """
        # 1. Check for single root node (using cached property)
        # Accessing root_node will raise if no root or multiple roots
        root = self.root_node
        
        # 2. Check leaf nodes: in-degree 1, out-degree 0 (using cached property)
        leaves = self.leaves
        for leaf in leaves:
            indeg = self._graph.indegree(leaf)
            if indeg != 1:
                raise PhyloZooNetworkDegreeError(
                    f"Leaf node {leaf} has in-degree {indeg}, but must have in-degree 1"
                )
        
        # 3. Check internal nodes (non-root, non-leaf)
        # Use cached properties: tree_nodes and hybrid_nodes cover all valid internal nodes
        # Any node that's not root, not leaf, and not in tree_nodes or hybrid_nodes is invalid
        tree_nodes = set(self.tree_nodes)
        hybrid_nodes = set(self.hybrid_nodes)
        
        # Use internal_nodes cached property (already excludes root and leaves)
        all_internal = set(self.internal_nodes)
        valid_internal = tree_nodes | hybrid_nodes
        
        invalid_nodes = all_internal - valid_internal
        if invalid_nodes:
            for node in invalid_nodes:
                indeg = self._graph.indegree(node)
                outdeg = self._graph.outdegree(node)
                raise PhyloZooNetworkDegreeError(
                    f"Internal node {node} has in-degree {indeg} and out-degree {outdeg}. "
                    f"Internal nodes must have either (in-degree 1 and out-degree >= 2) "
                    f"or (in-degree >= 2 and out-degree 1)."
                )
    
    def _validate_branchlength_constraints(self) -> None:
        """
        Validate branch length constraints for parallel edges.
        
        For each set of parallel edges between nodes u and v, this method ensures:
        1. If one edge has a branch_length attribute, all parallel edges must have branch_length
        2. All branch_length values must be the same across parallel edges
        
        Raises
        ------
        ValueError
            If branch length constraints are violated for any set of parallel edges.
        
        Notes
        -----
        This is an internal validation method called by ``validate()``.
        """
        # Group edges by (u, v) to find parallel edges
        edge_groups: dict[tuple[T, T], list[tuple[int, dict[str, Any]]]] = {}
        
        for u, v, key, data in self._graph.edges_iter(keys=True, data=True):
            edge_key = (u, v)
            if edge_key not in edge_groups:
                edge_groups[edge_key] = []
            edge_groups[edge_key].append((key, data or {}))
        
        # Check each group of parallel edges
        for (u, v), edges in edge_groups.items():
            if len(edges) <= 1:
                continue  # No parallel edges, skip
            
            # Collect branch_length values
            branch_lengths: list[float] = []
            missing_branch_lengths: list[int] = []
            
            for key, data in edges:
                bl = data.get('branch_length')
                if bl is None:
                    missing_branch_lengths.append(key)
                else:
                    if not isinstance(bl, (int, float)):
                        raise PhyloZooNetworkAttributeError(
                            f"Edge ({u}, {v}, key={key}) has branch_length of type {type(bl).__name__}, "
                            f"but must be numeric."
                        )
                    branch_lengths.append(float(bl))
            
            # Check constraint: if any edge has branch_length, all must have it
            if branch_lengths and missing_branch_lengths:
                missing_keys_str = ', '.join(str(k) for k in missing_branch_lengths)
                raise PhyloZooNetworkAttributeError(
                    f"Parallel edges between {u} and {v} have inconsistent branch_length attributes. "
                    f"Some edges have branch_length (keys: {[k for k, d in edges if d.get('branch_length') is not None]}), "
                    f"but others do not (keys: {missing_keys_str}). "
                    f"If one parallel edge has branch_length, all must have branch_length."
                )
            
            # Check constraint: all branch_length values must be the same
            if len(branch_lengths) > 1:
                first_bl = branch_lengths[0]
                for i, bl in enumerate(branch_lengths[1:], start=1):
                    if abs(bl - first_bl) > 1e-10:
                        keys_with_bl = [k for k, d in edges if d.get('branch_length') is not None]
                        raise PhyloZooNetworkAttributeError(
                            f"Parallel edges between {u} and {v} have different branch_length values. "
                            f"All parallel edges must have the same branch_length. "
                            f"Found values: {branch_lengths} for keys: {keys_with_bl}"
                        )
    
    def validate(self) -> None:
        """
        Validate the network structure and edge attributes.
        
        Checks:
        1. Network is connected (weakly connected)
        2. No self-loops
        3. Directed acyclic graph (no directed cycles)
        4. Single root node (exactly one node with in-degree 0)
        5. Leaf nodes: all have in-degree 1 and out-degree 0
        6. Internal nodes: all have either (in-degree 1 and out-degree >= 2) or
           (in-degree >= 2 and out-degree 1)
        7. Bootstrap values: all bootstrap values must be in [0.0, 1.0]
        8. Gamma constraints: gamma can only be set on hybrid edges, and if any gamma
           is specified for a hybrid node, all incoming edges must have gamma values
           summing to 1.0
        9. Branch length constraints: for each set of parallel edges, if one edge has
           branch_length, all must have branch_length, and all values must be the same
        
        Raises
        ------
        PhyloZooNetworkStructureError
            If connectivity, self-loop, or acyclicity constraints are violated.
        PhyloZooNetworkDegreeError
            If degree constraints are violated.
        PhyloZooNetworkAttributeError
            If bootstrap, gamma, or branch length constraints are violated.
        PhyloZooEmptyNetworkWarning
            If empty network is detected.
        PhyloZooSingleNodeNetworkWarning
            If single-node network is detected.
        Notes
        -----
        This method performs comprehensive validation of the network structure
        and edge attributes. See class docstring for detailed validation rules.
        Empty networks (no nodes) are considered valid but will raise a warning.
        Single-node networks (where root and leaf are the same node) are also valid
        but will raise a warning.
        """
        # Empty networks are valid
        if self.number_of_nodes() == 0:
            warnings.warn(
                "Empty network (no nodes) detected. While valid, this may not be useful for phylogenetic analysis.",
                PhyloZooEmptyNetworkWarning,
                stacklevel=2
            )
            return

        # Single-node networks are valid only if they have no self-loops
        if self.number_of_nodes() == 1:
            if has_self_loops(self._graph):
                raise PhyloZooNetworkStructureError("Self-loops are not allowed in DirectedPhyNetwork.")
            warnings.warn(
                "Single-node network detected. While valid, this may not be useful for phylogenetic analysis.",
                PhyloZooSingleNodeNetworkWarning,
                stacklevel=2
            )
            return

        self._validate_structural_constraints()
        
        # 3. Validate degree constraints (includes root check)
        self._validate_degree_constraints()
        
        # 4. Validate bootstrap constraints
        self._validate_bootstrap_constraints()
        
        # 5. Validate gamma constraints
        self._validate_gamma_constraints()
        
        # 6. Validate branch length constraints
        self._validate_branchlength_constraints()
    
    def _validate_structural_constraints(self) -> None:
        """
        Validate structural constraints (emptiness, connectivity, acyclicity, self-loops).
        
        Raises
        ------
        PhyloZooNetworkStructureError
            If connectivity, self-loop, or acyclicity constraints are violated.
        """
        # 1. Check that network is connected (weakly connected)
        if not is_connected(self._graph):
            raise PhyloZooNetworkStructureError(
                "Network is not connected. All nodes must be in a single weakly connected component."
            )
        
        # 2. Disallow self-loops
        if has_self_loops(self._graph):
            raise PhyloZooNetworkStructureError("Self-loops are not allowed in DirectedPhyNetwork.")
        
        # 3. Check for directed cycles (must be acyclic)
        if not nx.is_directed_acyclic_graph(self._graph._graph):
            cycles = list(nx.simple_cycles(self._graph._graph))
            raise PhyloZooNetworkStructureError(
                f"Network contains directed cycles. Found {len(cycles)} cycle(s). "
                f"First cycle: {cycles[0] if cycles else 'unknown'}"
            )
    
    # ========== Label Operations ==========
    
    def get_label(self, node_id: T) -> str | None:
        """
        Get the label for a node.
        
        Parameters
        ----------
        node_id : T
            Node identifier.
        
        Returns
        -------
        str | None
            Label for the node. Returns None if node has no label.
            Leaves always have labels (taxa), but internal nodes may be unlabeled.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(edges=[(3, 1)], nodes=[(1, {'label': 'A'})])
        >>> net.get_label(1)
        'A'
        >>> net.get_label(3) is None
        True
        """
        return self._node_to_label.get(node_id)
    
    def get_node_attribute(self, node_id: T, attr: str | None = None) -> dict[str, Any] | Any | None:
        """
        Get node attribute(s).
        
        Parameters
        ----------
        node_id : T
            Node identifier.
        attr : str | None, optional
            Attribute name. If None, returns all attributes as a dict.
            If specified, returns the value of that specific attribute.
            By default None.
        
        Returns
        -------
        dict[str, Any] | Any | None
            If attr is None: dict of all node attributes (empty dict if no attributes).
            If attr is specified: attribute value, or None if not set.
        
        Raises
        ------
        PhyloZooValueError
            If the node does not exist in the network.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(
        ...     edges=[(3, 1)],
        ...     nodes=[(1, {'label': 'A'}), (3, {'label': 'root', 'custom': 42})]
        ... )
        >>> net.get_node_attribute(1, 'label')
        'A'
        >>> net.get_node_attribute(3, 'label')
        'root'
        >>> net.get_node_attribute(3, 'custom')
        42
        >>> net.get_node_attribute(1, 'nonexistent') is None
        True
        >>> net.get_node_attribute(1)  # Get all attributes
        {'label': 'A'}
        >>> net.get_node_attribute(3)  # Get all attributes
        {'label': 'root', 'custom': 42}
        """
        if node_id not in self._graph:
            raise PhyloZooValueError(f"Node {node_id} does not exist in the network.")
        
        if attr is None:
            return self._graph._graph.nodes[node_id].copy()
        else:
            return self._graph._graph.nodes[node_id].get(attr)
    
    def get_node_id(self, label: str) -> T | None:
        """
        Get the node ID for a label.
        
        Parameters
        ----------
        label : str
            Node label.
        
        Returns
        -------
        T | None
            Node ID if found, None otherwise.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(edges=[(3, 1)], nodes=[(1, {'label': 'A'})])
        >>> net.get_node_id("A")
        1
        """
        return self._label_to_node.get(label)
    
    # ========== Network Attribute Access (Read-Only) ==========
    
    def get_network_attribute(self, key: str | None = None) -> dict[str, Any] | Any | None:
        """
        Get network-level attribute(s).
        
        Parameters
        ----------
        key : str | None, optional
            Attribute key. If None, returns all attributes as a dict.
            If specified, returns the value of that specific attribute.
            By default None.
        
        Returns
        -------
        dict[str, Any] | Any | None
            If key is None: dict of all network attributes (empty dict if no attributes).
            If key is specified: attribute value, or None if not set.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(
        ...     edges=[(3, 1)],
        ...     nodes=[(1, {'label': 'A'})],
        ...     attributes={'source': 'file.nex', 'version': '1.0'}
        ... )
        >>> net.get_network_attribute('source')
        'file.nex'
        >>> net.get_network_attribute('nonexistent') is None
        True
        >>> net.get_network_attribute()  # Get all attributes
        {'source': 'file.nex', 'version': '1.0'}
        """
        if key is None:
            return self._graph._graph.graph.copy()
        else:
            return self._graph._graph.graph.get(key)
    
    # ========== Edge Attribute Access (Read-Only) ==========
    
    def get_edge_attribute(self, u: T, v: T, key: int | None = None, attr: str | None = None) -> dict[str, Any] | Any | None:
        """
        Get edge attribute(s).
        
        Parameters
        ----------
        u, v : T
            Edge endpoints.
        key : int | None, optional
            Edge key for parallel edges. If None and multiple parallel edges exist,
            raises ValueError. Must specify key when parallel edges exist.
        attr : str | None, optional
            Attribute name. If None, returns all attributes as a dict.
            If specified, returns the value of that specific attribute.
            By default None.
        
        Returns
        -------
        dict[str, Any] | Any | None
            If attr is None: dict of all edge attributes (empty dict if no attributes).
            If attr is specified: attribute value, or None if not set.
        
        Raises
        ------
        PhyloZooValueError
            If the edge does not exist, or if key is None and multiple parallel edges exist.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(
        ...     edges=[
        ...         {'u': 3, 'v': 2, 'key': 0, 'branch_length': 0.5, 'bootstrap': 0.95},
        ...         {'u': 3, 'v': 2, 'key': 1, 'branch_length': 0.7},  # Parallel edge
        ...         (2, 1)  # Tree node 2 to leaf
        ...     ],
        ...     nodes=[(1, {'label': 'A'})]
        ... )
        >>> net.get_edge_attribute(3, 2, key=0, attr='branch_length')
        0.5
        >>> net.get_edge_attribute(3, 2, key=1, attr='branch_length')
        0.7
        >>> net.get_edge_attribute(3, 2, key=0)  # Get all attributes
        {'branch_length': 0.5, 'bootstrap': 0.95}
        >>> net.get_edge_attribute(2, 1)  # Get all attributes
        {}
        """
        if not self.has_edge(u, v):
            return {} if attr is None else None
        
        # Check if parallel edges exist
        edges_data = self._graph._graph[u][v]
        num_edges = len(edges_data)
        
        if num_edges == 0:
            return {} if attr is None else None
        elif num_edges == 1:
            # Single edge - key not needed
            first_key = next(iter(edges_data))
            edge_attrs = edges_data[first_key]
            if attr is None:
                return edge_attrs.copy()
            else:
                return edge_attrs.get(attr)
        else:
            # Multiple parallel edges - key is required
            if key is None:
                raise PhyloZooValueError(
                    f"Multiple parallel edges exist between {u} and {v}. "
                    f"Must specify 'key' parameter to get attributes from a specific edge."
                )
            if key not in edges_data:
                raise PhyloZooValueError(f"Edge ({u}, {v}, {key}) does not exist in the network.")
            edge_attrs = edges_data[key]
            if attr is None:
                return edge_attrs.copy()
            else:
                return edge_attrs.get(attr)
    
    def get_branch_length(self, u: T, v: T, key: int | None = None) -> float | None:
        """
        Get branch length for an edge.
        
        Parameters
        ----------
        u, v : T
            Edge endpoints.
        key : int | None, optional
            Edge key for parallel edges. Required if multiple parallel edges exist.
        
        Returns
        -------
        float | None
            Branch length, or None if not set.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(
        ...     edges=[{'u': 3, 'v': 1, 'branch_length': 0.5}],
        ...     nodes=[(1, {'label': 'A'})]
        ... )
        >>> net.get_branch_length(3, 1)
        0.5
        """
        return self.get_edge_attribute(u, v, key, 'branch_length')
    
    def get_bootstrap(self, u: T, v: T, key: int | None = None) -> float | None:
        """
        Get bootstrap support for an edge.
        
        Bootstrap values are typically in the range 0.0 to 1.0.
        
        Parameters
        ----------
        u, v : T
            Edge endpoints.
        key : int | None, optional
            Edge key for parallel edges. Required if multiple parallel edges exist.
        
        Returns
        -------
        float | None
            Bootstrap support value (typically 0.0 to 1.0), or None if not set.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(
        ...     edges=[{'u': 3, 'v': 1, 'bootstrap': 0.95}],
        ...     nodes=[(1, {'label': 'A'})]
        ... )
        >>> net.get_bootstrap(3, 1)
        0.95
        """
        return self.get_edge_attribute(u, v, key, 'bootstrap')
    
    def get_gamma(self, u: T, v: T, key: int | None = None) -> float | None:
        """
        Get gamma value for a hybrid edge.
        
        Gamma values can only be set on hybrid edges (edges pointing into hybrid nodes).
        If ANY gamma value is specified for edges entering a hybrid node, then
        ALL edges entering that hybrid node must have gamma values, and they must
        sum to exactly 1.0.
        
        Parameters
        ----------
        u, v : T
            Edge endpoints (v must be a hybrid node).
        key : int | None, optional
            Edge key for parallel edges. Required if multiple parallel edges exist.
        
        Returns
        -------
        float | None
            Gamma value, or None if not set.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(
        ...     edges=[
        ...         (7, 5), (7, 6),  # Root to tree nodes
        ...         {'u': 5, 'v': 4, 'gamma': 0.6},  # Hybrid edge
        ...         {'u': 6, 'v': 4, 'gamma': 0.4},  # Hybrid edge (Sum = 1.0)
        ...         (5, 8), (6, 9),  # Tree nodes also have other children
        ...         (4, 1)  # Hybrid to leaf
        ...     ],
        ...     nodes=[(1, {'label': 'A'}), (8, {'label': 'B'}), (9, {'label': 'C'})]
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
    
    def has_edge(self, u: T, v: T, key: int | None = None) -> bool:
        """
        Check if edge exists.
        
        Parameters
        ----------
        u : T
            Source node.
        v : T
            Target node.
        key : int | None, optional
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
    
    def incident_parent_edges(self, v: T, keys: bool = False, data: bool = False) -> Iterator[tuple[T, T] | tuple[T, T, int] | tuple[T, T, int, dict[str, Any]]]:
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
        ...     nodes=[(2, {'label': 'A'})]
        ... )
        >>> list(net.incident_parent_edges(2))
        [(1, 2), (3, 2)]
        >>> list(net.incident_parent_edges(2, data=True))
        [(1, 2, {'branch_length': 0.5}), (3, 2, {'branch_length': 0.3})]
        """
        return self._graph.incident_parent_edges(v, keys=keys, data=data)
    
    def incident_child_edges(self, v: T, keys: bool = False, data: bool = False) -> Iterator[tuple[T, T] | tuple[T, T, int] | tuple[T, T, int, dict[str, Any]]]:
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
        ...     nodes=[(2, {'label': 'A'}), (3, {'label': 'B'})]
        ... )
        >>> list(net.incident_child_edges(1))
        [(1, 2), (1, 3)]
        >>> list(net.incident_child_edges(1, data=True))
        [(1, 2, {'branch_length': 0.5}), (1, 3, {'branch_length': 0.3})]
        """
        return self._graph.incident_child_edges(v, keys=keys, data=data)
    
    # ========== Phylogenetic-Specific Methods ==========
    
    @cached_property
    def leaves(self) -> set[T]:
        """
        Get the set of leaf node IDs (nodes with no outgoing edges).
        
        Returns
        -------
        set[T]
            Set of leaf node identifiers. Returns a new set (which is mutable).
        """
        return {node for node in self._graph.nodes if self._graph.outdegree(node) == 0}
    
    @cached_property
    def taxa(self) -> set[str]:
        """
        Get the set of taxon labels (labels of leaves).
        
        Returns
        -------
        set[str]
            Set of taxon labels.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        >>> net.taxa
        {'A', 'B'}
        """
        return {self._node_to_label[leaf] for leaf in self.leaves}
    
    @cached_property
    def internal_nodes(self) -> set[T]:
        """
        Get the set of all internal (non-root, non-leaf) nodes.
        
        Returns
        -------
        set[T]
            Set of internal node identifiers. Returns a new set (which is mutable).
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        >>> net.internal_nodes
        {3}
        """
        # Handle empty network
        if len(self._graph) == 0:
            return set()
        root = self.root_node
        leaves = self.leaves
        return {v for v in self._graph.nodes if v != root and v not in leaves}
    
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
        PhyloZooValueError
            If there is no root node or multiple root nodes.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        >>> net.root_node
        3
        """
        roots = [v for v in self._graph.nodes if self._graph.indegree(v) == 0]
        if len(roots) == 0:
            raise PhyloZooValueError("Network has no root node")
        if len(roots) > 1:
            raise PhyloZooValueError(f"Network has multiple root nodes: {roots}")
        return roots[0]
    
    @cached_property
    def hybrid_nodes(self) -> set[T]:
        """
        Get the set of all hybrid nodes.
        
        A hybrid node is a node with in-degree >= 2 and out-degree 1.
        
        Returns
        -------
        set[T]
            Set of hybrid node identifiers. Returns a new set (which is mutable).
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(edges=[(5, 4), (6, 4), (4, 1)], nodes=[(1, {'label': 'A'})])
        >>> net.hybrid_nodes
        {4}
        """
        return {
            v for v in self._graph.nodes
            if self._graph.indegree(v) >= 2 and self._graph.outdegree(v) == 1
        }
    
    @cached_property
    def LSA_node(self) -> T:
        """
        Return the Least Stable Ancestor (LSA) node of the network.
        
        The LSA is the lowest node through which all paths from the root to the leaves pass.
        In other words, it is the unique node that is an ancestor of all leaves and is
        the lowest such node (has maximum depth from the root).
        
        Returns
        -------
        T
            The LSA node identifier.
        
        Raises
        ------
        PhyloZooValueError
            If the network is empty.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        >>> net.LSA_node
        3
        """
        from .features import lsa_node
        return lsa_node(self)
    
    @cached_property
    def tree_nodes(self) -> set[T]:
        """
        Get the set of all tree nodes.
        
        A tree node is an internal node (non-root, non-leaf) with in-degree 1
        and out-degree >= 2.
        
        Returns
        -------
        set[T]
            Set of tree node identifiers. Returns a new set (which is mutable).
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        >>> net.tree_nodes
        {3}
        """
        # Handle empty network
        if len(self._graph) == 0:
            return set()
        root = self.root_node
        leaves = self.leaves
        return {
            v for v in self._graph.nodes
            if v != root
            and v not in leaves
            and self._graph.indegree(v) == 1
            and self._graph.outdegree(v) >= 2
        }
    
    @cached_property
    def hybrid_edges(self) -> set[tuple[T, T, int]]:
        """
        Get the set of all hybrid edges with keys.
        
        Hybrid edges are edges that point into hybrid nodes.
        
        Returns
        -------
        set[tuple[T, T, int]]
            Set of (source, target, key) tuples for hybrid edges. Returns a new set (which is mutable).
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(edges=[(3, 2), (4, 2)], nodes=[(2, {'label': 'A'})])
        >>> net.hybrid_edges
        {(3, 2, 0), (4, 2, 0)}
        """
        res = set()
        for v in self.hybrid_nodes:
            for p, _, key in self._graph.incident_parent_edges(v, keys=True):
                res.add((p, v, key))
        return res
    
    @cached_property
    def tree_edges(self) -> set[tuple[T, T, int]]:
        """
        Get the set of all tree edges with keys.
        
        Tree edges are all edges that are not hybrid edges.
        
        Returns
        -------
        set[tuple[T, T, int]]
            Set of (source, target, key) tuples for tree edges. Returns a new set (which is mutable).
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        >>> len(net.tree_edges)
        2
        """
        hybrid_edges = self.hybrid_edges
        # Get all edges with keys
        return {(u, v, key) for u, v, key in self._graph.edges(keys=True) if (u, v, key) not in hybrid_edges}
    
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
        >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
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
        >>> net = DirectedPhyNetwork(edges=[(3, 1)], nodes=[(1, {'label': 'A'})])
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
            String representation showing nodes, edges, taxa count, and taxon list.
        
        Examples
        --------
        >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        >>> repr(net)
        'DirectedPhyNetwork(nodes=3, edges=2, taxa=2, taxa_list=[A, B])'
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
            f"taxa={n_taxa}, "
            f"taxa_list=[{taxa_list_str}])"
        )
    
    def __len__(self) -> int:
        """
        Return the number of nodes in the network.
        
        Returns
        -------
        int
            Number of nodes.
        """
        return self.number_of_nodes()
    
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

