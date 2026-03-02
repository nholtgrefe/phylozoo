"""
Mixed network module.

This module provides classes and functions for working with mixed phylogenetic networks.
"""

import math
import warnings
from functools import cached_property
from typing import Any, Iterator, TypeVar

from ....utils.exceptions import (
    PhyloZooNetworkDegreeError,
    PhyloZooNetworkStructureError,
    PhyloZooNetworkAttributeError,
    PhyloZooValueError,
    PhyloZooTypeError,
    PhyloZooEmptyNetworkWarning,
    PhyloZooSingleNodeNetworkWarning,
)

from ...primitives.m_multigraph import MixedMultiGraph
from ...primitives.m_multigraph.features import is_connected, has_self_loops
from ....utils.validation import validation_aware

T = TypeVar('T')


@validation_aware(allowed=["validate", "_validate_*"], default=["validate"])
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
    directed_edges : list[tuple[T, T] | tuple[T, T, int] | dict[str, Any]] | None, optional
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
        attributes.
        
        Can be empty or None for empty/single-node networks. By default None.
    undirected_edges : list[tuple[T, T] | tuple[T, T, int] | dict[str, Any]] | None, optional
        List of undirected edges. Formats:
        - (u, v) tuples (key auto-generated)
        - (u, v, key) tuples (explicit key)
        - Dict with 'u', 'v' and optional 'key' (for parallel edges) plus edge attributes
        
        Edge attributes (validated):
        - branch_length (float; for set of parallel edges, all must have equal branch_length)
        - bootstrap (float in [0.0, 1.0])
        
        Note: Undirected edges cannot have gamma values.
        
        Can be empty or None for empty/single-node networks. By default None.
    nodes : list[T | tuple[T | dict[str | Any | None]]], optional
        List of nodes. Formats:
        - Simple node IDs: `1`, `"node1"`, etc.
        - Tuples: `(node_id, {'label': '...','attr': ...})`
        
        Node attributes (validated):
        - label: string, unique across all nodes; use another key for non-string data.
        
        Leaves without labels are auto-labeled. Leaf-labels are referred to as `taxa`.
        Use a different attribute name (e.g., 'label2') for non-validated and/or additional
        attributes.

        Can be empty or None. By default None.
    attributes : dict[str, Any] | None, optional
        Optional dictionary of graph-level attributes to store with the network.
        These attributes are stored in the underlying graph's `.graph` attribute
        and are preserved through copy operations. Can be used to store metadata
        like provenance, source file, creation date, etc.
        By default None.
    
    Attributes
    ----------
    nodes
        Cached view of all node IDs (same interface as the underlying graph's ``nodes``).
    edges
        Cached view of all edges (directed and undirected; same interface as the underlying graph's ``edges``).
    _graph : MixedMultiGraph[T]
        Internal graph structure using MixedMultiGraph.
        **Warning:** Do not modify directly. Use class methods instead.
    _node_to_label : dict[T, str]
        Mapping from node IDs to labels. Only nodes with explicit labels are included.
        Leaves always have labels (taxa), but internal nodes may be unlabeled.
    _label_to_node : dict[str, T]
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
        directed_edges: list[tuple[T, T] | tuple[T, T, int] | dict[str, Any]] | None = None,
        undirected_edges: list[tuple[T, T] | tuple[T, T, int] | dict[str, Any]] | None = None,
        nodes: list[T | tuple[T | dict[str | Any | None]]] = None,
        attributes: dict[str, Any] | None = None,
    ) -> None:
        """
        Initialize a mixed phylogenetic network.
        
        Parameters
        ----------
        directed_edges : list[tuple[T, T] | tuple[T, T, int] | dict[str, Any]] | None, optional
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
        undirected_edges : list[tuple[T, T] | tuple[T, T, int] | dict[str, Any]] | None, optional
            List of undirected edges. Formats:
            - (u, v) tuples (key auto-generated)
            - (u, v, key) tuples (explicit key)
            - Dict with 'u', 'v' and optional 'key' plus edge attributes
            
            Edge attributes (validated):
            - branch_length (float)
            - bootstrap (float in [0.0, 1.0])
            
            Can be empty list or None for empty/single-node networks. By default None.
        nodes : list[T | tuple[T | dict[str | Any | None]]], optional
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
        
        self._graph: MixedMultiGraph[T] = MixedMultiGraph(
            directed_edges=directed_edges,
            undirected_edges=undirected_edges,
            attributes=attributes
        )
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
        nodes: list[T | tuple[T | dict[str | Any | None]]],
    ) -> None:
        """
        Add all nodes to the underlying graph with their attributes.
        
        If a node has a 'label' attribute, it is immediately validated and added to
        the label dictionaries using _add_label_to_dicts. Label validation (string
        type and uniqueness) is performed during this step.
        
        Parameters
        ----------
        nodes : list[T | tuple[T | dict[str | Any | None]]]
            Node specifications as simple IDs or (node_id, attr_dict) tuples.
        
        Raises
        ------
        PhyloZooTypeError
            If a node tuple does not provide a dict of attributes.
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
        PhyloZooNetworkAttributeError
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
                    raise PhyloZooNetworkAttributeError(
                        f"Bootstrap value on directed edge ({u}, {v}, key={key}) must be numeric, "
                        f"got {type(bootstrap).__name__}"
                    )
                if math.isnan(bootstrap) or bootstrap < 0.0 or bootstrap > 1.0:
                    raise PhyloZooNetworkAttributeError(
                        f"Bootstrap value on directed edge ({u}, {v}, key={key}) is {bootstrap}, "
                        f"but must be in [0.0, 1.0]"
                    )
        
        # Check undirected edges
        for u, v, key, data in self._graph._undirected.edges(keys=True, data=True):
            if 'bootstrap' in data:
                bootstrap = data['bootstrap']
                if not isinstance(bootstrap, (int, float)):
                    raise PhyloZooNetworkAttributeError(
                        f"Bootstrap value on undirected edge ({u}, {v}, key={key}) must be numeric, "
                        f"got {type(bootstrap).__name__}"
                    )
                if math.isnan(bootstrap) or bootstrap < 0.0 or bootstrap > 1.0:
                    raise PhyloZooNetworkAttributeError(
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
        PhyloZooNetworkAttributeError
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
        hybrid_edges_set = self.hybrid_edges  # Now contains (u, v, key) tuples
        
        # Check directed edges
        for u, v, key, data in self._graph._directed.edges(keys=True, data=True):
            if 'gamma' in data:
                # Check if this edge is a hybrid edge
                if (u, v, key) not in hybrid_edges_set:
                    raise PhyloZooNetworkAttributeError(
                        f"Gamma value can only be set on hybrid edges (edges pointing into "
                        f"hybrid nodes). Directed edge ({u}, {v}, key={key}) is not a hybrid edge."
                    )
        
        # Check that gamma is not set on undirected edges
        for u, v, key, data in self._graph._undirected.edges(keys=True, data=True):
            if 'gamma' in data:
                raise PhyloZooNetworkAttributeError(
                    f"Gamma values cannot be set on undirected edges. "
                    f"Undirected edge ({u}, {v}, key={key}) cannot have gamma values."
                )
        
        # Then validate gamma constraints for hybrid nodes
        for hybrid_node in self.hybrid_nodes:
            gamma_values: list[float] = []
            incoming_edges: list[tuple[T, T, int]] = []
            
            # Use incident_parent_edges to get all incoming directed edges (including parallel edges)
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
        PhyloZooNetworkDegreeError
            If any degree constraints are violated.
        """
        # 1. Check that all internal nodes have degree >= 3
        internal_nodes = self.internal_nodes
        for node in internal_nodes:
            degree = self._graph.degree(node)
            if degree < 3:
                raise PhyloZooNetworkDegreeError(
                    f"Internal node {node} has degree {degree}, but all internal nodes "
                    f"must have degree >= 3."
                )
        
        # 2. Check that each node has indegree either 0 or total_degree-1
        # This constraint applies to all nodes (including leaves)
        for node in self._graph.nodes:
            indegree = self._graph.indegree(node)
            total_degree = self._graph.degree(node)
            if indegree != 0 and indegree != total_degree - 1:
                raise PhyloZooNetworkDegreeError(
                    f"Node {node} has indegree {indegree} and total degree {total_degree}. "
                    f"Each node must have indegree either 0 or total_degree-1."
                )
    
    def _validate_branchlength_constraints(self) -> None:
        """
        Validate branch length constraints for parallel edges.
        
        For each set of parallel edges between nodes u and v (for both directed
        and undirected edges separately), this method ensures:
        1. If one edge has a branch_length attribute, all parallel edges must have branch_length
        2. All branch_length values must be the same across parallel edges
        
        Raises
        ------
        PhyloZooNetworkAttributeError
            If branch length constraints are violated for any set of parallel edges.
        
        Notes
        -----
        This is an internal validation method called by ``validate()``.
        """
        # Check directed edges
        directed_edge_groups: dict[tuple[T, T], list[tuple[int, dict[str, Any]]]] = {}
        
        for u, v, key, data in self._graph.directed_edges_iter(keys=True, data=True):
            edge_key = (u, v)
            if edge_key not in directed_edge_groups:
                directed_edge_groups[edge_key] = []
            directed_edge_groups[edge_key].append((key, data or {}))
        
        # Check undirected edges (normalize for consistency)
        undirected_edge_groups: dict[tuple[T, T], list[tuple[int, dict[str, Any]]]] = {}
        
        for u, v, key, data in self._graph.undirected_edges_iter(keys=True, data=True):
            edge_key = self._graph.normalize_undirected_edge(u, v)
            if edge_key not in undirected_edge_groups:
                undirected_edge_groups[edge_key] = []
            undirected_edge_groups[edge_key].append((key, data or {}))
        
        # Validate directed edges
        for (u, v), edges in directed_edge_groups.items():
            if len(edges) <= 1:
                continue  # No parallel edges, skip
            
            branch_lengths: list[float] = []
            missing_branch_lengths: list[int] = []
            
            for key, data in edges:
                bl = data.get('branch_length')
                if bl is None:
                    missing_branch_lengths.append(key)
                else:
                    if not isinstance(bl, (int, float)):
                        raise PhyloZooNetworkAttributeError(
                            f"Directed edge ({u}, {v}, key={key}) has branch_length of type {type(bl).__name__}, "
                            f"but must be numeric."
                        )
                    branch_lengths.append(float(bl))
            
            if branch_lengths and missing_branch_lengths:
                missing_keys_str = ', '.join(str(k) for k in missing_branch_lengths)
                raise PhyloZooNetworkAttributeError(
                    f"Parallel directed edges between {u} and {v} have inconsistent branch_length attributes. "
                    f"Some edges have branch_length (keys: {[k for k, d in edges if d.get('branch_length') is not None]}), "
                    f"but others do not (keys: {missing_keys_str}). "
                    f"If one parallel edge has branch_length, all must have branch_length."
                )
            
            if len(branch_lengths) > 1:
                first_bl = branch_lengths[0]
                for i, bl in enumerate(branch_lengths[1:], start=1):
                    if abs(bl - first_bl) > 1e-10:
                        keys_with_bl = [k for k, d in edges if d.get('branch_length') is not None]
                        raise PhyloZooNetworkAttributeError(
                            f"Parallel directed edges between {u} and {v} have different branch_length values. "
                            f"All parallel edges must have the same branch_length. "
                            f"Found values: {branch_lengths} for keys: {keys_with_bl}"
                        )
        
        # Validate undirected edges
        for (u, v), edges in undirected_edge_groups.items():
            if len(edges) <= 1:
                continue  # No parallel edges, skip
            
            branch_lengths: list[float] = []
            missing_branch_lengths: list[int] = []
            
            for key, data in edges:
                bl = data.get('branch_length')
                if bl is None:
                    missing_branch_lengths.append(key)
                else:
                    if not isinstance(bl, (int, float)):
                        raise PhyloZooNetworkAttributeError(
                            f"Undirected edge ({u}, {v}, key={key}) has branch_length of type {type(bl).__name__}, "
                            f"but must be numeric."
                        )
                    branch_lengths.append(float(bl))
            
            if branch_lengths and missing_branch_lengths:
                missing_keys_str = ', '.join(str(k) for k in missing_branch_lengths)
                raise PhyloZooNetworkAttributeError(
                    f"Parallel undirected edges between {u} and {v} have inconsistent branch_length attributes. "
                    f"Some edges have branch_length (keys: {[k for k, d in edges if d.get('branch_length') is not None]}), "
                    f"but others do not (keys: {missing_keys_str}). "
                    f"If one parallel edge has branch_length, all must have branch_length."
                )
            
            if len(branch_lengths) > 1:
                first_bl = branch_lengths[0]
                for i, bl in enumerate(branch_lengths[1:], start=1):
                    if abs(bl - first_bl) > 1e-10:
                        keys_with_bl = [k for k, d in edges if d.get('branch_length') is not None]
                        raise PhyloZooNetworkAttributeError(
                            f"Parallel undirected edges between {u} and {v} have different branch_length values. "
                            f"All parallel edges must have the same branch_length. "
                            f"Found values: {branch_lengths} for keys: {keys_with_bl}"
                        )
    
    def validate(self) -> None:
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
        8. Branch length constraints: for each set of parallel edges, if one edge has
           branch_length, all must have branch_length, and all values must be the same
        
        Raises
        ------
        PhyloZooNetworkStructureError
            If connectivity or self-loop constraints are violated.
        PhyloZooNetworkDegreeError
            If degree constraints are violated.
        PhyloZooNetworkAttributeError
            If bootstrap, gamma, or branch length constraints are violated.
        PhyloZooEmptyNetworkWarning
            If empty network is detected.
        PhyloZooSingleNodeNetworkWarning
            If single-node network is detected.
        
        Warns
        -----
        UserWarning
            Always issued to indicate that additional validation checks may be added later.
        
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
                PhyloZooEmptyNetworkWarning,
                stacklevel=2
            )
            return

        # Single-node networks are valid only if they have no self-loops
        if self.number_of_nodes() == 1:
            if has_self_loops(self._graph):
                raise PhyloZooNetworkStructureError("Self-loops are not allowed in MixedPhyNetwork.")
            warnings.warn(
                "Single-node network detected. While valid, this may not be useful for phylogenetic analysis.",
                PhyloZooSingleNodeNetworkWarning,
                stacklevel=2
            )
            return

        self._validate_structural_constraints()
        
        # 3. Validate degree constraints
        self._validate_degree_constraints()
        
        # 4. Validate mixed network constraint (issues warning)
        self._validate_mixednetwork_constraint()
        
        # 5. Validate bootstrap constraints
        self._validate_bootstrap_constraints()
        
        # 6. Validate gamma constraints
        self._validate_gamma_constraints()
        
        # 7. Validate branch length constraints
        self._validate_branchlength_constraints()
    
    def _validate_structural_constraints(self) -> None:
        """
        Validate structural constraints (emptiness, connectivity, self-loops).
        
        Raises
        ------
        PhyloZooNetworkStructureError
            If connectivity or self-loop constraints are violated.
        """
        # 1. Check that network is connected (weakly connected)
        if not is_connected(self._graph):
            raise PhyloZooNetworkStructureError(
                "Network is not connected. All nodes must be in a single connected component."
            )
        
        # 2. Disallow self-loops
        if has_self_loops(self._graph):
            raise PhyloZooNetworkStructureError("Self-loops are not allowed in MixedPhyNetwork.")
    
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
        >>> net = MixedPhyNetwork(undirected_edges=[(3, 1)], nodes=[(1, {'label': 'A'})])
        >>> net.get_label(1)
        'A'
        >>> net.get_label(3) is None
        True
        """
        return self._node_to_label.get(node_id)
    
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
        >>> net = MixedPhyNetwork(undirected_edges=[(3, 1)], nodes=[(1, {'label': 'A'})])
        >>> net.get_node_id("A")
        1
        """
        return self._label_to_node.get(label)
    
    # ========== Node Attribute Access (Read-Only) ==========
    
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
        >>> net = MixedPhyNetwork(
        ...     undirected_edges=[(3, 1)],
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
        
        # Use undirected graph for node attributes (nodes exist in both)
        if attr is None:
            return self._graph._undirected.nodes[node_id].copy()
        else:
            return self._graph._undirected.nodes[node_id].get(attr)
    
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
        >>> net = MixedPhyNetwork(
        ...     undirected_edges=[(3, 1)],
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
            return self._graph._directed.graph.copy()
        else:
            return self._graph._directed.graph.get(key)
    
    # ========== Edge Attribute Access (Read-Only) ==========
    
    def get_edge_attribute(
        self,
        u: T,
        v: T,
        key: int | None = None,
        attr: str | None = None
    ) -> dict[str, Any] | Any | None:
        """
        Get edge attribute(s).
        
        Parameters
        ----------
        u, v : T
            Edge endpoints.
        key : int | None, optional
            Edge key for parallel edges. If None and multiple parallel edges exist,
            raises PhyloZooValueError. Must specify key when parallel edges exist.
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
                return {} if attr is None else None
            
            num_edges = len(edges_data)
            if num_edges == 0:
                return {} if attr is None else None
            elif num_edges == 1:
                first_key = next(iter(edges_data))
                edge_attrs = edges_data[first_key]
                if attr is None:
                    return edge_attrs.copy()
                else:
                    return edge_attrs.get(attr)
            else:
                if key is None:
                    raise PhyloZooValueError(
                        f"Multiple parallel directed edges exist between {u} and {v}. "
                        f"Must specify 'key' parameter to get attributes from a specific edge."
                    )
                if key not in edges_data:
                    raise PhyloZooValueError(f"Edge ({u}, {v}, {key}) does not exist in the network.")
                edge_attrs = edges_data[key]
                if attr is None:
                    return edge_attrs.copy()
                else:
                    return edge_attrs.get(attr)
        
        # Try undirected edges
        if self._graph._undirected.has_edge(u, v, key=key):
            # Get all undirected edges between u and v
            edges_data = self._graph._undirected[u].get(v, {})
            if not edges_data:
                return {} if attr is None else None
            
            num_edges = len(edges_data)
            if num_edges == 0:
                return {} if attr is None else None
            elif num_edges == 1:
                first_key = next(iter(edges_data))
                edge_attrs = edges_data[first_key]
                if attr is None:
                    return edge_attrs.copy()
                else:
                    return edge_attrs.get(attr)
            else:
                if key is None:
                    raise PhyloZooValueError(
                        f"Multiple parallel undirected edges exist between {u} and {v}. "
                        f"Must specify 'key' parameter to get attributes from a specific edge."
                    )
                if key not in edges_data:
                    raise PhyloZooValueError(f"Edge ({u}, {v}, {key}) does not exist in the network.")
                edge_attrs = edges_data[key]
                if attr is None:
                    return edge_attrs.copy()
                else:
                    return edge_attrs.get(attr)
        
        return {} if attr is None else None
    
    def get_branch_length(
        self,
        u: T,
        v: T,
        key: int | None = None
    ) -> float | None:
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
        """
        return self.get_edge_attribute(u, v, key, 'branch_length')
    
    def get_bootstrap(
        self,
        u: T,
        v: T,
        key: int | None = None
    ) -> float | None:
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
        """
        return self.get_edge_attribute(u, v, key, 'bootstrap')
    
    def get_gamma(
        self,
        u: T,
        v: T,
        key: int | None = None
    ) -> float | None:
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
        key : int | None, optional
            Edge key for parallel edges. Required if multiple parallel edges exist.
        
        Returns
        -------
        float | None
            Gamma value, or None if not set.
        
        Raises
        ------
        PhyloZooValueError
            If the edge is undirected (gamma can only be on directed edges).
        """
        # Gamma is only on directed hybrid edges - check if edge is undirected
        if self._graph._undirected.has_edge(u, v, key=key):
            raise PhyloZooValueError(
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
        key: int | None = None,
        directed: bool | None = None
    ) -> bool:
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
        directed : bool | None, optional
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
    ) -> Iterator[tuple[T, T] | tuple[T, T, int] | tuple[T, T, int, dict[str, Any]]]:
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
    ) -> Iterator[tuple[T, T] | tuple[T, T, int] | tuple[T, T, int, dict[str, Any]]]:
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
    ) -> Iterator[tuple[T, T] | tuple[T, T, int] | tuple[T, T, int, dict[str, Any]]]:
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
    def nodes(self):
        """
        View of all node IDs in the network.

        Cached. Same callable/view interface as the underlying graph's ``nodes``
        (e.g. ``net.nodes()``, ``net.nodes(data=True)``).

        Returns
        -------
        NodeView
            Set-like, callable view of node identifiers.
        """
        return self._graph.nodes

    @cached_property
    def edges(self):
        """
        View of all edges (directed and undirected) in the network.

        Cached. Same callable/view interface as the underlying graph's ``edges``
        (e.g. ``net.edges()``, ``net.edges(keys=True)``, ``net.edges(keys=True, data=True)``).

        Returns
        -------
        EdgeView
            Callable view of edges (u, v) or (u, v, key) or with data.
        """
        return self._graph.edges

    @cached_property
    def leaves(self) -> set[T]:
        """
        Get the set of leaf node IDs (nodes with degree 1, or degree 0 for single-node networks).
        
        Returns
        -------
        set[T]
            Set of leaf node identifiers. Returns a new set (which is mutable).
        """
        # Single-node network: the node is a leaf even though it has degree 0
        if self.number_of_nodes() == 1:
            return set(self._graph.nodes)
        return {node for node in self._graph.nodes if self._graph.degree(node) == 1}
    
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
        >>> net = MixedPhyNetwork(undirected_edges=[(3, 1), (3, 2)], nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})])
        >>> net.taxa
        {'A', 'B'}
        """
        return {self._node_to_label[leaf] for leaf in self.leaves}
    
    @cached_property
    def hybrid_nodes(self) -> set[T]:
        """
        Get the set of all hybrid nodes.
        
        A hybrid node is a node with in-degree >= 2 and total degree = in-degree + 1.
        
        Returns
        -------
        set[T]
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
    def hybrid_edges(self) -> set[tuple[T, T, int]]:
        """
        Get the set of all hybrid edges with keys.
        
        Hybrid edges are all directed edges.
        
        Returns
        -------
        set[tuple[T, T, int]]
            Set of (source, target, key) tuples for hybrid edges. Returns a new set (which is mutable).
        
        Examples
        --------
        >>> net = MixedPhyNetwork(
        ...     directed_edges=[(3, 2), (4, 2)],
        ...     undirected_edges=[(2, 1)],
        ...     nodes=[(1, {'label': 'A'})]
        ... )
        >>> net.hybrid_edges
        {(3, 2, 0), (4, 2, 0)}
        """
        return set(self._graph._directed.edges(keys=True))
    
    @cached_property
    def tree_nodes(self) -> set[T]:
        """
        Get the set of all tree nodes.
        
        A tree node is a node with in-degree 0 and total degree >= 3.
        
        Returns
        -------
        set[T]
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
    def internal_nodes(self) -> set[T]:
        """
        Get the set of all internal nodes.
        
        Internal nodes are all nodes that are not leaves.
        
        Returns
        -------
        set[T]
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
    def tree_edges(self) -> set[tuple[T, T, int]]:
        """
        Get the set of all tree edges with keys.
        
        Tree edges are simply all undirected edges.
        
        Returns
        -------
        set[tuple[T, T, int]]
            Set of (source, target, key) tuples for tree edges. Returns a new set (which is mutable).
        
        Examples
        --------
        >>> net = MixedPhyNetwork(
        ...     directed_edges=[(5, 4), (6, 4)],
        ...     undirected_edges=[(4, 1)],
        ...     nodes=[(1, {'label': 'A'})]
        ... )
        >>> net.tree_edges
        {(4, 1, 0)}
        """
        return set(self._graph._undirected.edges(keys=True))

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
