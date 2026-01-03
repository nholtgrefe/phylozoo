"""
Network I/O module.

This module provides functions for reading and writing directed phylogenetic networks
to/from various file formats (e.g., Newick, Nexus, eNewick, etc.).
"""

from __future__ import annotations

from typing import Any, TypeVar

from .base import DirectedPhyNetwork
from ...primitives.d_multigraph import DirectedMultiGraph

T = TypeVar('T')


def to_enewick(network: DirectedPhyNetwork) -> str:
    """
    Convert a DirectedPhyNetwork to Extended Newick (eNewick) format.
    
    This function serializes a directed phylogenetic network to the Extended Newick
    format, which supports hybrid nodes (reticulations) using the #H marker notation.
    
    Features:
    - Branch lengths are encoded as :length (e.g., :0.5)
    - Gamma and bootstrap values are encoded as comments: [&gamma=0.6,bootstrap=0.95]
    - Internal node labels are included when present
    - Hybrid nodes use Extended Newick #Hn markers
    - Output is deterministic (children sorted by node ID then label)
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to serialize.
    
    Returns
    -------
    str
        The eNewick string representation of the network.
    
    Raises
    ------
    ValueError
        If the network is empty (has no nodes).
    
    Examples
    --------
    >>> # Simple tree
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> to_enewick(net)
    '(A,B);'
    
    >>> # Tree with branch lengths
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         {'u': 3, 'v': 1, 'branch_length': 0.5},
    ...         {'u': 3, 'v': 2, 'branch_length': 0.3}
    ...     ],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> to_enewick(net)
    '(A:0.5,B:0.3);'
    
    >>> # Network with hybrid node
    >>> net = DirectedPhyNetwork(
    ...     edges=[
    ...         (5, 3), (5, 4),
    ...         {'u': 3, 'v': 2, 'gamma': 0.6},
    ...         {'u': 4, 'v': 2, 'gamma': 0.4},
    ...         (2, 1)
    ...     ],
    ...     nodes=[(1, {'label': 'A'})]
    ... )
    >>> to_enewick(net)  # doctest: +SKIP
    '((A)#H1,#H1);'
    
    Notes
    -----
    - For parallel edges between the same nodes, only the first edge is used
    - Labels containing special characters (spaces, parentheses, colons, etc.) 
      are automatically quoted with single quotes
    - The output is deterministic: multiple calls on the same network produce 
      the same string
    """
    if network.number_of_nodes() == 0:
        raise ValueError("Cannot convert empty network to eNewick")
    
    # Single node network
    if network.number_of_nodes() == 1:
        node = next(iter(network._graph.nodes))
        label = network.get_label(node)
        if label:
            return f"{_quote_label_if_needed(label)};"
        else:
            return f"{node};"
    
    # Initialize hybrid tracking
    hybrid_nodes_set = set(network.hybrid_nodes)
    hybrid_to_id: dict[T, int] = {}
    hybrid_counter = 1
    
    # Assign IDs to hybrid nodes deterministically
    for hybrid in sorted(hybrid_nodes_set, key=lambda n: (n, network.get_label(n) or '')):
        hybrid_to_id[hybrid] = hybrid_counter
        hybrid_counter += 1
    
    # Track which hybrids have been defined
    defined_hybrids: set[T] = set()
    
    def build_subtree(node: T, parent_edge_data: dict[str, Any] | None = None) -> str:
        """
        Recursively build eNewick string for subtree rooted at node.
        
        Parameters
        ----------
        node : T
            Current node to process.
        parent_edge_data : dict[str, Any] | None
            Edge data from parent to this node (for branch length, gamma, bootstrap).
        
        Returns
        -------
        str
            eNewick substring for this subtree.
        """
        # Check if this is a hybrid node
        is_hybrid = node in hybrid_nodes_set
        
        # If hybrid and already defined, return reference only
        if is_hybrid and node in defined_hybrids:
            # Just return the reference #Hn with parent edge attributes
            hybrid_id = hybrid_to_id[node]
            result = f"#H{hybrid_id}"
            result += _format_edge_attributes(parent_edge_data)
            return result
        
        # Mark hybrid as defined if this is first occurrence
        if is_hybrid:
            defined_hybrids.add(node)
        
        # Get children (sorted deterministically)
        children_list = list(network.children(node))
        children_sorted = sorted(children_list, key=lambda n: (n, network.get_label(n) or ''))
        
        # If this is a leaf, just return the label
        if len(children_sorted) == 0:
            label = network.get_label(node)
            if label:
                result = _quote_label_if_needed(label)
            else:
                result = str(node)
            
            # Add hybrid marker if this is a hybrid leaf
            if is_hybrid:
                hybrid_id = hybrid_to_id[node]
                result += f"#H{hybrid_id}"
            
            # Add parent edge attributes
            result += _format_edge_attributes(parent_edge_data)
            return result
        
        # Internal node: process children
        child_strings = []
        for child in children_sorted:
            # Get edge data from node to child
            # Handle parallel edges by using the first edge
            edge_data = _get_first_edge_data(network, node, child)
            child_str = build_subtree(child, edge_data)
            child_strings.append(child_str)
        
        # Build the result string
        result = f"({','.join(child_strings)})"
        
        # Add internal node label if present
        label = network.get_label(node)
        if label:
            result += _quote_label_if_needed(label)
        
        # Add hybrid marker if this is a hybrid
        if is_hybrid:
            hybrid_id = hybrid_to_id[node]
            result += f"#H{hybrid_id}"
        
        # Add parent edge attributes
        result += _format_edge_attributes(parent_edge_data)
        
        return result
    
    # Start from root
    root = network.root_node
    enewick_str = build_subtree(root, None)
    
    return enewick_str + ";"


def _get_first_edge_data(network: DirectedPhyNetwork, u: T, v: T) -> dict[str, Any]:
    """
    Get edge data for the first edge from u to v.
    
    For parallel edges, returns data from the first edge encountered.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The network.
    u, v : T
        Edge endpoints.
    
    Returns
    -------
    dict[str, Any]
        Edge data dictionary.
    """
    # Use incident_child_edges to get all edges from u
    for edge_tuple in network.incident_child_edges(u, keys=True, data=True):
        if len(edge_tuple) == 4:
            edge_u, edge_v, key, data = edge_tuple
            if edge_v == v:
                return dict(data) if data else {}
    
    return {}


def _format_edge_attributes(edge_data: dict[str, Any] | None) -> str:
    """
    Format edge attributes (branch_length, gamma, bootstrap) for eNewick.
    
    Parameters
    ----------
    edge_data : dict[str, Any] | None
        Edge data dictionary.
    
    Returns
    -------
    str
        Formatted attribute string (e.g., "[&gamma=0.6,bootstrap=0.95]:0.5").
    """
    if edge_data is None:
        return ""
    
    result = ""
    
    # Build comment section for gamma and bootstrap
    comment_parts = []
    if 'gamma' in edge_data:
        gamma_val = edge_data['gamma']
        comment_parts.append(f"gamma={gamma_val}")
    if 'bootstrap' in edge_data:
        bootstrap_val = edge_data['bootstrap']
        comment_parts.append(f"bootstrap={bootstrap_val}")
    
    if comment_parts:
        result += f"[&{','.join(comment_parts)}]"
    
    # Add branch length
    if 'branch_length' in edge_data:
        branch_length = edge_data['branch_length']
        result += f":{branch_length}"
    
    return result


def _quote_label_if_needed(label: str) -> str:
    """
    Quote a label if it contains special characters.
    
    Special characters that require quoting: spaces, parentheses, colons,
    semicolons, commas, brackets, quotes.
    
    Parameters
    ----------
    label : str
        The label to potentially quote.
    
    Returns
    -------
    str
        The label, quoted if necessary.
    """
    special_chars = {' ', '(', ')', ':', ';', ',', '[', ']', "'", '"', '#'}
    
    if any(char in label for char in special_chars):
        # Escape single quotes in the label
        escaped_label = label.replace("'", "''")
        return f"'{escaped_label}'"
    
    return label
