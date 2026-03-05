"""
Network I/O module.

This module registers format handlers for DirectedPhyNetwork with FormatRegistry
for use with the IOMixin system.

All eNewick functionality (parser, writer, reader) is in the _enewick module.
DOT format support uses the underlying DirectedMultiGraph functions, but uses
labels as node names where possible.
"""

from __future__ import annotations

import re
from typing import Any

from phylozoo.utils.io import FormatRegistry
from .base import DirectedPhyNetwork
from ._enewick import to_enewick, from_enewick
from ...primitives.d_multigraph.io import (
    _escape_dot_string,
    _format_dot_attributes,
    _parse_dot_attributes,
    _convert_node_id,
    from_dot as dmgraph_from_dot,
)
from ....utils.exceptions import PhyloZooFormatError

# Register eNewick format handlers with FormatRegistry
FormatRegistry.register(
    DirectedPhyNetwork, 'enewick',
    reader=from_enewick,
    writer=to_enewick,
    extensions=['.enewick', '.eNewick', 'enwk', '.nwk', '.newick'],
    default=True
)


def to_dot(network: DirectedPhyNetwork, **kwargs: Any) -> str:
    """
    Convert a DirectedPhyNetwork to a DOT format string.
    
    This function uses labels as node names in DOT where available, falling back
    to node IDs if no label is present. This makes the DOT output more readable
    for phylogenetic networks where leaves have meaningful taxon labels.
    
    Parameters
    ----------
    network : DirectedPhyNetwork
        The directed phylogenetic network to convert.
    **kwargs
        Additional arguments:
        - graph_name (str): Optional name for the graph (default: '').
    
    Returns
    -------
    str
        The DOT format string representation of the network.
    
    Examples
    --------
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> dot_str = to_dot(net)
    >>> 'digraph' in dot_str
    True
    >>> 'A ->' in dot_str or '-> A' in dot_str or 'A [' in dot_str
    True
    
    Notes
    -----

    - Nodes with labels use the label as the node name in DOT
    - Nodes without labels use the node ID as the node name
    - All node and edge attributes are preserved
    - Graph attributes are included
    - Parallel edges are supported
    """
    lines = []
    
    # Graph name (optional)
    graph_name = kwargs.get('graph_name', '')
    if graph_name:
        lines.append(f'digraph {_escape_dot_string(graph_name)} {{')
    else:
        lines.append('digraph {')
    
    # Graph attributes
    graph_attrs = network.get_network_attribute()
    if graph_attrs:
        for key, value in graph_attrs.items():
            if isinstance(value, str):
                value_str = _escape_dot_string(value)
            else:
                value_str = str(value)
            lines.append(f'    {key}={value_str};')
    
    # Node declarations with attributes
    # Use label as node name if available, otherwise use node ID
    for node in network._graph.nodes():
        node_attrs = {}
        if node in network._graph._graph:
            node_data = network._graph._graph.nodes[node]
            if node_data:
                node_attrs = dict(node_data)
        
        # Get label if available
        label = network.get_label(node)
        
        # Use label as node name in DOT if available, otherwise use node ID
        if label:
            node_name = label
            # Remove label from attributes if it matches the node name (to avoid duplication)
            if 'label' in node_attrs and node_attrs['label'] == label:
                # Create a copy without the label attribute
                node_attrs = {k: v for k, v in node_attrs.items() if k != 'label'}
        else:
            node_name = str(node)
            # Use node_id as label attribute if no label attribute exists
            if 'label' not in node_attrs:
                node_attrs['label'] = str(node)
        
        node_name_str = _escape_dot_string(node_name)
        attrs_str = _format_dot_attributes(node_attrs)
        
        if attrs_str:
            lines.append(f'    {node_name_str} {attrs_str};')
        else:
            lines.append(f'    {node_name_str};')
    
    # Edge declarations with attributes
    # Map node IDs to their DOT names (labels or IDs)
    node_to_dot_name: dict[Any, str] = {}
    for node in network._graph.nodes():
        label = network.get_label(node)
        if label:
            node_to_dot_name[node] = label
        else:
            node_to_dot_name[node] = str(node)
    
    for u, v, key, data in network._graph.edges(keys=True, data=True):
        u_name = node_to_dot_name[u]
        v_name = node_to_dot_name[v]
        
        u_str = _escape_dot_string(u_name)
        v_str = _escape_dot_string(v_name)
        
        # Include key in edge attributes if there are parallel edges
        edge_attrs = dict(data) if data else {}
        
        # Add key as attribute if there are multiple edges between u and v
        if network._graph._graph.number_of_edges(u, v) > 1:
            edge_attrs['key'] = key
        
        attrs_str = _format_dot_attributes(edge_attrs)
        
        if attrs_str:
            lines.append(f'    {u_str} -> {v_str} {attrs_str};')
        else:
            lines.append(f'    {u_str} -> {v_str};')
    
    lines.append('}')
    
    return '\n'.join(lines) + '\n'


def from_dot(dot_string: str, **kwargs: Any) -> DirectedPhyNetwork:
    """
    Parse a DOT format string and create a DirectedPhyNetwork.
    
    This function parses DOT format and creates a DirectedPhyNetwork. Node names
    in DOT are used as labels if they are valid labels (strings), otherwise they
    are used as node IDs.
    
    Parameters
    ----------
    dot_string : str
        DOT format string containing graph data.
    **kwargs
        Additional arguments (currently unused, for compatibility).
    
    Returns
    -------
    DirectedPhyNetwork
        Parsed directed phylogenetic network.
    
    Raises
    ------
    PhyloZooFormatError
        If the DOT string is malformed or cannot be parsed, or if the resulting
        network is invalid according to DirectedPhyNetwork validation rules.
    
    Examples
    --------
    >>> dot_str = '''digraph {
    ...     A [label="Species A"];
    ...     B [label="Species B"];
    ...     root -> A;
    ...     root -> B;
    ... }'''
    >>> net = from_dot(dot_str)
    >>> net.number_of_nodes()
    3
    >>> 'A' in net.taxa or net.get_label(list(net.leaves)[0]) == 'A'
    True
    
    Notes
    -----

    - Node names in DOT become labels if they are strings
    - If a node has a 'label' attribute in DOT, it overrides the node name
    - Edge attributes are preserved
    - Graph attributes are preserved
    - The network is validated after creation
    """
    # First parse as DirectedMultiGraph to get the structure
    dmgraph = dmgraph_from_dot(dot_string, **kwargs)
    
    # Extract edges with attributes
    edges: list[dict[str, Any]] = []
    for u, v, key, data in dmgraph.edges(keys=True, data=True):
        edge_dict: dict[str, Any] = {'u': u, 'v': v}
        if key != 0:
            edge_dict['key'] = key
        if data:
            edge_dict.update(data)
        edges.append(edge_dict)
    
    # Extract nodes with labels
    # In DOT, node names are used as labels if they look like labels (strings)
    # If a node has a 'label' attribute, use that instead
    nodes: list[tuple[Any, dict[str, Any]]] = []
    for node in dmgraph.nodes():
        node_attrs = {}
        if node in dmgraph._graph:
            node_data = dmgraph._graph.nodes[node]
            if node_data:
                node_attrs = dict(node_data)
        
        # If node name is a string and looks like a label, use it as label
        # unless there's already a 'label' attribute
        if isinstance(node, str) and 'label' not in node_attrs:
            # Use node name as label (it's already a string)
            node_attrs['label'] = node
        elif 'label' in node_attrs:
            # Label attribute already exists, ensure it's a string
            if not isinstance(node_attrs['label'], str):
                node_attrs['label'] = str(node_attrs['label'])
        # If node is not a string and has no label attribute, don't add a label
        # DirectedPhyNetwork will auto-label leaves if needed
        
        nodes.append((node, node_attrs))
    
    # Extract graph attributes
    graph_attrs = dmgraph._graph.graph.copy() if dmgraph._graph.graph else None
    
    # Create DirectedPhyNetwork
    return DirectedPhyNetwork(
        edges=edges,
        nodes=nodes if nodes else None,
        attributes=graph_attrs if graph_attrs else None
    )


# Register DOT format handlers with FormatRegistry
FormatRegistry.register(
    DirectedPhyNetwork, 'dot',
    reader=from_dot,
    writer=to_dot,
    extensions=['.dot', '.gv'],
)
