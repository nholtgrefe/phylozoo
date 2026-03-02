"""
Mixed multi-graph I/O module.

This module provides format handlers for reading and writing mixed multi-graphs
to/from files. Format handlers are registered with FormatRegistry for use with
the IOMixin system.

The following format handlers are defined and registered:

- **phylozoo-dot**: PhyloZoo DOT format (extensions: .pzdot).
  Writer: `to_phylozoo_dot()` converts MixedMultiGraph to phylozoo-dot string.
  Reader: `from_phylozoo_dot()` parses phylozoo-dot string to MixedMultiGraph.

These handlers are automatically registered when this module is imported.
MixedMultiGraph inherits from IOMixin, so you can use:
- `graph.save('file.pzdot')` - Save to file (auto-detects format)
- `graph.load('file.pzdot')` - Load from file (auto-detects format)
- `graph.to_string(format='phylozoo-dot')` - Convert to string
- `graph.from_string(string, format='phylozoo-dot')` - Parse from string
- `MixedMultiGraph.convert('in.pzdot', 'out.pzdot')` - Convert between formats

Notes
-----
PhyloZoo DOT format:
- Similar to Graphviz DOT format but supports both directed and undirected edges
- Uses `--` for undirected edges
- Uses `->` for directed edges
- Supports node attributes (label, shape, color, etc.)
- Supports edge attributes (label, weight, color, etc.)
- Supports graph attributes
- Supports parallel edges (multigraph)
- Uses node_id as the label if no explicit label is provided
"""

from __future__ import annotations

import re
from typing import Any

from phylozoo.utils.io import FormatRegistry
from ....utils.exceptions import PhyloZooParseError
from .base import MixedMultiGraph


def _escape_dot_string(s: str) -> str:
    """
    Escape a string for use in DOT format.
    
    Parameters
    ----------
    s : str
        String to escape.
    
    Returns
    -------
    str
        Escaped string.
    """
    # If string contains special characters or spaces, wrap in quotes
    if any(c in s for c in [' ', '\t', '\n', '"', '\\', '[', ']', '{', '}', '-', '>']):
        # Escape backslashes and quotes
        s = s.replace('\\', '\\\\')
        s = s.replace('"', '\\"')
        return f'"{s}"'
    return s


def _format_dot_attributes(attrs: dict[str, Any]) -> str:
    """
    Format attributes for DOT format.
    
    Parameters
    ----------
    attrs : dict[str, Any]
        Dictionary of attributes.
    
    Returns
    -------
    str
        Formatted attribute string like '[key1=value1, key2=value2]'.
    """
    if not attrs:
        return ''
    
    parts = []
    for key, value in attrs.items():
        # Convert value to string and escape if needed
        if isinstance(value, str):
            value_str = _escape_dot_string(value)
        elif isinstance(value, (int, float)):
            value_str = str(value)
        elif isinstance(value, bool):
            value_str = 'true' if value else 'false'
        else:
            value_str = _escape_dot_string(str(value))
        
        parts.append(f'{key}={value_str}')
    
    return '[' + ', '.join(parts) + ']'


def to_phylozoo_dot(graph: MixedMultiGraph, **kwargs: Any) -> str:
    """
    Convert a MixedMultiGraph to a phylozoo-dot format string.
    
    Parameters
    ----------
    graph : MixedMultiGraph
        The mixed multi-graph to convert.
    **kwargs
        Additional arguments (currently unused, for compatibility).
    
    Returns
    -------
    str
        The phylozoo-dot format string representation of the graph.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
    >>> from phylozoo.core.primitives.m_multigraph.io import to_phylozoo_dot
    >>> 
    >>> G = MixedMultiGraph()
    >>> G.add_directed_edge(1, 2, weight=1.0)
    0
    >>> G.add_undirected_edge(2, 3, weight=2.0)
    0
    >>> pzdot_str = to_phylozoo_dot(G)
    >>> 'graph' in pzdot_str
    True
    >>> '1 -> 2' in pzdot_str
    True
    >>> '2 -- 3' in pzdot_str
    True
    
    Notes
    -----
    The phylozoo-dot format includes:
    - graph declaration (not digraph, since we have both types)
    - Node declarations with attributes
    - Undirected edge declarations with `--`
    - Directed edge declarations with `->`
    - Graph attributes (if any)
    - Support for parallel edges (multigraph)
    """
    lines = []
    
    # Graph name (optional, use empty string)
    graph_name = kwargs.get('graph_name', '')
    if graph_name:
        lines.append(f'graph {_escape_dot_string(graph_name)} {{')
    else:
        lines.append('graph {')
    
    # Graph attributes (from either subgraph, they should be the same)
    if hasattr(graph, '_directed') and hasattr(graph._directed, 'graph'):
        graph_attrs = graph._directed.graph
        if graph_attrs:
            for key, value in graph_attrs.items():
                if isinstance(value, str):
                    value_str = _escape_dot_string(value)
                else:
                    value_str = str(value)
                lines.append(f'    {key}={value_str};')
    
    # Get all nodes (from combined graph)
    all_nodes = set(graph.nodes())
    
    # Node declarations with attributes
    for node in all_nodes:
        node_attrs = {}
        # Check both directed and undirected graphs for node attributes
        if hasattr(graph, '_directed') and node in graph._directed:
            node_data = graph._directed.nodes[node]
            if node_data:
                node_attrs.update(node_data)
        if hasattr(graph, '_undirected') and node in graph._undirected:
            node_data = graph._undirected.nodes[node]
            if node_data:
                node_attrs.update(node_data)
        
        # Use node_id as label if no label attribute
        if 'label' not in node_attrs:
            node_attrs['label'] = str(node)
        
        node_id_str = _escape_dot_string(str(node))
        attrs_str = _format_dot_attributes(node_attrs)
        
        if attrs_str:
            lines.append(f'    {node_id_str} {attrs_str};')
        else:
            lines.append(f'    {node_id_str};')
    
    # Undirected edge declarations with attributes
    for u, v, key, data in graph.undirected_edges_iter(keys=True, data=True):
        u_str = _escape_dot_string(str(u))
        v_str = _escape_dot_string(str(v))
        
        # Include key in edge attributes if there are parallel edges
        edge_attrs = dict(data) if data else {}
        
        # Add key as attribute if there are multiple edges between u and v
        if graph._undirected.number_of_edges(u, v) > 1:
            edge_attrs['key'] = key
        
        attrs_str = _format_dot_attributes(edge_attrs)
        
        if attrs_str:
            lines.append(f'    {u_str} -- {v_str} {attrs_str};')
        else:
            lines.append(f'    {u_str} -- {v_str};')
    
    # Directed edge declarations with attributes
    for u, v, key, data in graph.directed_edges_iter(keys=True, data=True):
        u_str = _escape_dot_string(str(u))
        v_str = _escape_dot_string(str(v))
        
        # Include key in edge attributes if there are parallel edges
        edge_attrs = dict(data) if data else {}
        
        # Add key as attribute if there are multiple edges between u and v
        if graph._directed.number_of_edges(u, v) > 1:
            edge_attrs['key'] = key
        
        attrs_str = _format_dot_attributes(edge_attrs)
        
        if attrs_str:
            lines.append(f'    {u_str} -> {v_str} {attrs_str};')
        else:
            lines.append(f'    {u_str} -> {v_str};')
    
    lines.append('}')
    
    return '\n'.join(lines) + '\n'


def from_phylozoo_dot(pzdot_string: str, **kwargs: Any) -> MixedMultiGraph:
    """
    Parse a phylozoo-dot format string and create a MixedMultiGraph.
    
    Parameters
    ----------
    pzdot_string : str
        PhyloZoo DOT format string containing graph data.
    **kwargs
        Additional arguments (currently unused, for compatibility).
    
    Returns
    -------
    MixedMultiGraph
        Parsed mixed multi-graph.
    
    Raises
    ------
    PhyloZooParseError
        If the phylozoo-dot string is malformed or cannot be parsed.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.io import from_phylozoo_dot
    >>> 
    >>> pzdot_str = '''graph {
    ...     1 [label="Node1"];
    ...     2 [label="Node2"];
    ...     1 -> 2 [weight=1.0];
    ...     2 -- 3 [weight=2.0];
    ... }'''
    >>> 
    >>> G = from_phylozoo_dot(pzdot_str)
    >>> G.number_of_nodes()
    3
    >>> G.number_of_edges()
    2
    
    Notes
    -----
    This parser expects:
    - graph declaration (not digraph)
    - Node declarations (optional attributes)
    - Undirected edge declarations with `--` (optional attributes)
    - Directed edge declarations with `->` (optional attributes)
    - Graph attributes (optional)
    - Support for parallel edges
    """
    # Remove comments
    lines = []
    for line in pzdot_string.split('\n'):
        # Remove C-style comments (// and /* */)
        if '//' in line:
            line = line[:line.index('//')]
        lines.append(line)
    
    content = '\n'.join(lines)
    
    # Extract graph name and body
    graph_match = re.search(r'graph\s+(\w+|"[^"]+")?\s*\{', content, re.IGNORECASE)
    if not graph_match:
        raise PhyloZooParseError("Could not find graph declaration in phylozoo-dot string")
    
    # Extract graph body (between { and })
    brace_count = 0
    start_idx = content.index('{')
    end_idx = start_idx
    
    for i, char in enumerate(content[start_idx:], start=start_idx):
        if char == '{':
            brace_count += 1
        elif char == '}':
            brace_count -= 1
            if brace_count == 0:
                end_idx = i
                break
    
    if brace_count != 0:
        raise PhyloZooParseError("Unmatched braces in phylozoo-dot string")
    
    body = content[start_idx + 1:end_idx]
    
    # Parse graph attributes, nodes, and edges
    graph_attrs = {}
    nodes_data: dict[Any, dict[str, Any]] = {}
    undirected_edges_data: list[tuple[Any, Any, int | None, dict[str, Any]]] = []
    directed_edges_data: list[tuple[Any, Any, int | None, dict[str, Any]]] = []
    
    # Node pattern: node_id [attributes];
    # Matches: word characters, quoted strings, or numeric values (including floats)
    node_pattern = r'([\w.]+|"[^"]+")\s*(?:\[([^\]]+)\])?\s*;'
    
    # Undirected edge pattern: u -- v [attributes];
    # Matches: word characters, quoted strings, or numeric values (including floats)
    undirected_edge_pattern = r'([\w.]+|"[^"]+")\s*--\s*([\w.]+|"[^"]+")\s*(?:\[([^\]]+)\])?\s*;'
    
    # Directed edge pattern: u -> v [attributes];
    # Matches: word characters, quoted strings, or numeric values (including floats)
    directed_edge_pattern = r'([\w.]+|"[^"]+")\s*->\s*([\w.]+|"[^"]+")\s*(?:\[([^\]]+)\])?\s*;'
    
    # Graph attribute pattern: key=value; (standalone, not in brackets)
    graph_attr_pattern = r'^(\w+)\s*=\s*([^;]+);$'
    
    for line in body.split('\n'):
        line = line.strip()
        if not line or line.startswith('//') or line.startswith('#'):
            continue
        
        # Try to match graph attribute first (standalone key=value;)
        graph_attr_match = re.match(graph_attr_pattern, line)
        if graph_attr_match:
            key = graph_attr_match.group(1).strip()
            value = graph_attr_match.group(2).strip().strip('"\'')
            graph_attrs[key] = value
            continue
        
        # Try to match undirected edge (contains --)
        undirected_edge_match = re.search(undirected_edge_pattern, line)
        if undirected_edge_match:
            u_str = undirected_edge_match.group(1).strip('"\'')
            v_str = undirected_edge_match.group(2).strip('"\'')
            attrs_str = undirected_edge_match.group(3) if undirected_edge_match.group(3) else ''
            
            # Parse attributes
            edge_attrs = _parse_dot_attributes(attrs_str)
            
            # Extract key if present
            key = None
            if 'key' in edge_attrs:
                try:
                    key = int(edge_attrs.pop('key'))
                except (ValueError, TypeError):
                    pass
            
            # Convert node strings to appropriate types
            u = _convert_node_id(u_str)
            v = _convert_node_id(v_str)
            
            undirected_edges_data.append((u, v, key, edge_attrs))
            continue
        
        # Try to match directed edge (contains ->)
        directed_edge_match = re.search(directed_edge_pattern, line)
        if directed_edge_match:
            u_str = directed_edge_match.group(1).strip('"\'')
            v_str = directed_edge_match.group(2).strip('"\'')
            attrs_str = directed_edge_match.group(3) if directed_edge_match.group(3) else ''
            
            # Parse attributes
            edge_attrs = _parse_dot_attributes(attrs_str)
            
            # Extract key if present
            key = None
            if 'key' in edge_attrs:
                try:
                    key = int(edge_attrs.pop('key'))
                except (ValueError, TypeError):
                    pass
            
            # Convert node strings to appropriate types
            u = _convert_node_id(u_str)
            v = _convert_node_id(v_str)
            
            directed_edges_data.append((u, v, key, edge_attrs))
            continue
        
        # Try to match node
        node_match = re.search(node_pattern, line)
        if node_match:
            node_str = node_match.group(1).strip('"\'')
            attrs_str = node_match.group(2) if node_match.group(2) else ''
            
            # Parse attributes
            node_attrs = _parse_dot_attributes(attrs_str)
            
            # Convert node string to appropriate type
            node_id = _convert_node_id(node_str)
            
            nodes_data[node_id] = node_attrs
            continue
    
    # Create graph
    graph = MixedMultiGraph(
        attributes=graph_attrs if graph_attrs else None,
        directed_edges=None,
        undirected_edges=None
    )
    
    # Add nodes with attributes
    for node_id, attrs in nodes_data.items():
        graph.add_node(node_id, **attrs)
    
    # Add undirected edges with attributes
    for u, v, key, attrs in undirected_edges_data:
        # Ensure nodes exist
        if u not in graph:
            graph.add_node(u)
        if v not in graph:
            graph.add_node(v)
        
        graph.add_undirected_edge(u, v, key=key, **attrs)
    
    # Add directed edges with attributes
    for u, v, key, attrs in directed_edges_data:
        # Ensure nodes exist
        if u not in graph:
            graph.add_node(u)
        if v not in graph:
            graph.add_node(v)
        
        graph.add_directed_edge(u, v, key=key, **attrs)
    
    return graph


def _parse_dot_attributes(attrs_str: str) -> dict[str, Any]:
    """
    Parse DOT attribute string like 'key1=value1, key2=value2'.
    
    Parameters
    ----------
    attrs_str : str
        Attribute string.
    
    Returns
    -------
    dict[str, Any]
        Dictionary of attributes.
    """
    attrs = {}
    if not attrs_str.strip():
        return attrs
    
    # Split by comma, but respect quoted strings
    parts = []
    current = ''
    in_quotes = False
    escape_next = False
    
    for char in attrs_str:
        if escape_next:
            current += char
            escape_next = False
            continue
        
        if char == '\\':
            escape_next = True
            current += char
            continue
        
        if char == '"' or char == "'":
            in_quotes = not in_quotes
            current += char
            continue
        
        if char == ',' and not in_quotes:
            parts.append(current.strip())
            current = ''
        else:
            current += char
    
    if current.strip():
        parts.append(current.strip())
    
    # Parse each key=value pair
    for part in parts:
        if '=' not in part:
            continue
        
        key, value = part.split('=', 1)
        key = key.strip()
        value = value.strip().strip('"\'')
        
        # Try to convert to appropriate type
        if value.lower() == 'true':
            attrs[key] = True
        elif value.lower() == 'false':
            attrs[key] = False
        else:
            # Try numeric conversion
            try:
                if '.' in value:
                    attrs[key] = float(value)
                else:
                    attrs[key] = int(value)
            except ValueError:
                attrs[key] = value
    
    return attrs


def _convert_node_id(node_str: str) -> Any:
    """
    Convert node string to appropriate type (int, float, or str).
    
    Parameters
    ----------
    node_str : str
        Node string.
    
    Returns
    -------
    Any
        Converted node ID.
    """
    # Try int first
    try:
        return int(node_str)
    except ValueError:
        pass
    
    # Try float
    try:
        return float(node_str)
    except ValueError:
        pass
    
    # Keep as string
    return node_str


# Register format handlers with FormatRegistry
FormatRegistry.register(
    MixedMultiGraph, 'phylozoo-dot',
    reader=from_phylozoo_dot,
    writer=to_phylozoo_dot,
    extensions=['.pzdot'],
    default=True
)

