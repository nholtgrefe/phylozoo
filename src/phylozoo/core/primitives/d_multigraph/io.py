"""
Directed multi-graph I/O module.

This module provides format handlers for reading and writing directed multi-graphs
to/from files. Format handlers are registered with FormatRegistry for use with
the IOMixin system.

The following format handlers are defined and registered:
- **dot**: DOT format (Graphviz) (extensions: .dot, .gv)
  - Writer: `to_dot()` - Converts DirectedMultiGraph to DOT string
  - Reader: `from_dot()` - Parses DOT string to DirectedMultiGraph
- **edgelist**: Edge-list format (extensions: .el)
  - Writer: `to_edgelist()` - Converts DirectedMultiGraph to edge-list string
  - Reader: `from_edgelist()` - Parses edge-list string to DirectedMultiGraph

These handlers are automatically registered when this module is imported.
DirectedMultiGraph inherits from IOMixin, so you can use:
- `graph.save('file.dot')` - Save to file (auto-detects format)
- `graph.load('file.dot')` - Load from file (auto-detects format)
- `graph.to_string(format='dot')` - Convert to string
- `graph.from_string(string, format='edgelist')` - Parse from string
- `DirectedMultiGraph.convert('in.dot', 'out.el')` - Convert between formats

Notes
-----
DOT format supports:
- Node attributes (label, shape, color, etc.)
- Edge attributes (label, weight, color, etc.)
- Graph attributes
- Parallel edges (multigraph support)

Edge-list format:
- Simple text format: one edge per line
- Format: `u v` or `u v key` or `u v key attr1=value1 attr2=value2`
- Uses node_id as the label/name
"""

from __future__ import annotations

import re
from typing import Any

from phylozoo.utils.exceptions import PhyloZooParseError

from phylozoo.utils.io import FormatRegistry
from .base import DirectedMultiGraph


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


def to_dot(graph: DirectedMultiGraph, **kwargs: Any) -> str:
    """
    Convert a DirectedMultiGraph to a DOT format string.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The directed multi-graph to convert.
    **kwargs
        Additional arguments (currently unused, for compatibility).
    
    Returns
    -------
    str
        The DOT format string representation of the graph.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
    >>> from phylozoo.core.primitives.d_multigraph.io import to_dot
    >>> 
    >>> G = DirectedMultiGraph()
    >>> G.add_edge(1, 2, weight=1.0)
    0
    >>> G.add_edge(2, 3, weight=2.0)
    0
    >>> dot_str = to_dot(G)
    >>> 'digraph' in dot_str
    True
    >>> '1 -> 2' in dot_str
    True
    
    Notes
    -----
    The DOT format includes:
    - digraph declaration
    - Node declarations with attributes
    - Edge declarations with attributes
    - Graph attributes (if any)
    - Support for parallel edges (multigraph)
    """
    lines = []
    
    # Graph name (optional, use empty string)
    graph_name = kwargs.get('graph_name', '')
    if graph_name:
        lines.append(f'digraph {_escape_dot_string(graph_name)} {{')
    else:
        lines.append('digraph {')
    
    # Graph attributes
    if hasattr(graph, '_graph') and hasattr(graph._graph, 'graph'):
        graph_attrs = graph._graph.graph
        if graph_attrs:
            for key, value in graph_attrs.items():
                if isinstance(value, str):
                    value_str = _escape_dot_string(value)
                else:
                    value_str = str(value)
                lines.append(f'    {key}={value_str};')
    
    # Node declarations with attributes
    for node in graph.nodes():
        node_attrs = {}
        if hasattr(graph, '_graph') and node in graph._graph:
            node_data = graph._graph.nodes[node]
            if node_data:
                node_attrs = dict(node_data)
        
        # Use node_id as label if no label attribute (as per user requirement)
        if 'label' not in node_attrs:
            node_attrs['label'] = str(node)
        
        node_id_str = _escape_dot_string(str(node))
        attrs_str = _format_dot_attributes(node_attrs)
        
        if attrs_str:
            lines.append(f'    {node_id_str} {attrs_str};')
        else:
            lines.append(f'    {node_id_str};')
    
    # Edge declarations with attributes
    for u, v, key, data in graph.edges_iter(keys=True, data=True):
        u_str = _escape_dot_string(str(u))
        v_str = _escape_dot_string(str(v))
        
        # Include key in edge attributes if there are parallel edges
        edge_attrs = dict(data) if data else {}
        
        # Add key as attribute if there are multiple edges between u and v
        if graph._graph.number_of_edges(u, v) > 1:
            edge_attrs['key'] = key
        
        attrs_str = _format_dot_attributes(edge_attrs)
        
        if attrs_str:
            lines.append(f'    {u_str} -> {v_str} {attrs_str};')
        else:
            lines.append(f'    {u_str} -> {v_str};')
    
    lines.append('}')
    
    return '\n'.join(lines) + '\n'


def from_dot(dot_string: str, **kwargs: Any) -> DirectedMultiGraph:
    """
    Parse a DOT format string and create a DirectedMultiGraph.
    
    Parameters
    ----------
    dot_string : str
        DOT format string containing graph data.
    **kwargs
        Additional arguments (currently unused, for compatibility).
    
    Returns
    -------
    DirectedMultiGraph
        Parsed directed multi-graph.
    
    Raises
    ------
    PhyloZooParseError
        If the DOT string is malformed or cannot be parsed.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph.io import from_dot
    >>> 
    >>> dot_str = '''digraph {
    ...     1 [label="Node1"];
    ...     2 [label="Node2"];
    ...     1 -> 2 [weight=1.0];
    ...     2 -> 3 [weight=2.0];
    ... }'''
    >>> 
    >>> G = from_dot(dot_str)
    >>> G.number_of_nodes()
    3
    >>> G.number_of_edges()
    2
    
    Notes
    -----
    This parser expects:
    - digraph declaration
    - Node declarations (optional attributes)
    - Edge declarations (optional attributes)
    - Graph attributes (optional)
    - Support for parallel edges
    """
    # Remove comments
    lines = []
    for line in dot_string.split('\n'):
        # Remove C-style comments (// and /* */)
        # Remove # comments
        if '//' in line:
            line = line[:line.index('//')]
        if '#' in line and not line.strip().startswith('#'):
            # Only remove # if it's not part of a string
            pass  # Keep for now, will handle in parsing
        lines.append(line)
    
    content = '\n'.join(lines)
    
    # Extract graph name and body
    digraph_match = re.search(r'digraph\s+(\w+|"[^"]+")?\s*\{', content, re.IGNORECASE)
    if not digraph_match:
        raise PhyloZooParseError("Could not find digraph declaration in DOT string")
    
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
        raise PhyloZooParseError("Unmatched braces in DOT string")
    
    body = content[start_idx + 1:end_idx]
    
    # Parse graph attributes, nodes, and edges
    graph_attrs = {}
    nodes_data: dict[Any, dict[str, Any]] = {}
    edges_data: list[tuple[Any, Any, int | None, dict[str, Any]]] = []
    
    # Node pattern: node_id [attributes];
    node_pattern = r'(\w+|"[^"]+")\s*(?:\[([^\]]+)\])?\s*;'
    
    # Edge pattern: u -> v [attributes];
    edge_pattern = r'(\w+|"[^"]+")\s*->\s*(\w+|"[^"]+")\s*(?:\[([^\]]+)\])?\s*;'
    
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
        
        # Try to match edge (edges contain ->)
        edge_match = re.search(edge_pattern, line)
        if edge_match:
            u_str = edge_match.group(1).strip('"\'')
            v_str = edge_match.group(2).strip('"\'')
            attrs_str = edge_match.group(3) if edge_match.group(3) else ''
            
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
            
            edges_data.append((u, v, key, edge_attrs))
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
    graph = DirectedMultiGraph(attributes=graph_attrs if graph_attrs else None)
    
    # Add nodes with attributes
    for node_id, attrs in nodes_data.items():
        graph.add_node(node_id, **attrs)
    
    # Add edges with attributes
    for u, v, key, attrs in edges_data:
        # Ensure nodes exist
        if u not in graph:
            graph.add_node(u)
        if v not in graph:
            graph.add_node(v)
        
        graph.add_edge(u, v, key=key, **attrs)
    
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


def to_edgelist(graph: DirectedMultiGraph, **kwargs: Any) -> str:
    """
    Convert a DirectedMultiGraph to an edge-list format string.
    
    Parameters
    ----------
    graph : DirectedMultiGraph
        The directed multi-graph to convert.
    **kwargs
        Additional arguments (currently unused, for compatibility).
    
    Returns
    -------
    str
        The edge-list format string representation of the graph.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
    >>> from phylozoo.core.primitives.d_multigraph.io import to_edgelist
    >>> 
    >>> G = DirectedMultiGraph()
    >>> G.add_edge(1, 2, weight=1.0)
    0
    >>> G.add_edge(2, 3, weight=2.0)
    0
    >>> el_str = to_edgelist(G)
    >>> '1 2' in el_str
    True
    >>> '2 3' in el_str
    True
    
    Notes
    -----
    The edge-list format:
    - One edge per line
    - Format: `u v` or `u v key` or `u v key attr1=value1 attr2=value2`
    - Uses node_id as the label/name
    - Includes edge keys for parallel edges
    - Includes edge attributes if present
    """
    lines = []
    
    for u, v, key, data in graph.edges_iter(keys=True, data=True):
        u_str = str(u)
        v_str = str(v)
        
        # Build line: u v [key] [attributes]
        line_parts = [u_str, v_str]
        
        # Add key if there are parallel edges
        if graph._graph.number_of_edges(u, v) > 1:
            line_parts.append(str(key))
        
        # Add attributes
        if data:
            for attr_key, attr_value in data.items():
                if isinstance(attr_value, str):
                    # Escape spaces in string values
                    if ' ' in attr_value:
                        attr_value = f'"{attr_value}"'
                    line_parts.append(f'{attr_key}={attr_value}')
                else:
                    line_parts.append(f'{attr_key}={attr_value}')
        
        lines.append(' '.join(line_parts))
    
    return '\n'.join(lines) + '\n'


def from_edgelist(edgelist_string: str, **kwargs: Any) -> DirectedMultiGraph:
    """
    Parse an edge-list format string and create a DirectedMultiGraph.
    
    Parameters
    ----------
    edgelist_string : str
        Edge-list format string containing graph data.
    **kwargs
        Additional arguments (currently unused, for compatibility).
    
    Returns
    -------
    DirectedMultiGraph
        Parsed directed multi-graph.
    
    Raises
    ------
    PhyloZooParseError
        If the edge-list string is malformed or cannot be parsed.
    
    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph.io import from_edgelist
    >>> 
    >>> el_str = '''1 2
    ... 2 3 weight=2.0
    ... 3 4 0 key1=value1'''
    >>> 
    >>> G = from_edgelist(el_str)
    >>> G.number_of_nodes()
    4
    >>> G.number_of_edges()
    3
    
    Notes
    -----
    This parser expects:
    - One edge per line
    - Format: `u v` or `u v key` or `u v key attr1=value1 attr2=value2`
    - Uses node_id as the label/name
    """
    graph = DirectedMultiGraph()
    
    for line in edgelist_string.strip().split('\n'):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        
        parts = line.split()
        if len(parts) < 2:
            raise PhyloZooParseError(f"Invalid edge line (need at least 2 values): {line}")
        
        u_str = parts[0]
        v_str = parts[1]
        
        # Convert node strings to appropriate types
        u = _convert_node_id(u_str)
        v = _convert_node_id(v_str)
        
        # Parse key and attributes
        key = None
        attrs = {}
        
        if len(parts) > 2:
            # Check if third part is a key (integer) or an attribute
            third_part = parts[2]
            if '=' not in third_part:
                # It's a key
                try:
                    key = int(third_part)
                    start_idx = 3
                except ValueError:
                    start_idx = 2
            else:
                start_idx = 2
            
            # Parse attributes
            for part in parts[start_idx:]:
                if '=' not in part:
                    continue
                
                attr_key, attr_value = part.split('=', 1)
                attr_value = attr_value.strip('"\'')
                
                # Try to convert to appropriate type
                try:
                    if '.' in attr_value:
                        attrs[attr_key] = float(attr_value)
                    else:
                        attrs[attr_key] = int(attr_value)
                except ValueError:
                    attrs[attr_key] = attr_value
        
        # Add edge
        graph.add_edge(u, v, key=key, **attrs)
    
    return graph


# Register format handlers with FormatRegistry
FormatRegistry.register(
    DirectedMultiGraph, 'dot',
    reader=from_dot,
    writer=to_dot,
    extensions=['.dot', '.gv'],
    default=True
)

FormatRegistry.register(
    DirectedMultiGraph, 'edgelist',
    reader=from_edgelist,
    writer=to_edgelist,
    extensions=['.el']
)

