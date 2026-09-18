"""
Directed multi-graph I/O module.

This module provides format handlers for reading and writing directed multi-graphs
to/from files. Format handlers are registered with FormatRegistry for use with
the IOMixin system.

The following format handlers are defined and registered:

- **dot**: DOT format (Graphviz) (extensions: .dot, .gv)
  - Writer: `to_dot()` - Converts DirectedMultiGraph to a Graphviz ``digraph``
  - Reader: `from_dot()` - Parses a DOT string to DirectedMultiGraph
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
The DOT scaffolding (escaping, attribute formatting/parsing, node-id coercion and
document parsing) is shared with the mixed-multigraph formats and lives in
:mod:`phylozoo.utils.io.format_utils.dot`.

DOT format supports:

- Node attributes (label, shape, color, etc.)
- Edge attributes (label, weight, color, etc.)
- Graph attributes
- Parallel edges (multigraph support), encoded with an explicit ``key`` attribute,
  also written for a non-parallel edge whose key is not 0 so that keys round-trip

Edge-list format:

- Simple text format: one edge per line
- Format: `u v` or `u v key` or `u v key attr1=value1 attr2=value2`
- Uses node_id as the label/name
"""

from __future__ import annotations

from typing import Any

from phylozoo.utils.exceptions import PhyloZooParseError
from phylozoo.utils.io import FormatRegistry
from phylozoo.utils.io.format_utils.dot import (
    convert_node_id,
    dot_edge_line,
    dot_node_line,
    escape_dot_string,
    parse_dot_document,
)

from .base import DirectedMultiGraph


def to_dot(graph: DirectedMultiGraph, **kwargs: Any) -> str:
    """
    Convert a DirectedMultiGraph to a DOT format string.

    Parameters
    ----------
    graph : DirectedMultiGraph
        The directed multi-graph to convert.
    **kwargs
        Additional arguments (``graph_name`` is honoured).

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
    The DOT format includes a ``digraph`` declaration, node declarations with
    attributes, ``->`` edge declarations, graph attributes (if any) and parallel
    edges (encoded with a ``key`` attribute, also written whenever a key is not 0).
    """
    lines = []

    graph_name = kwargs.get("graph_name", "")
    lines.append(f"digraph {escape_dot_string(graph_name)} {{" if graph_name else "digraph {")

    # Graph attributes
    if hasattr(graph, "_graph") and hasattr(graph._graph, "graph"):
        for key, value in graph._graph.graph.items():
            value_str = escape_dot_string(value) if isinstance(value, str) else str(value)
            lines.append(f"    {key}={value_str};")

    # Node declarations
    for node in graph.nodes():
        node_attrs = dict(graph._graph.nodes[node]) if node in graph._graph else {}
        lines.append(dot_node_line(node, node_attrs))

    # Edge declarations
    for u, v, key, data in graph.edges_iter(keys=True, data=True):
        edge_attrs = dict(data) if data else {}
        if key != 0 or graph._graph.number_of_edges(u, v) > 1:
            edge_attrs["key"] = key
        lines.append(dot_edge_line(u, v, edge_attrs, arrow="->"))

    lines.append("}")
    return "\n".join(lines) + "\n"


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
        Parsed directed multi-graph (every edge is treated as directed).

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
    """
    graph_attrs, nodes_data, edges_data = parse_dot_document(dot_string)

    graph: Any = DirectedMultiGraph(attributes=graph_attrs if graph_attrs else None)

    for node_id, attrs in nodes_data.items():
        graph.add_node(node_id, **attrs)

    for u, v, key, attrs, _directed in edges_data:
        if u not in graph:
            graph.add_node(u)
        if v not in graph:
            graph.add_node(v)
        graph.add_edge(u, v, key=key, **attrs)

    return graph  # type: ignore[no-any-return]


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
    - Includes edge keys for parallel edges, and for any edge whose key is not 0
    - Includes edge attributes if present
    """
    lines = []

    for u, v, key, data in graph.edges_iter(keys=True, data=True):
        line_parts = [str(u), str(v)]

        # Add key for parallel edges, and for a lone edge carrying a non-zero key
        if key != 0 or graph._graph.number_of_edges(u, v) > 1:
            line_parts.append(str(key))

        # Add attributes
        if data:
            for attr_key, attr_value in data.items():
                if isinstance(attr_value, str) and " " in attr_value:
                    attr_value = f'"{attr_value}"'
                line_parts.append(f"{attr_key}={attr_value}")

        lines.append(" ".join(line_parts))

    return "\n".join(lines) + "\n"


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
    graph: Any = DirectedMultiGraph()

    for line in edgelist_string.strip().split("\n"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        parts = line.split()
        if len(parts) < 2:
            raise PhyloZooParseError(f"Invalid edge line (need at least 2 values): {line}")

        u = convert_node_id(parts[0])
        v = convert_node_id(parts[1])

        key = None
        attrs: dict[str, Any] = {}

        if len(parts) > 2:
            # The third token is a key (a bare integer) or the first attribute.
            third_part = parts[2]
            if "=" not in third_part:
                try:
                    key = int(third_part)
                    start_idx = 3
                except ValueError:
                    start_idx = 2
            else:
                start_idx = 2

            for part in parts[start_idx:]:
                if "=" not in part:
                    continue
                attr_key, attr_value = part.split("=", 1)
                attr_value = attr_value.strip("\"'")
                try:
                    attrs[attr_key] = float(attr_value) if "." in attr_value else int(attr_value)
                except ValueError:
                    attrs[attr_key] = attr_value

        graph.add_edge(u, v, key=key, **attrs)

    return graph  # type: ignore[no-any-return]


# Register format handlers with FormatRegistry
FormatRegistry.register(
    DirectedMultiGraph,
    "dot",
    reader=from_dot,
    writer=to_dot,
    extensions=[".dot", ".gv"],
    default=True,
)

FormatRegistry.register(
    DirectedMultiGraph, "edgelist", reader=from_edgelist, writer=to_edgelist, extensions=[".el"]
)
