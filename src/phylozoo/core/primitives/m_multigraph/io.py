"""
Mixed multi-graph I/O module.

This module provides format handlers for reading and writing mixed multi-graphs
to/from files. Format handlers are registered with FormatRegistry for use with
the IOMixin system.

The following format handlers are defined and registered:

- **phylozoo-dot**: PhyloZoo DOT dialect (extensions: .pzdot, *default*).
  A ``graph`` block with both ``->`` (directed) and ``--`` (undirected) edges.
  Compact and human-readable, but **not** valid Graphviz DOT (a single ``graph``
  block may not contain ``->``).
- **dot**: standard Graphviz DOT (extensions: .dot, .gv).
  A ``digraph`` whose undirected edges carry ``dir=none`` (Graphviz's marker for
  an edge drawn without arrowheads). This *is* valid DOT, so the file opens in any
  Graphviz tool, and it round-trips losslessly (``dir=none`` <-> undirected).

Both formats encode edge keys with an explicit ``key`` attribute -- written for
parallel edges and for any edge whose key is not 0, so keys survive a round-trip -- and
share their reader and all low-level scaffolding with
:mod:`phylozoo.utils.io.format_utils.dot` (which is also used by
:class:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph`).

MixedMultiGraph inherits from IOMixin, so you can use:

- `graph.save('file.pzdot')` / `graph.save('file.dot', format='dot')`
- `graph.load('file.pzdot')` / `graph.load('file.dot')`
- `graph.to_string(format='phylozoo-dot')` / `graph.to_string(format='dot')`
- `graph.from_string(string, format='dot')`
- `MixedMultiGraph.convert('in.pzdot', 'out.dot')` - Convert between formats
"""

from __future__ import annotations

from typing import Any

from phylozoo.utils.io import FormatRegistry
from phylozoo.utils.io.format_utils.dot import (
    dot_edge_line,
    dot_node_line,
    escape_dot_string,
    parse_dot_document,
)

from .base import MixedMultiGraph


def _graph_attr_lines(graph: MixedMultiGraph) -> list[str]:
    """Render the graph-level attribute declarations."""
    lines: list[str] = []
    if hasattr(graph, "_directed") and hasattr(graph._directed, "graph"):
        for key, value in graph._directed.graph.items():
            value_str = escape_dot_string(value) if isinstance(value, str) else str(value)
            lines.append(f"    {key}={value_str};")
    return lines


def _node_attrs(graph: MixedMultiGraph, node: Any) -> dict[str, Any]:
    """Collect a node's attributes from both the directed and undirected subgraphs."""
    attrs: dict[str, Any] = {}
    if hasattr(graph, "_directed") and node in graph._directed:
        attrs.update(graph._directed.nodes[node])
    if hasattr(graph, "_undirected") and node in graph._undirected:
        attrs.update(graph._undirected.nodes[node])
    return attrs


def to_phylozoo_dot(graph: MixedMultiGraph, **kwargs: Any) -> str:
    """
    Convert a MixedMultiGraph to a phylozoo-dot string.

    Parameters
    ----------
    graph : MixedMultiGraph
        The mixed multi-graph to convert.
    **kwargs
        Additional arguments (``graph_name`` is honoured).

    Returns
    -------
    str
        The phylozoo-dot representation: a ``graph`` block using ``--`` for
        undirected edges and ``->`` for directed edges.

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
    """
    graph_name = kwargs.get("graph_name", "")
    lines = [f"graph {escape_dot_string(graph_name)} {{" if graph_name else "graph {"]
    lines.extend(_graph_attr_lines(graph))

    for node in set(graph.nodes()):
        lines.append(dot_node_line(node, _node_attrs(graph, node)))

    for u, v, key, data in graph.undirected_edges_iter(keys=True, data=True):
        edge_attrs = dict(data) if data else {}
        if key != 0 or graph._undirected.number_of_edges(u, v) > 1:
            edge_attrs["key"] = key
        lines.append(dot_edge_line(u, v, edge_attrs, arrow="--"))

    for u, v, key, data in graph.directed_edges_iter(keys=True, data=True):
        edge_attrs = dict(data) if data else {}
        if key != 0 or graph._directed.number_of_edges(u, v) > 1:
            edge_attrs["key"] = key
        lines.append(dot_edge_line(u, v, edge_attrs, arrow="->"))

    lines.append("}")
    return "\n".join(lines) + "\n"


def to_dot(graph: MixedMultiGraph, **kwargs: Any) -> str:
    """
    Convert a MixedMultiGraph to a standard Graphviz DOT string.

    Undirected edges are written as ``u -> v [dir=none]`` inside a ``digraph``,
    which any Graphviz tool parses and draws without arrowheads. The result
    round-trips losslessly (``dir=none`` is read back as an undirected edge).

    Parameters
    ----------
    graph : MixedMultiGraph
        The mixed multi-graph to convert.
    **kwargs
        Additional arguments (``graph_name`` is honoured).

    Returns
    -------
    str
        The DOT representation.

    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
    >>> from phylozoo.core.primitives.m_multigraph.io import to_dot
    >>>
    >>> G = MixedMultiGraph()
    >>> G.add_directed_edge(1, 2)
    0
    >>> G.add_undirected_edge(2, 3)
    0
    >>> dot_str = to_dot(G)
    >>> 'digraph' in dot_str
    True
    >>> '1 -> 2' in dot_str
    True
    >>> '2 -> 3 [dir=none]' in dot_str
    True
    """
    graph_name = kwargs.get("graph_name", "")
    lines = [f"digraph {escape_dot_string(graph_name)} {{" if graph_name else "digraph {"]
    lines.extend(_graph_attr_lines(graph))

    for node in set(graph.nodes()):
        lines.append(dot_node_line(node, _node_attrs(graph, node)))

    for u, v, key, data in graph.directed_edges_iter(keys=True, data=True):
        edge_attrs = dict(data) if data else {}
        if key != 0 or graph._directed.number_of_edges(u, v) > 1:
            edge_attrs["key"] = key
        lines.append(dot_edge_line(u, v, edge_attrs, arrow="->"))

    for u, v, key, data in graph.undirected_edges_iter(keys=True, data=True):
        edge_attrs = dict(data) if data else {}
        if key != 0 or graph._undirected.number_of_edges(u, v) > 1:
            edge_attrs["key"] = key
        edge_attrs["dir"] = "none"
        lines.append(dot_edge_line(u, v, edge_attrs, arrow="->"))

    lines.append("}")
    return "\n".join(lines) + "\n"


def _read_mixed_dot(dot_string: str) -> MixedMultiGraph:
    """Build a MixedMultiGraph from any DOT/phylozoo-dot document.

    The shared :func:`~phylozoo.utils.io.format_utils.dot.parse_dot_document`
    recognises undirected edges written either as ``--`` (phylozoo-dot) or as
    ``-> [dir=none]`` (standard DOT), so a single reader serves both formats.
    """
    graph_attrs, nodes_data, edges_data = parse_dot_document(dot_string)

    graph: Any = MixedMultiGraph(
        attributes=graph_attrs if graph_attrs else None, directed_edges=None, undirected_edges=None
    )

    for node_id, attrs in nodes_data.items():
        graph.add_node(node_id, **attrs)

    for u, v, key, attrs, directed in edges_data:
        if u not in graph:
            graph.add_node(u)
        if v not in graph:
            graph.add_node(v)
        if directed:
            graph.add_directed_edge(u, v, key=key, **attrs)
        else:
            graph.add_undirected_edge(u, v, key=key, **attrs)

    return graph  # type: ignore[no-any-return]


def from_phylozoo_dot(pzdot_string: str, **kwargs: Any) -> MixedMultiGraph:
    """
    Parse a phylozoo-dot string and create a MixedMultiGraph.

    Parameters
    ----------
    pzdot_string : str
        PhyloZoo DOT string (a ``graph`` block with ``->`` and ``--`` edges).
    **kwargs
        Additional arguments (currently unused, for compatibility).

    Returns
    -------
    MixedMultiGraph
        Parsed mixed multi-graph.

    Raises
    ------
    PhyloZooParseError
        If the string is malformed or cannot be parsed.

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
    """
    return _read_mixed_dot(pzdot_string)


def from_dot(dot_string: str, **kwargs: Any) -> MixedMultiGraph:
    """
    Parse a standard Graphviz DOT string and create a MixedMultiGraph.

    Edges carrying ``dir=none`` become undirected; all other ``->`` edges become
    directed.

    Parameters
    ----------
    dot_string : str
        DOT string (a ``digraph`` block; undirected edges marked ``dir=none``).
    **kwargs
        Additional arguments (currently unused, for compatibility).

    Returns
    -------
    MixedMultiGraph
        Parsed mixed multi-graph.

    Raises
    ------
    PhyloZooParseError
        If the string is malformed or cannot be parsed.

    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph.io import from_dot
    >>>
    >>> dot_str = '''digraph {
    ...     1 -> 2 [weight=1.0];
    ...     2 -> 3 [dir=none];
    ... }'''
    >>>
    >>> G = from_dot(dot_str)
    >>> len(list(G.directed_edges_iter()))
    1
    >>> len(list(G.undirected_edges_iter()))
    1
    """
    return _read_mixed_dot(dot_string)


# Register format handlers with FormatRegistry
FormatRegistry.register(
    MixedMultiGraph,
    "phylozoo-dot",
    reader=from_phylozoo_dot,
    writer=to_phylozoo_dot,
    extensions=[".pzdot"],
    default=True,
)

FormatRegistry.register(
    MixedMultiGraph,
    "dot",
    reader=from_dot,
    writer=to_dot,
    extensions=[".dot", ".gv"],
)
