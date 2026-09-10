"""
Shared DOT-format scaffolding for the graph primitives.

Three closely related serialisations are built on the same DOT syntax and share
the helpers below (mirroring how :mod:`~phylozoo.utils.io.format_utils.nexus` and
:mod:`~phylozoo.utils.io.format_utils.phylip` are shared between classes):

* ``dot`` for :class:`~phylozoo.core.primitives.d_multigraph.base.DirectedMultiGraph`
  — a standard Graphviz ``digraph`` with ``->`` edges;
* ``dot`` for :class:`~phylozoo.core.primitives.m_multigraph.base.MixedMultiGraph`
  — a standard Graphviz ``digraph`` whose undirected edges carry ``dir=none``
  (Graphviz's own marker for an edge drawn without arrowheads), so the file is
  valid DOT and opens in any Graphviz tool;
* ``phylozoo-dot`` for ``MixedMultiGraph`` — the legacy PhyloZoo dialect using a
  ``graph`` block with both ``->`` (directed) and ``--`` (undirected) edges.

Edge keys are encoded with an explicit ``key=<int>`` attribute. It is written
whenever an edge is parallel *or* its key is not 0, so that a lone edge left
holding a non-zero key (e.g. after a parallel sibling was removed) keeps it
across a round-trip. Edges with the default key 0 are written without it.

:func:`parse_dot_document` parses all three: it accepts either header, treats an
edge as **undirected** when it is written with ``--`` or carries ``dir=none``, and
returns the structural ``dir``/``key`` markers consumed (not as leftover
attributes). The class-specific readers in ``core/primitives/*/io.py`` then build
the appropriate graph from its output.
"""

from __future__ import annotations

import re
from typing import Any

from phylozoo.utils.exceptions import PhyloZooParseError

__all__ = [
    "escape_dot_string",
    "format_dot_attributes",
    "parse_dot_attributes",
    "convert_node_id",
    "dot_node_line",
    "dot_edge_line",
    "parse_dot_document",
]

# Node ids may be word characters / dots (e.g. floats), or a quoted string.
_NODE = r'([\w.]+|"[^"]+")'
_NODE_LINE = re.compile(rf"{_NODE}\s*(?:\[([^\]]+)\])?\s*;")
_UNDIRECTED_EDGE = re.compile(rf"{_NODE}\s*--\s*{_NODE}\s*(?:\[([^\]]+)\])?\s*;")
_DIRECTED_EDGE = re.compile(rf"{_NODE}\s*->\s*{_NODE}\s*(?:\[([^\]]+)\])?\s*;")
_GRAPH_ATTR = re.compile(r"^(\w+)\s*=\s*([^;]+);$")
_HEADER = re.compile(r'(?:di)?graph\s+(?:\w+|"[^"]+")?\s*\{', re.IGNORECASE)


def escape_dot_string(s: str) -> str:
    """Escape a string for DOT, quoting it when it contains special characters."""
    if any(c in s for c in [" ", "\t", "\n", '"', "\\", "[", "]", "{", "}", "-", ">"]):
        s = s.replace("\\", "\\\\")
        s = s.replace('"', '\\"')
        return f'"{s}"'
    return s


def format_dot_attributes(attrs: dict[str, Any]) -> str:
    """Format an attribute mapping as ``[key1=value1, key2=value2]`` (``""`` if empty)."""
    if not attrs:
        return ""

    parts = []
    for key, value in attrs.items():
        if isinstance(value, bool):
            value_str = "true" if value else "false"
        elif isinstance(value, str):
            value_str = escape_dot_string(value)
        elif isinstance(value, (int, float)):
            value_str = str(value)
        else:
            value_str = escape_dot_string(str(value))
        parts.append(f"{key}={value_str}")

    return "[" + ", ".join(parts) + "]"


def parse_dot_attributes(attrs_str: str) -> dict[str, Any]:
    """
    Parse a DOT attribute string like ``key1=value1, key2=value2`` into a dict.

    Values are coerced to ``bool``/``int``/``float`` where possible, except a
    ``label``, which is always kept as a **string** — a label is a display string,
    and coercing e.g. ``label=1`` to an ``int`` would leave nodes with non-string
    labels that downstream label validation rejects.
    """
    attrs: dict[str, Any] = {}
    if not attrs_str.strip():
        return attrs

    # Split on commas, but respect quoted strings.
    parts: list[str] = []
    current = ""
    in_quotes = False
    escape_next = False
    for char in attrs_str:
        if escape_next:
            current += char
            escape_next = False
            continue
        if char == "\\":
            escape_next = True
            current += char
            continue
        if char in ('"', "'"):
            in_quotes = not in_quotes
            current += char
            continue
        if char == "," and not in_quotes:
            parts.append(current.strip())
            current = ""
        else:
            current += char
    if current.strip():
        parts.append(current.strip())

    for part in parts:
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        key = key.strip()
        value = value.strip().strip("\"'")

        if key == "label":
            attrs[key] = value
        elif value.lower() == "true":
            attrs[key] = True
        elif value.lower() == "false":
            attrs[key] = False
        else:
            try:
                attrs[key] = float(value) if "." in value else int(value)
            except ValueError:
                attrs[key] = value

    return attrs


def convert_node_id(node_str: str) -> Any:
    """Coerce a node-id string to ``int``, then ``float``, else keep it a string."""
    try:
        return int(node_str)
    except ValueError:
        pass
    try:
        return float(node_str)
    except ValueError:
        pass
    return node_str


def dot_node_line(node: Any, attrs: dict[str, Any], *, indent: str = "    ") -> str:
    """Render one ``node [attrs];`` declaration (bare ``node;`` if no attributes)."""
    node_str = escape_dot_string(str(node))
    attrs_str = format_dot_attributes(attrs)
    return f"{indent}{node_str} {attrs_str};" if attrs_str else f"{indent}{node_str};"


def dot_edge_line(
    u: Any, v: Any, attrs: dict[str, Any], *, arrow: str = "->", indent: str = "    "
) -> str:
    """Render one ``u <arrow> v [attrs];`` edge declaration (``arrow`` is ``->`` or ``--``)."""
    u_str = escape_dot_string(str(u))
    v_str = escape_dot_string(str(v))
    attrs_str = format_dot_attributes(attrs)
    if attrs_str:
        return f"{indent}{u_str} {arrow} {v_str} {attrs_str};"
    return f"{indent}{u_str} {arrow} {v_str};"


def _pop_key(attrs: dict[str, Any]) -> int | None:
    """Pop a parallel-edge ``key`` from ``attrs`` and return it as an int (or None)."""
    if "key" in attrs:
        try:
            return int(attrs.pop("key"))
        except (ValueError, TypeError):
            attrs.pop("key", None)
    return None


def parse_dot_document(
    dot_string: str,
) -> tuple[
    dict[str, Any],
    dict[Any, dict[str, Any]],
    list[tuple[Any, Any, int | None, dict[str, Any], bool]],
]:
    """
    Parse a DOT / phylozoo-dot document into its components.

    Parameters
    ----------
    dot_string : str
        A ``graph`` or ``digraph`` document.

    Returns
    -------
    graph_attrs : dict
        Graph-level attributes.
    nodes_data : dict
        Maps each declared node id to its attribute dict.
    edges_data : list
        One ``(u, v, key, attrs, directed)`` tuple per edge. ``directed`` is
        ``False`` when the edge used ``--`` or carried ``dir=none``; the ``dir``
        and ``key`` markers are consumed and never appear in ``attrs``.

    Raises
    ------
    PhyloZooParseError
        If no graph declaration is found or the braces are unbalanced.
    """
    # Strip // line comments.
    content = "\n".join(
        line[: line.index("//")] if "//" in line else line for line in dot_string.split("\n")
    )

    if not _HEADER.search(content):
        raise PhyloZooParseError("Could not find graph declaration in the DOT string")

    # Extract the body between the outermost matching braces.
    start_idx = content.index("{")
    brace = 0
    end_idx = start_idx
    for i, char in enumerate(content[start_idx:], start=start_idx):
        if char == "{":
            brace += 1
        elif char == "}":
            brace -= 1
            if brace == 0:
                end_idx = i
                break
    if brace != 0:
        raise PhyloZooParseError("Unmatched braces in the DOT string")
    body = content[start_idx + 1 : end_idx]

    graph_attrs: dict[str, Any] = {}
    nodes_data: dict[Any, dict[str, Any]] = {}
    edges_data: list[tuple[Any, Any, int | None, dict[str, Any], bool]] = []

    for line in body.split("\n"):
        line = line.strip()
        if not line or line.startswith("//") or line.startswith("#"):
            continue

        graph_attr = _GRAPH_ATTR.match(line)
        if graph_attr:
            graph_attrs[graph_attr.group(1).strip()] = graph_attr.group(2).strip().strip("\"'")
            continue

        undirected = _UNDIRECTED_EDGE.search(line)
        if undirected:
            u = convert_node_id(undirected.group(1).strip("\"'"))
            v = convert_node_id(undirected.group(2).strip("\"'"))
            attrs = parse_dot_attributes(undirected.group(3) or "")
            key = _pop_key(attrs)
            attrs.pop("dir", None)
            edges_data.append((u, v, key, attrs, False))
            continue

        directed = _DIRECTED_EDGE.search(line)
        if directed:
            u = convert_node_id(directed.group(1).strip("\"'"))
            v = convert_node_id(directed.group(2).strip("\"'"))
            attrs = parse_dot_attributes(directed.group(3) or "")
            key = _pop_key(attrs)
            is_directed = attrs.pop("dir", None) != "none"
            edges_data.append((u, v, key, attrs, is_directed))
            continue

        node = _NODE_LINE.search(line)
        if node:
            node_id = convert_node_id(node.group(1).strip("\"'"))
            nodes_data[node_id] = parse_dot_attributes(node.group(2) or "")
            continue

    return graph_attrs, nodes_data, edges_data
