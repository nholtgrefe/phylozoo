"""
Type definitions for visualize2 module.

This module provides backend-agnostic type definitions for edge routing,
edge types, and other visualization primitives.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeVar

T = TypeVar('T')


@dataclass(frozen=True)
class EdgeType:
    """
    Edge type with three dimensions.

    Attributes
    ----------
    is_directed : bool
        True if edge is directed, False if undirected.
    is_hybrid : bool
        True if edge is a hybrid edge, False if tree edge.
    is_parallel : bool
        True if edge is part of a parallel edge set, False if unique.

    Examples
    --------
    >>> edge_type = EdgeType(is_directed=True, is_hybrid=False, is_parallel=False)
    >>> edge_type.is_directed
    True
    """

    is_directed: bool
    is_hybrid: bool
    is_parallel: bool


@dataclass(frozen=True)
class EdgeRoute:
    """
    Edge routing information for layout.

    This class stores geometric routing information for edges, independent
    of the rendering backend. It contains the path points and optional
    control points for curved edges.

    Attributes
    ----------
    edge_type : EdgeType
        Type of edge with three dimensions (directed/undirected, hybrid/tree, parallel/unique).
    points : tuple[tuple[float, float], ...]
        Polyline points for the edge route. For tree edges, this is
        typically [parent_pos, child_pos]. For hybrid edges, includes
        curve control points.
    curve_control : tuple[float, float] | None
        Control point for curved edges (hybrid edges only).
        None for straight tree edges.

    Examples
    --------
    >>> route = EdgeRoute(
    ...     edge_type=EdgeType(is_directed=True, is_hybrid=False, is_parallel=False),
    ...     points=((0.0, 0.0), (1.0, 1.0)),
    ...     curve_control=None
    ... )
    >>> route.edge_type.is_directed
    True
    """

    edge_type: EdgeType
    points: tuple[tuple[float, float], ...]
    curve_control: tuple[float, float] | None = None


# Type aliases for clarity
NodePosition = tuple[float, float]
Positions = dict[T, NodePosition]
EdgeKey = tuple[T, T, int]
EdgeRoutes = dict[EdgeKey, EdgeRoute]

