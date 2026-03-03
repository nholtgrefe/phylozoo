"""
Internal type definitions for viz module.

This module provides type definitions for edge routing, edge types,
layout base classes, and other visualization primitives.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Generic, TypeVar

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
    of the rendering backend. It contains the path points. Curved edges
    (e.g. parallel edges) are computed at render time from the points.

    Attributes
    ----------
    edge_type : EdgeType
        Type of edge with three dimensions (directed/undirected, hybrid/tree, parallel/unique).
    points : tuple[tuple[float, float], ...]
        Polyline points for the edge route. Typically (parent_pos, child_pos)
        for straight edges.

    Examples
    --------
    >>> route = EdgeRoute(
    ...     edge_type=EdgeType(is_directed=True, is_hybrid=False, is_parallel=False),
    ...     points=((0.0, 0.0), (1.0, 1.0)),
    ... )
    >>> route.edge_type.is_directed
    True
    """

    edge_type: EdgeType
    points: tuple[tuple[float, float], ...]


# Type aliases for clarity
NodePosition = tuple[float, float]
Positions = dict[T, NodePosition]
EdgeKey = tuple[T, T, int]
EdgeRoutes = dict[EdgeKey, EdgeRoute]


@dataclass(frozen=True)
class Layout(Generic[T]):
    """
    Base layout class for all network and graph layouts.

    This is the common base class that all specific layout classes inherit from.
    It provides the core interface for storing computed positions and edge routes.

    Attributes
    ----------
    network : Any
        The network or graph this layout is for (read-only).
    positions : dict[T, tuple[float, float]]
        Node positions mapping node ID to (x, y) coordinates (read-only).
    edge_routes : dict[tuple[T, T, int], EdgeRoute]
        Edge routing information mapping (u, v, key) to EdgeRoute (read-only).
    algorithm : str
        Name of the layout algorithm used (read-only).
    parameters : dict[str, Any]
        Parameters used to generate this layout (read-only).

    Examples
    --------
    >>> from phylozoo.viz._types import Layout
    >>> layout = Layout(
    ...     network=network,
    ...     positions={1: (0.0, 0.0), 2: (1.0, 1.0)},
    ...     edge_routes={},
    ...     algorithm='test'
    ... )
    >>> layout.get_position(1)
    (0.0, 0.0)
    """

    network: Any
    positions: Positions[T]
    edge_routes: EdgeRoutes[T]
    algorithm: str = 'unknown'
    parameters: dict[str, Any] = field(default_factory=dict)

    def get_position(self, node: T) -> tuple[float, float]:
        """
        Get position of a node.

        Parameters
        ----------
        node : T
            Node ID.

        Returns
        -------
        tuple[float, float]
            (x, y) position.

        Raises
        ------
        KeyError
            If node is not in the layout.

        Examples
        --------
        >>> layout = Layout(network=net, positions={1: (0.0, 0.0)}, edge_routes={})
        >>> x, y = layout.get_position(1)
        >>> isinstance(x, float)
        True
        """
        return self.positions[node]

    def get_edge_route(
        self, u: T, v: T, key: int = 0
    ) -> EdgeRoute | None:
        """
        Get routing information for an edge.

        Parameters
        ----------
        u : T
            Source node.
        v : T
            Target node.
        key : int, optional
            Edge key (for parallel edges). By default 0.

        Returns
        -------
        EdgeRoute | None
            Edge routing information, or None if edge not found.

        Examples
        --------
        >>> layout = Layout(network=net, positions={}, edge_routes={})
        >>> route = layout.get_edge_route(1, 2, 0)
        >>> route is None
        True
        """
        return self.edge_routes.get((u, v, key))

    def __repr__(self) -> str:
        """Return string representation."""
        return (
            f"{self.__class__.__name__}(network={self.network}, "
            f"nodes={len(self.positions)}, edges={len(self.edge_routes)}, "
            f"algorithm={self.algorithm})"
        )
