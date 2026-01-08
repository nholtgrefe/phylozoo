"""
Base layout classes for phylogenetic networks.

This module provides the base Layout class that can be used for any graph type.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, TypeVar

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
    """
    is_directed: bool
    is_hybrid: bool
    is_parallel: bool


@dataclass(frozen=True)
class EdgeRoute:
    """
    Edge routing information for layout.

    Attributes
    ----------
    edge_type : EdgeType
        Type of edge with three dimensions (directed/undirected, hybrid/tree, parallel/unique).
    points : tuple[tuple[float, float], ...]
        Polyline points for the edge route. For tree edges, this is
        typically [parent_pos, (mid_x, parent_y), (mid_x, child_y), child_pos].
        For hybrid edges, includes curve control points.
    curve_control : tuple[float, float] | None
        Control point for curved inlet (hybrid edges only).
        None for tree edges.

    Examples
    --------
    >>> route = EdgeRoute(
    ...     edge_type=EdgeType(is_directed=True, is_hybrid=False, is_parallel=False),
    ...     points=((0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 2.0)),
    ...     curve_control=None
    ... )
    >>> route.edge_type.is_directed
    True
    """

    edge_type: EdgeType
    points: tuple[tuple[float, float], ...]
    curve_control: tuple[float, float] | None = None


@dataclass(frozen=True)
class Layout:
    """
    Base layout class for any graph type.

    Immutable value object storing node positions and edge routing information
    in a renderer-agnostic format. Can be used with matplotlib, pyqtgraph,
    ASCII, SVG, TikZ, etc.

    The layout separates geometry from graph structure - all routing
    information is pure geometry and does not modify the graph.

    Attributes
    ----------
    graph : Any
        The graph object this layout is for (read-only).
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
    >>> from phylozoo.visualize.layout import Layout, EdgeRoute, EdgeType
    >>> graph = ...  # Some graph object
    >>> positions = {1: (0.0, 0.0), 2: (1.0, 1.0)}
    >>> edge_routes = {(1, 2, 0): EdgeRoute(EdgeType(is_directed=True, is_hybrid=False, is_parallel=False), ((0.0, 0.0), (1.0, 1.0)))}
    >>> layout = Layout(graph, positions, edge_routes)
    >>> layout.get_position(1)
    (0.0, 0.0)
    """

    graph: Any
    positions: dict[T, tuple[float, float]]
    edge_routes: dict[tuple[T, T, int], EdgeRoute]
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
        >>> layout = Layout(graph, positions, edge_routes)
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
        >>> layout = Layout(graph, positions, edge_routes)
        >>> route = layout.get_edge_route(1, 2, 0)
        >>> route is not None
        True
        """
        return self.edge_routes.get((u, v, key))

    def to_dict(self) -> dict[str, Any]:
        """
        Convert layout to dictionary (for serialization).

        Returns
        -------
        dict[str, Any]
            Dictionary representation of the layout.

        Examples
        --------
        >>> layout = Layout(graph, positions, edge_routes)
        >>> data = layout.to_dict()
        >>> 'positions' in data
        True
        """
        return {
            'positions': self.positions,
            'edge_routes': {
                f"{u}_{v}_{key}": {
                    'edge_type': {
                        'is_directed': route.edge_type.is_directed,
                        'is_hybrid': route.edge_type.is_hybrid,
                        'is_parallel': route.edge_type.is_parallel,
                    },
                    'points': route.points,
                    'curve_control': route.curve_control,
                }
                for (u, v, key), route in self.edge_routes.items()
            },
            'algorithm': self.algorithm,
            'parameters': self.parameters,
        }

    def __repr__(self) -> str:
        """Return string representation."""
        return (
            f"Layout(graph={self.graph}, "
            f"nodes={len(self.positions)}, edges={len(self.edge_routes)}, "
            f"algorithm={self.algorithm})"
        )

