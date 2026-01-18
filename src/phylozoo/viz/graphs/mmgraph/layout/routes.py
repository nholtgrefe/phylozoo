"""
Edge route computation for MixedMultiGraph.

This module computes edge routing information for MixedMultiGraph,
handling parallel edges with curved paths.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, TypeVar

from phylozoo.viz._types import EdgeRoute, EdgeType

if TYPE_CHECKING:
    from phylozoo.core.primitives.m_multigraph import MixedMultiGraph

T = TypeVar('T')


def compute_mmgraph_routes(
    graph: 'MixedMultiGraph[T]',
    positions: dict[T, tuple[float, float]],
) -> dict[tuple[T, T, int], EdgeRoute]:
    """
    Compute edge routes for a MixedMultiGraph, handling parallel edges.

    Parameters
    ----------
    graph : MixedMultiGraph
        The graph.
    positions : dict[T, tuple[float, float]]
        Node positions.

    Returns
    -------
    dict[tuple[T, T, int], EdgeRoute]
        Edge routing information mapping (u, v, key) to EdgeRoute.

    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
    >>> G = MixedMultiGraph(directed_edges=[(1, 2)], undirected_edges=[(2, 3)])
    >>> positions = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (2.0, 0.0)}
    >>> routes = compute_mmgraph_routes(G, positions)
    >>> len(routes)
    2
    """
    edge_routes: dict[tuple[T, T, int], EdgeRoute] = {}

    # Get edges from both directed and undirected
    edges = []
    for u, v, key in graph._directed.edges(keys=True):
        edges.append((u, v, key, True))  # True = directed
    for u, v, key in graph._undirected.edges(keys=True):
        edges.append((u, v, key, False))  # False = undirected

    # Group edges by (u, v) to detect parallel edges
    # For undirected edges, use (min, max) as key
    edge_groups: dict[tuple[T, T], list[tuple[T, T, int, bool]]] = {}
    for u, v, key, directed in edges:
        edge_key = (u, v) if directed else (min(u, v), max(u, v))
        if edge_key not in edge_groups:
            edge_groups[edge_key] = []
        edge_groups[edge_key].append((u, v, key, directed))

    # Compute routes for each edge
    for (u, v), edge_list in edge_groups.items():
        num_parallel = len(edge_list)
        
        for idx, (u_edge, v_edge, key, directed) in enumerate(edge_list):
            if u_edge not in positions or v_edge not in positions:
                continue

            is_parallel = num_parallel > 1
            
            # For parallel edges, the backend will compute the curve
            # So we just provide start and end points
            x1, y1 = positions[u_edge]
            x2, y2 = positions[v_edge]
            points = ((x1, y1), (x2, y2))
            curve_control = None

            edge_routes[(u_edge, v_edge, key)] = EdgeRoute(
                edge_type=EdgeType(
                    is_directed=directed,
                    is_hybrid=False,  # Graphs don't have hybrid edges
                    is_parallel=is_parallel,
                ),
                points=points,
                curve_control=curve_control,
            )

    return edge_routes
