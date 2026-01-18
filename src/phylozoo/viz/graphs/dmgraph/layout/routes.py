"""
Edge route computation for DirectedMultiGraph.

This module computes edge routing information for DirectedMultiGraph,
handling parallel edges with curved paths.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, TypeVar

from phylozoo.viz._types import EdgeRoute, EdgeType

if TYPE_CHECKING:
    from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph

T = TypeVar('T')


def compute_dmgraph_routes(
    graph: 'DirectedMultiGraph[T]',
    positions: dict[T, tuple[float, float]],
) -> dict[tuple[T, T, int], EdgeRoute]:
    """
    Compute edge routes for a DirectedMultiGraph, handling parallel edges.

    Parameters
    ----------
    graph : DirectedMultiGraph
        The graph.
    positions : dict[T, tuple[float, float]]
        Node positions.

    Returns
    -------
    dict[tuple[T, T, int], EdgeRoute]
        Edge routing information mapping (u, v, key) to EdgeRoute.

    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
    >>> G = DirectedMultiGraph(edges=[(1, 2), (1, 2)])  # Parallel edges
    >>> positions = {1: (0.0, 0.0), 2: (1.0, 0.0)}
    >>> routes = compute_dmgraph_routes(G, positions)
    >>> len(routes)
    2
    """
    edge_routes: dict[tuple[T, T, int], EdgeRoute] = {}

    # Get edges
    edges = list(graph._graph.edges(keys=True))

    # Group edges by (u, v) to detect parallel edges
    edge_groups: dict[tuple[T, T], list[tuple[T, T, int]]] = {}
    for u, v, key in edges:
        edge_key = (u, v)
        if edge_key not in edge_groups:
            edge_groups[edge_key] = []
        edge_groups[edge_key].append((u, v, key))

    # Compute routes for each edge
    for (u, v), edge_list in edge_groups.items():
        num_parallel = len(edge_list)
        
        for idx, (u_edge, v_edge, key) in enumerate(edge_list):
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
                    is_directed=True,
                    is_hybrid=False,  # Graphs don't have hybrid edges
                    is_parallel=is_parallel,
                ),
                points=points,
                curve_control=curve_control,
            )

    return edge_routes
