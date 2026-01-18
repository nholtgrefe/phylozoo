"""
Edge route computation for DirectedPhyNetwork.

This module computes edge routing information (paths, control points) for
both backbone tree edges and hybrid edges. The routes are backend-agnostic
and contain only geometric information.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, TypeVar

from phylozoo.viz._types import EdgeRoute, EdgeType

if TYPE_CHECKING:
    from phylozoo.core.network.dnetwork import DirectedPhyNetwork

T = TypeVar('T')


def compute_backbone_routes(
    network: 'DirectedPhyNetwork',
    positions: dict[T, tuple[float, float]],
    backbone_edges: set[tuple[T, T, int]],
) -> dict[tuple[T, T, int], EdgeRoute]:
    """
    Route backbone tree edges as straight lines.

    Parameters
    ----------
    network : DirectedPhyNetwork
        The network.
    positions : dict[T, tuple[float, float]]
        Node positions.
    backbone_edges : set[tuple[T, T, int]]
        Backbone tree edges to route.

    Returns
    -------
    dict[tuple[T, T, int], EdgeRoute]
        Edge routing information for backbone edges.

    Examples
    --------
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)])
    >>> positions = {3: (0.0, 0.0), 1: (0.0, 1.0), 2: (1.0, 1.0)}
    >>> routes = compute_backbone_routes(net, positions, {(3, 1, 0), (3, 2, 0)})
    >>> len(routes)
    2
    """
    edge_routes: dict[tuple[T, T, int], EdgeRoute] = {}

    # Count parallel edges
    parallel_counts: dict[tuple[T, T], int] = {}
    for u, v, key in backbone_edges:
        edge_pair = (u, v)
        parallel_counts[edge_pair] = parallel_counts.get(edge_pair, 0) + 1

    for u, v, key in backbone_edges:
        if u in positions and v in positions:
            is_parallel = parallel_counts[(u, v)] > 1

            # Straight line route
            points = (positions[u], positions[v])

            edge_routes[(u, v, key)] = EdgeRoute(
                edge_type=EdgeType(
                    is_directed=True,
                    is_hybrid=False,
                    is_parallel=is_parallel,
                ),
                points=points,
                curve_control=None,
            )

    return edge_routes


def compute_hybrid_routes(
    network: 'DirectedPhyNetwork',
    positions: dict[T, tuple[float, float]],
    reticulate_edges: set[tuple[T, T, int]],
) -> dict[tuple[T, T, int], EdgeRoute]:
    """
    Route hybrid edges as straight lines.

    Parameters
    ----------
    network : DirectedPhyNetwork
        The network.
    positions : dict[T, tuple[float, float]]
        Node positions.
    reticulate_edges : set[tuple[T, T, int]]
        Hybrid edges to route.

    Returns
    -------
    dict[tuple[T, T, int], EdgeRoute]
        Edge routing information for hybrid edges.

    Examples
    --------
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> net = DirectedPhyNetwork(edges=[(3, 2), (4, 2)])
    >>> positions = {3: (0.0, 0.0), 4: (1.0, 0.0), 2: (0.5, 1.0)}
    >>> routes = compute_hybrid_routes(net, positions, {(4, 2, 0)})
    >>> len(routes)
    1
    """
    edge_routes: dict[tuple[T, T, int], EdgeRoute] = {}

    # Count parallel edges
    parallel_counts: dict[tuple[T, T], int] = {}
    for u, v, key in reticulate_edges:
        edge_pair = (u, v)
        parallel_counts[edge_pair] = parallel_counts.get(edge_pair, 0) + 1

    for u, v, key in reticulate_edges:
        if u not in positions or v not in positions:
            continue

        is_parallel = parallel_counts[(u, v)] > 1

        # Straight line route (no curves)
        points = (positions[u], positions[v])

        edge_routes[(u, v, key)] = EdgeRoute(
            edge_type=EdgeType(
                is_directed=True,
                is_hybrid=True,
                is_parallel=is_parallel,
            ),
            points=points,
            curve_control=None,
        )

    return edge_routes
