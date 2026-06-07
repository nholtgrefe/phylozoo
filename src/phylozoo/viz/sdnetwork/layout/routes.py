"""
Edge route computation for SemiDirectedPhyNetwork.

This module computes edge routing information for radial layouts.
The routes are backend-agnostic and contain only geometric information.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, TypeVar

from phylozoo.viz._types import EdgeRoute, EdgeType

if TYPE_CHECKING:
    from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork

T = TypeVar("T")


def compute_radial_routes(
    network: "SemiDirectedPhyNetwork",
    positions: dict[T, tuple[float, float]],
) -> dict[tuple[T, T, int], EdgeRoute]:
    """
    Route edges for radial layout as straight lines.

    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The network.
    positions : dict[T, tuple[float, float]]
        Node positions.

    Returns
    -------
    dict[tuple[T, T, int], EdgeRoute]
        Edge routing information for all edges.

    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> positions = {3: (0.0, 0.0), 1: (1.0, 0.0), 2: (0.0, 1.0)}
    >>> routes = compute_radial_routes(net, positions)
    >>> len(routes)
    2
    """
    edge_routes: dict[tuple[T, T, int], EdgeRoute] = {}

    # Only route edges between nodes that are in the original network
    # (to avoid routing to/from subdivision nodes created during conversion)
    original_nodes = set(network._graph.nodes)

    # Count parallel edges
    parallel_counts: dict[tuple[T, T], int] = {}
    for u, v, key in network._graph.edges(keys=True):
        edge_pair = (u, v)
        parallel_counts[edge_pair] = parallel_counts.get(edge_pair, 0) + 1

    # Process all edges (both directed and undirected)
    # Only route edges between nodes in the original network
    for u, v, key in network._graph.edges(keys=True):
        # Skip edges involving nodes not in the original network
        if u not in original_nodes or v not in original_nodes:
            continue
        # Skip edges where positions aren't available
        if u not in positions or v not in positions:
            continue

        is_parallel = parallel_counts[(u, v)] > 1

        # The combined (undirected) MultiGraph may yield a directed edge in
        # the reverse of the direction stored in _directed (e.g. (6,5) when
        # _directed has (5,6)).  Check both orientations.
        is_directed_fwd = network._graph._directed.has_edge(u, v, key=key)
        is_directed_rev = network._graph._directed.has_edge(v, u, key=key)
        is_directed = is_directed_fwd or is_directed_rev

        # Ensure points run parent → hybrid so the arrowhead points correctly.
        if is_directed_rev:
            src, dst = v, u
        else:
            src, dst = u, v

        is_hybrid = dst in network.hybrid_nodes if is_directed else False

        # Straight line route
        points = (positions[src], positions[dst])

        edge_routes[(u, v, key)] = EdgeRoute(
            edge_type=EdgeType(
                is_directed=is_directed,
                is_hybrid=is_hybrid,
                is_parallel=is_parallel,
            ),
            points=points,
        )

    return edge_routes
