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

T = TypeVar("T")


def compute_backbone_routes(
    network: "DirectedPhyNetwork",
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
            # Edges into hybrid nodes are hybrid edges (one may be in backbone)
            is_hybrid = v in network.hybrid_nodes

            points = (positions[u], positions[v])

            edge_routes[(u, v, key)] = EdgeRoute(
                edge_type=EdgeType(
                    is_directed=True,
                    is_hybrid=is_hybrid,
                    is_parallel=is_parallel,
                ),
                points=points,
            )

    return edge_routes


def compute_hybrid_routes(
    network: "DirectedPhyNetwork",
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
    >>> net = DirectedPhyNetwork(
    ...     edges=[(5, 3), (5, 4), (3, 2), (4, 2), (2, 1), (3, 6), (4, 7)],
    ...     nodes=[(1, {'label': 'A'}), (6, {'label': 'B'}), (7, {'label': 'C'})]
    ... )
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
        )

    return edge_routes


def compute_rectangular_routes(
    network: "DirectedPhyNetwork",
    positions: dict[T, tuple[float, float]],
    backbone_edges: set[tuple[T, T, int]],
    reticulate_edges: set[tuple[T, T, int]],
    direction: str = "TD",
) -> dict[tuple[T, T, int], EdgeRoute]:
    """
    Route edges orthogonally (elbows), for ``pz-cladogram`` with ``rectangular=True``.

    A backbone edge runs along the parent's layer to the child's position and
    then along the layer axis into the child (the classic cladogram elbow).

    A reticulate edge enters the hybrid node sideways, along the hybrid's
    layer, so it never overlaps the backbone edge arriving from above. If the
    parent sits on the hybrid's layer (see ``horizontal_reticulations`` in
    :func:`~phylozoo.viz.dnetwork.layout.cladogram.compute_pz_cladogram_layout`) the edge
    is a single horizontal segment; otherwise it first runs along the layer
    axis at the parent's position, shifted slightly sideways if that column
    already carries a backbone edge.

    Parameters
    ----------
    network : DirectedPhyNetwork
        The network.
    positions : dict[T, tuple[float, float]]
        Node positions.
    backbone_edges : set[tuple[T, T, int]]
        Backbone tree edges.
    reticulate_edges : set[tuple[T, T, int]]
        Remaining (reticulate) edges.
    direction : str, optional
        'TD' (layers are horizontal) or 'LR' (layers are vertical). By default 'TD'.

    Returns
    -------
    dict[tuple[T, T, int], EdgeRoute]
        Edge routing information; backbone routes have three points, reticulate
        routes two to four.

    Examples
    --------
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)])
    >>> positions = {3: (0.5, 1.0), 1: (0.0, 0.0), 2: (1.0, 0.0)}
    >>> routes = compute_rectangular_routes(net, positions, {(3, 1, 0), (3, 2, 0)}, set())
    >>> routes[(3, 1, 0)].points
    ((0.5, 1.0), (0.0, 1.0), (0.0, 0.0))
    """
    top_down = direction.upper() == "TD"

    # Work in (across, depth) coordinates: depth grows away from the root.
    def split(point: tuple[float, float]) -> tuple[float, float]:
        return (point[0], -point[1]) if top_down else (point[1], point[0])

    def join(across: float, depth: float) -> tuple[float, float]:
        return (across, -depth) if top_down else (depth, across)

    kids: dict[T, list[T]] = {}
    for u, v, _ in backbone_edges:
        kids.setdefault(u, []).append(v)
    leaves = network.leaves
    columns = sorted({split(pt)[0] for v, pt in positions.items() if v in leaves})
    unit = min((b - a for a, b in zip(columns, columns[1:]) if b - a > 1e-9), default=1.0)

    parallel_counts: dict[tuple[T, T], int] = {}
    for u, v, _ in backbone_edges | reticulate_edges:
        parallel_counts[(u, v)] = parallel_counts.get((u, v), 0) + 1

    def route(u: T, v: T, key: int, points: tuple[tuple[float, float], ...]) -> None:
        edge_routes[(u, v, key)] = EdgeRoute(
            edge_type=EdgeType(
                is_directed=True,
                is_hybrid=v in network.hybrid_nodes,
                is_parallel=parallel_counts[(u, v)] > 1,
            ),
            points=points,
        )

    edge_routes: dict[tuple[T, T, int], EdgeRoute] = {}
    for u, v, key in backbone_edges:
        if u in positions and v in positions:
            (ua, ud), (va, vd) = split(positions[u]), split(positions[v])
            route(u, v, key, (join(ua, ud), join(va, ud), join(va, vd)))

    for p, h, key in reticulate_edges:
        if p not in positions or h not in positions:
            continue
        (pa, pd), (ha, hd) = split(positions[p]), split(positions[h])
        children = [split(positions[c]) for c in kids.get(p, []) if c in positions]
        if abs(pd - hd) < 1e-9:
            # Parent on the hybrid's layer: a single horizontal segment.
            route(p, h, key, (join(pa, pd), join(ha, hd)))
        elif any(abs(ca - pa) < 1e-9 for ca, _ in children):
            # The parent's column is taken by a backbone edge: step aside first.
            da = pa + 0.35 * unit * (1.0 if ha >= pa else -1.0)
            route(p, h, key, (join(pa, pd), join(da, pd), join(da, hd), join(ha, hd)))
        else:
            route(p, h, key, (join(pa, pd), join(pa, hd), join(ha, hd)))
    return edge_routes
