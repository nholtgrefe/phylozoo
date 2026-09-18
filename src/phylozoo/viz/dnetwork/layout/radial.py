"""
Radial layout for DirectedPhyNetwork (pz-radial).

A circular cladogram: the root sits at the centre, every leaf on the outer
circle, and a node's radius grows with its layer. The angular order of the
leaves and the tree backbone come from the layered ``pz-cladogram`` computation, so
reticulate edges are kept short and drawn as straight chords, as in the radial
network view of Dendroscope.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Any

import numpy as np

from phylozoo.utils.exceptions import PhyloZooValueError
from phylozoo.viz._layout_utils import count_crossings

from .base import DNetLayout
from .cladogram import _local_search, _preorder, compute_pz_cladogram_layout
from .routes import compute_backbone_routes, compute_hybrid_routes

if TYPE_CHECKING:
    from phylozoo.core.network.dnetwork import DirectedPhyNetwork


def compute_pz_radial_layout(
    network: "DirectedPhyNetwork",
    radius: float = 1.0,
    start_angle: float = 0.0,
    angle_direction: str = "clockwise",
    **kwargs: Any,
) -> DNetLayout:
    """
    Compute a radial (circular cladogram) layout for a DirectedPhyNetwork.

    This is a custom PhyloZoo layout algorithm (pz-radial). The layered
    ``pz-cladogram`` layout is computed first (tree backbone, hybrid-aware leaf
    order, leaves on the bottom layer) and then mapped to polar coordinates:
    the leaf order becomes the angular order (leaves evenly spaced on the
    outer circle), a node's layer becomes its radius, and every internal node
    sits at the mean angle of its backbone children. The sibling order is then
    re-optimised by adjacent swaps with the crossings counted in the circle,
    so reticulate edges (straight chords) cross as little as possible there.

    Parameters
    ----------
    network : DirectedPhyNetwork
        The network to layout.
    radius : float, optional
        Radius of the leaf circle. By default 1.0.
    start_angle : float, optional
        Angle of the first leaf, in radians. By default 0.0.
    angle_direction : str, optional
        'clockwise' or 'counterclockwise'. By default 'clockwise'.
    **kwargs
        Ordering options of :func:`~phylozoo.viz.dnetwork.layout.cladogram.compute_pz_cladogram_layout`
        (``trials``, ``seed``, ``horizontal_reticulations``).

    Returns
    -------
    DNetLayout
        The computed layout, with ``algorithm='pz-radial'``.

    Raises
    ------
    PhyloZooLayoutError
        If network is empty.
    PhyloZooValueError
        If angle_direction is invalid.

    Examples
    --------
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> from phylozoo.viz.dnetwork.layout import compute_pz_radial_layout
    >>>
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> layout = compute_pz_radial_layout(net)
    >>> layout.get_position(3)
    (0.0, 0.0)
    >>> round(math.hypot(*layout.get_position(1)), 6)
    1.0
    """
    if angle_direction not in ("clockwise", "counterclockwise"):
        raise PhyloZooValueError(
            f"angle_direction must be 'clockwise' or 'counterclockwise', got '{angle_direction}'"
        )
    layered = compute_pz_cladogram_layout(network, direction="TD", align_leaves=True, **kwargs)
    kids: dict[Any, list[Any]] = {n: [] for n in layered.positions}
    for u, v, _ in layered.backbone_edges:
        kids[u].append(v)
    for cs in kids.values():  # start from the layered left-to-right order
        cs.sort(key=lambda n: layered.positions[n][0])
    root = network.root_node
    ys = [y for _, y in layered.positions.values()]
    top, bottom = max(ys), min(ys)
    radius_of = {
        n: radius * (top - y) / (top - bottom) if top > bottom else 0.0
        for n, (_, y) in layered.positions.items()
    }
    sign = -1.0 if angle_direction == "clockwise" else 1.0
    edges = [(u, v) for u, v, _ in layered.backbone_edges | layered.reticulate_edges]
    reticulate = [(u, v) for u, v, _ in layered.reticulate_edges]

    def polar_positions(order: dict[Any, list[Any]]) -> dict[Any, tuple[float, float]]:
        """Leaves evenly around the circle in backbone order, internal nodes at the mean angle."""
        leaves = [n for n in _preorder(root, order) if not order[n]]
        step = 2 * math.pi / len(leaves)
        angle = {leaf: start_angle + sign * i * step for i, leaf in enumerate(leaves)}

        def mean_angle(n: Any) -> float:
            if n not in angle:
                angle[n] = sum(mean_angle(c) for c in order[n]) / len(order[n])
            return angle[n]

        return {
            n: (
                radius_of[n] * math.cos(mean_angle(n)) + 0.0,
                radius_of[n] * math.sin(mean_angle(n)) + 0.0,
            )
            for n in order
        }

    def polar_score(order: dict[Any, list[Any]]) -> tuple[int, float]:
        pos = polar_positions(order)
        segments = np.array([(*pos[u], *pos[v]) for u, v in edges], dtype=float)
        chords = sum(math.dist(pos[u], pos[v]) for u, v in reticulate)
        return count_crossings(segments), round(chords, 9)

    # Re-run the sibling swap search with crossings measured in the circle, on
    # the nodes that can affect them (ancestors of reticulate-edge endpoints).
    parent_of = {c: p for p, cs in kids.items() for c in cs}
    relevant: set[Any] = set()
    for u, v in reticulate:
        for n in (u, v):
            while n in parent_of:
                n = parent_of[n]
                relevant.add(n)
    depth = {n: 0 for n in kids}
    _local_search(
        root,
        kids,
        edges,
        reticulate,
        depth,
        1.0,
        score=polar_score,
        nodes=[n for n in kids if n in relevant and len(kids[n]) > 1],
    )
    positions = polar_positions(kids)

    routes = compute_backbone_routes(network, positions, layered.backbone_edges)
    routes.update(compute_hybrid_routes(network, positions, layered.reticulate_edges))
    return DNetLayout(
        network=network,
        positions=positions,
        edge_routes=routes,
        backbone_edges=layered.backbone_edges,
        reticulate_edges=layered.reticulate_edges,
        algorithm="pz-radial",
        parameters={
            "radius": radius,
            "start_angle": start_angle,
            "angle_direction": angle_direction,
            **kwargs,
        },
    )
