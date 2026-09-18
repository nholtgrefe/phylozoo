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

from phylozoo.utils.exceptions import PhyloZooValueError

from .base import DNetLayout
from .cladogram import compute_pz_cladogram_layout
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
    sits at the mean angle of its backbone children. Reticulate edges are
    straight chords.

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

    # Leaves evenly around the circle in the layered left-to-right order.
    leaves = sorted((n for n in kids if not kids[n]), key=lambda n: layered.positions[n][0])
    sign = -1.0 if angle_direction == "clockwise" else 1.0
    step = 2 * math.pi / len(leaves)
    angle = {leaf: start_angle + sign * i * step for i, leaf in enumerate(leaves)}

    def mean_angle(n: Any) -> float:
        if n not in angle:
            angle[n] = sum(mean_angle(c) for c in kids[n]) / len(kids[n])
        return angle[n]

    ys = [y for _, y in layered.positions.values()]
    top, bottom = max(ys), min(ys)
    positions: dict[Any, tuple[float, float]] = {}
    for n, (_, y) in layered.positions.items():
        r = radius * (top - y) / (top - bottom) if top > bottom else 0.0
        positions[n] = (r * math.cos(mean_angle(n)) + 0.0, r * math.sin(mean_angle(n)) + 0.0)

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
