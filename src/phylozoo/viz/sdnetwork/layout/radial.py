"""
Radial layout for SemiDirectedPhyNetwork (pz-radial).

The network is rooted (at a chosen or automatically selected location) and
drawn as a circular cladogram with the root at the centre and the leaves on
the outer circle, reticulate edges as straight chords; see
:func:`phylozoo.viz.dnetwork.layout.radial.compute_pz_radial_layout`.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any

from phylozoo.core.network.sdnetwork.derivations import to_d_network
from phylozoo.utils.exceptions import (
    PhyloZooEmptyNetworkWarning,
    PhyloZooLayoutError,
    PhyloZooSingleNodeNetworkWarning,
)
from phylozoo.viz.dnetwork.layout.radial import compute_pz_radial_layout as _d_radial

from .base import SDNetLayout
from .routes import compute_radial_routes

if TYPE_CHECKING:
    from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork


def compute_pz_radial_layout(
    network: "SemiDirectedPhyNetwork",
    radius: float = 1.0,
    start_angle: float = 0.0,
    angle_direction: str = "clockwise",
    root_location: Any = None,
    **kwargs: Any,
) -> SDNetLayout:
    """
    Compute a radial (circular cladogram) layout for a SemiDirectedPhyNetwork.

    This is a custom PhyloZoo layout algorithm (pz-radial). The network is
    rooted with :func:`~phylozoo.core.network.sdnetwork.derivations.to_d_network`
    and the rooted network is drawn as a circular cladogram: root at the
    centre, leaves evenly spaced on the outer circle, radius growing with the
    depth, reticulate edges as straight chords. The leaf order comes from the
    hybrid-aware ordering of the layered ``pz-cladogram`` layout.

    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The network to layout (any rootable network, not only trees).
    radius : float, optional
        Radius of the leaf circle. By default 1.0.
    start_angle : float, optional
        Angle of the first leaf, in radians. By default 0.0.
    angle_direction : str, optional
        'clockwise' or 'counterclockwise'. By default 'clockwise'.
    root_location : node or edge, optional
        Where to root the network (a node, or an edge ``(u, v, key)``), as
        accepted by ``to_d_network``. None lets ``to_d_network`` choose.
        By default None.
    **kwargs
        Ordering options passed on to
        :func:`~phylozoo.viz.dnetwork.layout.cladogram.compute_pz_cladogram_layout`
        (``trials``, ``seed``).

    Returns
    -------
    SDNetLayout
        The computed layout.

    Raises
    ------
    PhyloZooLayoutError
        If network is empty.
    PhyloZooValueError
        If angle_direction is invalid.

    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> from phylozoo.viz.sdnetwork.layout import compute_pz_radial_layout
    >>>
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2), (3, 100)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'}), (100, {'label': 'C'})]
    ... )
    >>> layout = compute_pz_radial_layout(net)
    >>> len(layout.positions)
    4
    """
    if network.number_of_nodes() == 0:
        raise PhyloZooLayoutError("Cannot compute layout for empty network")
    if network.number_of_nodes() == 1:
        node = next(iter(network._graph.nodes))
        return SDNetLayout(
            network=network,
            positions={node: (0.0, 0.0)},
            edge_routes={},
            algorithm="pz-radial",
            parameters={"radius": radius, "start_angle": start_angle},
        )

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=PhyloZooEmptyNetworkWarning)
        warnings.filterwarnings("ignore", category=PhyloZooSingleNodeNetworkWarning)
        rooted = to_d_network(network, root_location=root_location)
    d_layout = _d_radial(
        rooted, radius=radius, start_angle=start_angle, angle_direction=angle_direction, **kwargs
    )
    # Rooting on an edge adds a root node that the semi-directed network does not have.
    original = set(network._graph.nodes)
    positions = {n: p for n, p in d_layout.positions.items() if n in original}
    return SDNetLayout(
        network=network,
        positions=positions,
        edge_routes=compute_radial_routes(network, positions),
        algorithm="pz-radial",
        parameters={
            "radius": radius,
            "start_angle": start_angle,
            "angle_direction": angle_direction,
            "root_location": root_location,
            **kwargs,
        },
    )
