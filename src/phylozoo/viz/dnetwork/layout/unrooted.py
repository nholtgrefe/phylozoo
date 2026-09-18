"""
Tree-of-blobs layout for DirectedPhyNetwork (pz-unrooted).

The directed network is drawn as its underlying undirected structure (the root
is an ordinary tree node of the tree of blobs); direction is shown by the
arrowheads. The placement lives in :mod:`phylozoo.viz._unrooted`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from phylozoo.viz._unrooted import simple_graph, tree_of_blobs_positions

from .base import DNetLayout
from .routes import compute_backbone_routes, compute_hybrid_routes

if TYPE_CHECKING:
    from phylozoo.core.network.dnetwork import DirectedPhyNetwork


def compute_pz_unrooted_layout(
    network: "DirectedPhyNetwork",
    edge_length: float = 1.0,
    node_spacing: float = 0.7,
    inner_scale: float = 0.6,
    start_angle: float = 0.0,
    refine: int = 50,
    daylight: int = 0,
) -> DNetLayout:
    """
    Compute a tree-of-blobs (unrooted-style) layout for a DirectedPhyNetwork.

    This is a custom PhyloZoo layout algorithm (pz-unrooted), identical to the
    semi-directed version (see
    :func:`phylozoo.viz.sdnetwork.layout.unrooted.compute_pz_unrooted_layout`)
    applied to the underlying undirected structure: the root is not
    distinguished and every edge carries an arrowhead when plotted.

    Parameters
    ----------
    network : DirectedPhyNetwork
        The network to layout.
    edge_length : float, optional
        Length of cut edges. By default 1.0.
    node_spacing : float, optional
        Target distance between consecutive nodes on a blob's circle. By default 0.7.
    inner_scale : float, optional
        Scale of subtrees hanging off interior blob nodes. By default 0.6.
    start_angle : float, optional
        Rotation of the whole drawing, in radians. By default 0.0.
    refine : int, optional
        Stress-majorization iterations (0 = raw geometry). By default 50.
    daylight : int, optional
        Equal-daylight rounds applied before the refinement. By default 0.

    Returns
    -------
    DNetLayout
        The computed layout, with ``algorithm='pz-unrooted'``.

    Examples
    --------
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> from phylozoo.viz.dnetwork.layout import compute_pz_unrooted_layout
    >>>
    >>> net = DirectedPhyNetwork(
    ...     edges=[(5, 3), (5, 4), (3, 2), (4, 2), (2, 1), (3, 6), (4, 7)],
    ...     nodes=[(1, {'label': 'A'}), (6, {'label': 'B'}), (7, {'label': 'C'})]
    ... )
    >>> layout = compute_pz_unrooted_layout(net)
    >>> layout.algorithm
    'pz-unrooted'
    >>> len(layout.positions) == net.number_of_nodes()
    True
    """
    positions = tree_of_blobs_positions(
        simple_graph(network),
        set(network.leaves),
        edge_length=edge_length,
        node_spacing=node_spacing,
        inner_scale=inner_scale,
        start_angle=start_angle,
        refine=refine,
        daylight=daylight,
    )
    tree_edges: set[tuple[Any, Any, int]] = set()
    hybrid_edges: set[tuple[Any, Any, int]] = set()
    for u, v, key in network._graph.edges(keys=True):
        (hybrid_edges if v in network.hybrid_nodes else tree_edges).add((u, v, key))
    routes = compute_backbone_routes(network, positions, tree_edges)
    routes.update(compute_hybrid_routes(network, positions, hybrid_edges))
    return DNetLayout(
        network=network,
        positions=positions,
        edge_routes=routes,
        backbone_edges=tree_edges,
        reticulate_edges=hybrid_edges,
        algorithm="pz-unrooted",
        parameters={
            "edge_length": edge_length,
            "node_spacing": node_spacing,
            "inner_scale": inner_scale,
            "start_angle": start_angle,
            "refine": refine,
            "daylight": daylight,
        },
    )
