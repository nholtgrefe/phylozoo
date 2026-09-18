"""
Tree-of-blobs layout for SemiDirectedPhyNetwork (pz-unrooted).

The placement itself lives in :mod:`phylozoo.viz._unrooted` and is shared with
the directed version.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from phylozoo.viz._unrooted import simple_graph, tree_of_blobs_positions

from .base import SDNetLayout
from .routes import compute_radial_routes

if TYPE_CHECKING:
    from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork


def compute_pz_unrooted_layout(
    network: "SemiDirectedPhyNetwork",
    edge_length: float = 1.0,
    node_spacing: float = 0.7,
    inner_scale: float = 0.6,
    start_angle: float = 0.0,
    refine: int = 50,
    daylight: int = 0,
) -> SDNetLayout:
    """
    Compute a tree-of-blobs layout for a SemiDirectedPhyNetwork.

    This is a custom PhyloZoo layout algorithm (pz-unrooted). It works for any
    connected semi-directed network; for trees it reduces to the classic
    equal-angle drawing of an unrooted tree.

    The algorithm:

    1. Decomposes the network into blobs and builds the tree of blobs
       (one node per blob, one edge per cut edge).
    2. Roots that tree at its centroid and assigns every subtree an angular
       wedge proportional to the number of leaves it contains (equal-angle).
    3. Draws each non-trivial blob as a polygon: the nodes on the outer face of
       a planar embedding lie evenly spaced on a circle, in face order, and the
       remaining nodes are placed at the barycentre of their neighbours (Tutte
       embedding).
       Subtrees attached to interior nodes are drawn inside the blob at a
       reduced scale.
    4. Optionally equalises the "daylight" around every node (``daylight``
       rounds, blobs kept rigid) and relaxes the drawing with iterations of stress
       majorization (the algorithm behind Graphviz ``neato``), started from the
       crossing-free drawing of step 3, with every blob edge targeting
       ``node_spacing``. This spreads cramped subtrees and keeps blobs evenly
       sided, and is deterministic.

    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The network to layout.
    edge_length : float, optional
        Length of cut edges (tree edges between blobs and pendant edges).
        By default 1.0.
    node_spacing : float, optional
        Target distance between consecutive nodes on a blob's circle.
        By default 0.7.
    inner_scale : float, optional
        Scale factor for subtrees hanging off interior blob nodes, which are
        drawn inside the blob. By default 0.6.
    start_angle : float, optional
        Rotation of the whole drawing, in radians. By default 0.0.
    refine : int, optional
        Number of stress-majorization iterations applied to the tree-of-blobs
        drawing. Use 0 for the raw geometric drawing (blobs as exact polygons).
        By default 50.
    daylight : int, optional
        Number of equal-daylight rounds applied before the stress refinement:
        around every tree node and around every blob (as a whole), the
        subtrees hanging off it through cut edges are rotated so that the empty
        angular gaps between them become equal (blobs stay rigid). Spreads the
        drawing while keeping every edge length exact (daylight=8 with refine=0 gives
        the classic SplitsTree-style unrooted drawing).
        By default 0.

    Returns
    -------
    SDNetLayout
        The computed layout. Positions are centred at the origin and scaled to
        unit range.

    Raises
    ------
    PhyloZooLayoutError
        If the network is empty or not connected.

    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> from phylozoo.viz.sdnetwork.layout import compute_pz_unrooted_layout
    >>>
    >>> net = SemiDirectedPhyNetwork(
    ...     directed_edges=[(5, 4), (6, 4)],
    ...     undirected_edges=[(5, 3), (5, 6), (6, 7), (4, 8), (8, 1), (8, 2)],
    ...     nodes=[(3, {'label': 'C'}), (7, {'label': 'D'}),
    ...            (1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> layout = compute_pz_unrooted_layout(net)
    >>> len(layout.positions) == net.number_of_nodes()
    True
    >>> layout.algorithm
    'pz-unrooted'
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
    return SDNetLayout(
        network=network,
        positions=positions,
        edge_routes=compute_radial_routes(network, positions),
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
