"""
NetworkX and Graphviz layout dispatcher for DirectedPhyNetwork.

This module provides layout computation using standard NetworkX and Graphviz
algorithms for DirectedPhyNetwork.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, TypeVar

import networkx as nx

from phylozoo.utils.exceptions import PhyloZooLayoutError

from phylozoo.viz._layout_utils import compute_nx_positions, normalize_positions

from .base import DNetLayout
from .routes import compute_backbone_routes, compute_hybrid_routes

if TYPE_CHECKING:
    from phylozoo.core.network.dnetwork import DirectedPhyNetwork

T = TypeVar('T')


def compute_nx_layout(
    network: 'DirectedPhyNetwork',
    layout: str = 'spring',
    **kwargs: Any,
) -> DNetLayout:
    """
    Compute layout for DirectedPhyNetwork using NetworkX or Graphviz.

    Parameters
    ----------
    network : DirectedPhyNetwork
        The network to layout.
    layout : str, optional
        Layout algorithm name. Supported values:
        - NetworkX: 'spring', 'circular', 'kamada_kawai', 'planar', 'random',
          'shell', 'spectral', 'spiral', 'bipartite'
        - Graphviz: 'dot', 'neato', 'fdp', 'sfdp', 'twopi', 'circo'
        By default 'spring'.
    **kwargs
        Additional parameters passed to the layout algorithm.

    Returns
    -------
    DNetLayout
        The computed layout.

    Raises
    ------
    PhyloZooLayoutError
        If layout algorithm is not supported, network is empty, or layout computation fails.
    PhyloZooImportError
        If Graphviz layout is requested but pygraphviz is not installed.

    Examples
    --------
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> from phylozoo.viz.dnetwork.layout.nx import compute_nx_layout
    >>>
    >>> net = DirectedPhyNetwork(edges=[(3, 1), (3, 2)])
    >>> layout = compute_nx_layout(net, layout='circular')
    >>> len(layout.positions)
    3
    """
    if network.number_of_nodes() == 0:
        raise PhyloZooLayoutError("Cannot compute layout for empty network")

    G_nx = nx.DiGraph()
    for node in network._graph.nodes:
        G_nx.add_node(node)
    for u, v, key in network._graph.edges(keys=True):
        G_nx.add_edge(u, v)

    pos = compute_nx_positions(G_nx, layout=layout, **kwargs)
    pos = normalize_positions(pos)

    # Compute edge routes
    # All edges are directed in DirectedPhyNetwork
    all_edges = set(network._graph.edges(keys=True))
    hybrid_edges = set(
        (u, v, key)
        for u, v, key in network._graph.edges(keys=True)
        if v in network.hybrid_nodes
    )
    backbone_edges = all_edges - hybrid_edges

    edge_routes = {}
    if backbone_edges:
        backbone_routes = compute_backbone_routes(network, pos, backbone_edges)
        edge_routes.update(backbone_routes)
    if hybrid_edges:
        hybrid_routes = compute_hybrid_routes(network, pos, hybrid_edges)
        edge_routes.update(hybrid_routes)

    return DNetLayout(
        network=network,
        positions=pos,
        edge_routes=edge_routes,
        backbone_edges=backbone_edges,
        reticulate_edges=hybrid_edges,
        algorithm=layout,
        parameters=kwargs,
    )
