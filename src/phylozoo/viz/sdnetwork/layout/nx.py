"""
NetworkX and Graphviz layout dispatcher for SemiDirectedPhyNetwork.

This module provides layout computation using standard NetworkX and Graphviz
algorithms for SemiDirectedPhyNetwork.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, TypeVar

import networkx as nx

from phylozoo.utils.exceptions import PhyloZooLayoutError

from phylozoo.viz._layout_utils import compute_nx_positions, normalize_positions

from .base import SDNetLayout
from .routes import compute_radial_routes

if TYPE_CHECKING:
    from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork

T = TypeVar('T')


def compute_nx_layout(
    network: 'SemiDirectedPhyNetwork',
    layout: str = 'spring',
    **kwargs: Any,
) -> SDNetLayout:
    """
    Compute layout for SemiDirectedPhyNetwork using NetworkX or Graphviz.

    Parameters
    ----------
    network : SemiDirectedPhyNetwork
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
    SDNetLayout
        The computed layout.

    Raises
    ------
    PhyloZooLayoutError
        If layout algorithm is not supported, network is empty, or layout computation fails.
    PhyloZooImportError
        If Graphviz layout is requested but pygraphviz is not installed.

    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> from phylozoo.viz.sdnetwork.layout.nx import compute_nx_layout
    >>>
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> layout = compute_nx_layout(net, layout='circular')
    >>> len(layout.positions)
    3
    """
    if network.number_of_nodes() == 0:
        raise PhyloZooLayoutError("Cannot compute layout for empty network")

    G_nx = nx.MultiGraph()
    original_nodes = set(network._graph.nodes)
    for node in original_nodes:
        G_nx.add_node(node)
    for u, v, key in network._graph.edges(keys=True):
        if u in original_nodes and v in original_nodes:
            G_nx.add_edge(u, v, key=key)

    pos = compute_nx_positions(G_nx, layout=layout, **kwargs)
    pos = normalize_positions(pos)

    # Filter positions to only include original network nodes
    filtered_positions: dict[T, tuple[float, float]] = {
        node: pos for node, pos in pos.items() if node in original_nodes
    }

    # Compute edge routes
    edge_routes = compute_radial_routes(network, filtered_positions)

    return SDNetLayout(
        network=network,
        positions=filtered_positions,
        edge_routes=edge_routes,
        algorithm=layout,
        parameters=kwargs,
    )
