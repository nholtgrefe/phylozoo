"""
NetworkX and Graphviz layout dispatcher for MixedMultiGraph.

This module provides layout computation using standard NetworkX and Graphviz
algorithms for MixedMultiGraph.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, TypeVar

import networkx as nx

from phylozoo.utils.exceptions import PhyloZooLayoutError

from phylozoo.viz._layout_utils import compute_nx_positions, normalize_positions

from .base import MGraphLayout
from .routes import compute_mmgraph_routes

if TYPE_CHECKING:
    from phylozoo.core.primitives.m_multigraph import MixedMultiGraph

T = TypeVar('T')


def compute_nx_layout(
    graph: 'MixedMultiGraph[T]',
    layout: str = 'spring',
    **kwargs: Any,
) -> MGraphLayout[T]:
    """
    Compute layout for MixedMultiGraph using NetworkX or Graphviz.

    Parameters
    ----------
    graph : MixedMultiGraph
        The graph to layout.
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
    MGraphLayout
        The computed layout.

    Raises
    ------
    PhyloZooLayoutError
        If layout algorithm is not supported, graph is empty, or layout computation fails.
    PhyloZooImportError
        If Graphviz layout is requested but pygraphviz is not installed.

    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
    >>> from phylozoo.viz.m_multigraph.layout.nx import compute_nx_layout
    >>>
    >>> G = MixedMultiGraph(directed_edges=[(1, 2)], undirected_edges=[(2, 3)])
    >>> layout = compute_nx_layout(G, layout='circular')
    >>> len(layout.positions)
    3
    """
    if graph.number_of_nodes() == 0:
        raise PhyloZooLayoutError("Cannot compute layout for empty graph")

    G_nx = nx.MultiGraph()
    for node in graph.nodes():
        G_nx.add_node(node)
    for u, v, key in graph._directed.edges(keys=True):
        G_nx.add_edge(u, v, key=key)
    for u, v, key in graph._undirected.edges(keys=True):
        G_nx.add_edge(u, v, key=key)

    pos = compute_nx_positions(G_nx, layout=layout, **kwargs)
    pos = normalize_positions(pos)

    # Compute edge routes
    edge_routes = compute_mmgraph_routes(graph, pos)

    return MGraphLayout(
        network=graph,
        positions=pos,
        edge_routes=edge_routes,
        algorithm=layout,
        parameters=kwargs,
    )
