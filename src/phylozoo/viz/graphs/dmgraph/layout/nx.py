"""
NetworkX and Graphviz layout dispatcher for DirectedMultiGraph.

This module provides layout computation using standard NetworkX and Graphviz
algorithms for DirectedMultiGraph.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, TypeVar

from phylozoo.utils.exceptions import PhyloZooLayoutError

from phylozoo.viz._layout_utils import compute_nx_positions, normalize_positions

from .base import DMGraphLayout
from .routes import compute_dmgraph_routes

if TYPE_CHECKING:
    from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph

T = TypeVar('T')


def compute_nx_layout(
    graph: 'DirectedMultiGraph[T]',
    layout: str = 'spring',
    **kwargs: Any,
) -> DMGraphLayout[T]:
    """
    Compute layout for DirectedMultiGraph using NetworkX or Graphviz.

    Parameters
    ----------
    graph : DirectedMultiGraph
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
    DMGraphLayout
        The computed layout.

    Raises
    ------
    PhyloZooLayoutError
        If layout algorithm is not supported, graph is empty, or layout computation fails.
    PhyloZooImportError
        If Graphviz layout is requested but pygraphviz is not installed.

    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
    >>> from phylozoo.viz.graphs.dmgraph.layout.nx import compute_nx_layout
    >>>
    >>> G = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
    >>> layout = compute_nx_layout(G, layout='circular')
    >>> len(layout.positions)
    3
    """
    if graph.number_of_nodes() == 0:
        raise PhyloZooLayoutError("Cannot compute layout for empty graph")

    G_nx = graph._graph
    pos = compute_nx_positions(G_nx, layout=layout, **kwargs)
    pos = normalize_positions(pos)

    # Compute edge routes
    edge_routes = compute_dmgraph_routes(graph, pos)

    return DMGraphLayout(
        network=graph,
        positions=pos,
        edge_routes=edge_routes,
        algorithm=layout,
        parameters=kwargs,
    )
