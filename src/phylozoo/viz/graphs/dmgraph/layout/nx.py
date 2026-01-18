"""
NetworkX and Graphviz layout dispatcher for DirectedMultiGraph.

This module provides layout computation using standard NetworkX and Graphviz
algorithms for DirectedMultiGraph.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, TypeVar

import networkx as nx

from phylozoo.utils.exceptions import (
    PhyloZooImportError,
    PhyloZooLayoutError,
)

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

    # Convert graph to NetworkX DiGraph for layout
    G_nx = graph._graph

    # Compute layout
    pos: dict[T, tuple[float, float]]

    # Graphviz layouts
    if layout in ('dot', 'neato', 'fdp', 'sfdp', 'twopi', 'circo'):
        try:
            pos = nx.nx_agraph.graphviz_layout(G_nx, prog=layout, **kwargs)
        except ImportError:
            raise PhyloZooImportError(
                f"Graphviz layout '{layout}' requires pygraphviz. "
                "Install with: pip install pygraphviz"
            )
        except Exception as e:
            raise PhyloZooLayoutError(
                f"Graphviz layout '{layout}' failed: {e}"
            ) from e

    # NetworkX layouts
    elif layout == 'spring':
        pos = nx.spring_layout(G_nx, **kwargs)
    elif layout == 'circular':
        pos = nx.circular_layout(G_nx, **kwargs)
    elif layout == 'kamada_kawai':
        pos = nx.kamada_kawai_layout(G_nx, **kwargs)
    elif layout == 'planar':
        pos = nx.planar_layout(G_nx, **kwargs)
    elif layout == 'random':
        pos = nx.random_layout(G_nx, **kwargs)
    elif layout == 'shell':
        pos = nx.shell_layout(G_nx, **kwargs)
    elif layout == 'spectral':
        pos = nx.spectral_layout(G_nx, **kwargs)
    elif layout == 'spiral':
        pos = nx.spiral_layout(G_nx, **kwargs)
    elif layout == 'bipartite':
        pos = nx.bipartite_layout(G_nx, **kwargs)
    else:
        raise PhyloZooLayoutError(
            f"Unsupported layout algorithm: '{layout}'. "
            "Supported: spring, circular, kamada_kawai, planar, random, "
            "shell, spectral, spiral, bipartite, dot, neato, fdp, sfdp, twopi, circo"
        )

    # Normalize positions to center at origin and scale appropriately
    if pos:
        xs = [x for x, _ in pos.values()]
        ys = [y for _, y in pos.values()]
        if xs and ys:
            min_x, max_x = min(xs), max(xs)
            min_y, max_y = min(ys), max(ys)
            width = max_x - min_x if max_x != min_x else 1.0
            height = max_y - min_y if max_y != min_y else 1.0

            # Center and normalize
            center_x = (min_x + max_x) / 2
            center_y = (min_y + max_y) / 2
            scale = 1.0 / max(width, height) if max(width, height) > 0 else 1.0

            normalized_pos: dict[T, tuple[float, float]] = {}
            for node, (x, y) in pos.items():
                normalized_pos[node] = (
                    (x - center_x) * scale,
                    (y - center_y) * scale,
                )
            pos = normalized_pos

    # Compute edge routes
    edge_routes = compute_dmgraph_routes(graph, pos)

    return DMGraphLayout(
        network=graph,
        positions=pos,
        edge_routes=edge_routes,
        algorithm=layout,
        parameters=kwargs,
    )
