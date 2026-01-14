"""
NetworkX and Graphviz layout computation for graphs.

This module provides layout computation using standard NetworkX and Graphviz
algorithms for DirectedMultiGraph and MixedMultiGraph.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, TypeVar

import networkx as nx

from phylozoo.utils.exceptions import (
    PhyloZooImportError,
    PhyloZooLayoutError,
    PhyloZooTypeError,
)

if TYPE_CHECKING:
    from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
    from phylozoo.core.primitives.m_multigraph import MixedMultiGraph

T = TypeVar('T')


@dataclass(frozen=True)
class GraphLayout:
    """
    Simple layout for graphs.

    This layout class stores computed positions for graph nodes using
    standard NetworkX or Graphviz layouts.

    Attributes
    ----------
    graph : DirectedMultiGraph | MixedMultiGraph
        The graph this layout is for.
    positions : dict[T, tuple[float, float]]
        Node positions mapping node ID to (x, y) coordinates.
    algorithm : str
        Name of the layout algorithm used.
    parameters : dict[str, Any]
        Parameters used to generate this layout.

    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
    >>> from phylozoo.viz.graphs.layout import compute_graph_layout
    >>> G = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
    >>> layout = compute_graph_layout(G, layout='spring')
    >>> len(layout.positions)
    3
    """

    graph: 'DirectedMultiGraph[T] | MixedMultiGraph[T]'
    positions: dict[T, tuple[float, float]]
    algorithm: str
    parameters: dict[str, Any]

    def get_position(self, node: T) -> tuple[float, float]:
        """
        Get position of a node.

        Parameters
        ----------
        node : T
            Node ID.

        Returns
        -------
        tuple[float, float]
            (x, y) position.

        Raises
        ------
        KeyError
            If node is not in the layout.
        """
        return self.positions[node]


def compute_graph_layout(
    graph: 'DirectedMultiGraph[T] | MixedMultiGraph[T]',
    layout: str = 'spring',
    **kwargs: Any,
) -> GraphLayout:
    """
    Compute layout for a graph using NetworkX or Graphviz.

    Parameters
    ----------
    graph : DirectedMultiGraph | MixedMultiGraph
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
    GraphLayout
        The computed layout.

    Raises
    ------
    PhyloZooLayoutError
        If layout algorithm is not supported, graph is empty, or layout computation fails.
    PhyloZooTypeError
        If graph type is unsupported.
    PhyloZooImportError
        If Graphviz layout is requested but pygraphviz is not installed.

    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
    >>> G = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
    >>> layout = compute_graph_layout(G, layout='circular')
    >>> len(layout.positions)
    3
    """
    if graph.number_of_nodes() == 0:
        raise PhyloZooLayoutError("Cannot compute layout for empty graph")

    # Convert graph to NetworkX format for layout
    if hasattr(graph, '_graph'):
        # DirectedMultiGraph
        G_nx = graph._graph
    elif hasattr(graph, '_directed') and hasattr(graph, '_undirected'):
        # MixedMultiGraph - combine both for layout
        G_nx = nx.MultiGraph()
        for node in graph.nodes():
            G_nx.add_node(node)
        for u, v, key in graph._directed.edges(keys=True):
            G_nx.add_edge(u, v, key=key)
        for u, v, key in graph._undirected.edges(keys=True):
            G_nx.add_edge(u, v, key=key)
    else:
        raise PhyloZooTypeError(f"Unsupported graph type: {type(graph)}")

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
                f"Graphviz layout '{layout}' failed: {e}. "
                "Falling back to spring layout."
            ) from e
            # Fallback to spring
            pos = nx.spring_layout(G_nx, **kwargs)
            layout = 'spring'
    
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

    return GraphLayout(
        graph=graph,
        positions=pos,
        algorithm=layout,
        parameters=kwargs,
    )

