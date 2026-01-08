"""
Public API for graph plotting.

This module provides the main plotting functions for DirectedMultiGraph
and MixedMultiGraph.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
    from phylozoo.core.primitives.m_multigraph import MixedMultiGraph

from ..dnetwork.backends import get_backend
from ..dnetwork.styling import NetworkStyle, default_style
from .layout import compute_graph_layout
from .rendering import compute_graph_routes


def plot_directed_multigraph(
    graph: 'DirectedMultiGraph',
    layout: str = 'spring',
    style: NetworkStyle | None = None,
    backend: str = 'matplotlib',
    ax: Any | None = None,
    **layout_kwargs: Any,
) -> Any:
    """
    Plot a DirectedMultiGraph.

    This function plots a directed multigraph using standard NetworkX or
    Graphviz layouts, with support for parallel edges displayed as curves.

    Parameters
    ----------
    graph : DirectedMultiGraph
        The graph to plot.
    layout : str, optional
        Layout algorithm. Supported: 'spring', 'circular', 'kamada_kawai',
        'planar', 'random', 'shell', 'spectral', 'spiral', 'bipartite',
        'dot', 'neato', 'fdp', 'sfdp', 'twopi', 'circo'.
        By default 'spring'.
    style : NetworkStyle, optional
        Styling configuration. If None, uses default style.
        By default None.
    backend : str, optional
        Backend to use ('matplotlib' or 'pyqtgraph').
        By default 'matplotlib'.
    ax : Any, optional
        Existing axes/figure to plot on (backend-specific).
        By default None.
    **layout_kwargs
        Additional parameters for layout computation (e.g., k, iterations, seed).

    Returns
    -------
    Any
        Backend-specific return type (e.g., matplotlib Axes or PyQtGraph window).

    Raises
    ------
    ValueError
        If layout algorithm is not supported or backend is not registered.

    Examples
    --------
    >>> from phylozoo.core.primitives.d_multigraph import DirectedMultiGraph
    >>> from phylozoo.visualize.graphs import plot_directed_multigraph
    >>> G = DirectedMultiGraph(edges=[(1, 2), (2, 3)])
    >>> ax = plot_directed_multigraph(G, layout='circular')
    """
    return _plot_graph(graph, layout, style, backend, ax, **layout_kwargs)


def plot_mixed_multigraph(
    graph: 'MixedMultiGraph',
    layout: str = 'spring',
    style: NetworkStyle | None = None,
    backend: str = 'matplotlib',
    ax: Any | None = None,
    **layout_kwargs: Any,
) -> Any:
    """
    Plot a MixedMultiGraph.

    This function plots a mixed multigraph (with both directed and undirected
    edges) using standard NetworkX or Graphviz layouts, with support for
    parallel edges displayed as curves.

    Parameters
    ----------
    graph : MixedMultiGraph
        The graph to plot.
    layout : str, optional
        Layout algorithm. Supported: 'spring', 'circular', 'kamada_kawai',
        'planar', 'random', 'shell', 'spectral', 'spiral', 'bipartite',
        'dot', 'neato', 'fdp', 'sfdp', 'twopi', 'circo'.
        By default 'spring'.
    style : NetworkStyle, optional
        Styling configuration. If None, uses default style.
        By default None.
    backend : str, optional
        Backend to use ('matplotlib' or 'pyqtgraph').
        By default 'matplotlib'.
    ax : Any, optional
        Existing axes/figure to plot on (backend-specific).
        By default None.
    **layout_kwargs
        Additional parameters for layout computation (e.g., k, iterations, seed).

    Returns
    -------
    Any
        Backend-specific return type (e.g., matplotlib Axes or PyQtGraph window).

    Raises
    ------
    ValueError
        If layout algorithm is not supported or backend is not registered.

    Examples
    --------
    >>> from phylozoo.core.primitives.m_multigraph import MixedMultiGraph
    >>> from phylozoo.visualize.graphs import plot_mixed_multigraph
    >>> G = MixedMultiGraph(directed_edges=[(1, 2)], undirected_edges=[(2, 3)])
    >>> ax = plot_mixed_multigraph(G, layout='spring')
    """
    return _plot_graph(graph, layout, style, backend, ax, **layout_kwargs)


def _plot_graph(
    graph: 'DirectedMultiGraph | MixedMultiGraph',
    layout: str,
    style: NetworkStyle | None,
    backend: str,
    ax: Any | None,
    **layout_kwargs: Any,
) -> Any:
    """
    Internal function to plot a graph.

    Parameters
    ----------
    graph : DirectedMultiGraph | MixedMultiGraph
        The graph to plot.
    layout : str
        Layout algorithm name.
    style : NetworkStyle | None
        Styling configuration.
    backend : str
        Backend name.
    ax : Any | None
        Existing axes/figure.
    **layout_kwargs
        Layout parameters.

    Returns
    -------
    Any
        Backend-specific return type.
    """
    # Get backend class
    BackendClass = get_backend(backend)

    # Create backend instance
    backend_instance = BackendClass()

    # Create or use existing axes
    if ax is not None:
        backend_instance._axes = ax
        if hasattr(ax, 'figure'):
            backend_instance._figure = ax.figure
    else:
        backend_instance.create_figure()
        backend_instance.create_axes()

    # Get or create style
    if style is None:
        style = default_style()

    # Compute layout
    computed_layout = compute_graph_layout(graph, layout=layout, **layout_kwargs)

    # Compute edge routes
    edge_routes = compute_graph_routes(graph, computed_layout.positions)

    # Create a simple layout-like object for backend compatibility
    # We'll create a minimal object that has the needed attributes
    class SimpleLayout:
        def __init__(self, graph: Any, positions: dict, edge_routes: dict):
            self.network = graph
            self.positions = positions
            self.edge_routes = edge_routes

    simple_layout = SimpleLayout(graph, computed_layout.positions, edge_routes)

    # Group parallel edges for offset calculation
    base_parallel_offset = 0.15
    parallel_offset_step = 0.1
    parallel_groups: dict[tuple[Any, Any], list[int]] = {}
    for (u, v, key), route in edge_routes.items():
        if route.edge_type.is_parallel:
            edge_key = (u, v) if route.edge_type.is_directed else (min(u, v), max(u, v))
            if edge_key not in parallel_groups:
                parallel_groups[edge_key] = []
            parallel_groups[edge_key].append(key)

    # Render all edges
    for (u, v, key), route in edge_routes.items():
        # Calculate parallel edge offset if needed
        parallel_offset = 0.0
        if route.edge_type.is_parallel:
            edge_key = (u, v) if route.edge_type.is_directed else (min(u, v), max(u, v))
            if edge_key in parallel_groups:
                parallel_keys = sorted(parallel_groups[edge_key])
                key_index = parallel_keys.index(key)
                direction = 1 if key_index % 2 == 0 else -1
                magnitude = base_parallel_offset + (key_index // 2) * parallel_offset_step
                parallel_offset = direction * magnitude

        backend_instance.render_edge(route, style, simple_layout, parallel_offset)

    # Render nodes
    # For graphs, all nodes are treated the same (no root/leaf/hybrid distinction)
    for node, position in computed_layout.positions.items():
        backend_instance.render_node(node, position, 'tree', style)

        # Add labels if enabled
        if style.with_labels:
            # Try to get label from graph node attributes
            label = None
            if hasattr(graph, '_graph'):
                # DirectedMultiGraph
                if node in graph._graph.nodes:
                    node_attrs = graph._graph.nodes[node]
                    label = node_attrs.get('label', str(node))
            elif hasattr(graph, '_directed'):
                # MixedMultiGraph
                if node in graph._directed.nodes:
                    node_attrs = graph._directed.nodes[node]
                    label = node_attrs.get('label', str(node))
                elif node in graph._undirected.nodes:
                    node_attrs = graph._undirected.nodes[node]
                    label = node_attrs.get('label', str(node))
            
            if label:
                backend_instance.add_label(position, label, style, 'tree')

    # Set axis properties (matplotlib-specific)
    if hasattr(backend_instance._axes, 'set_aspect'):
        backend_instance._axes.set_aspect('equal')
        backend_instance._axes.axis('off')

    # Return appropriate object based on backend
    if backend == 'pyqtgraph':
        return backend_instance._window
    else:
        return backend_instance._axes

