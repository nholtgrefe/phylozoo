"""
Public API for DirectedPhyNetwork plotting.

This module provides the main plotting function for users.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from ..utils.types import EdgeRoute

if TYPE_CHECKING:
    from phylozoo.core.network.dnetwork import DirectedPhyNetwork

from .backends import get_backend
from .layout import compute_dag_layout
from .styling import NetworkStyle, default_style


def plot_network(
    network: 'DirectedPhyNetwork',
    layout: str = 'dag',
    style: NetworkStyle | None = None,
    backend: str = 'matplotlib',
    ax: Any | None = None,
    **layout_kwargs: Any,
) -> Any:
    """
    Plot a DirectedPhyNetwork.

    This is the main public API function for plotting networks. It handles
    layout computation, styling, and rendering through the specified backend.

    Parameters
    ----------
    network : DirectedPhyNetwork
        The network to plot.
    layout : str, optional
        Layout algorithm. Currently only 'dag' is available.
        By default 'dag'.
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
        Additional parameters for layout computation (e.g., node_gap,
        layer_gap, iterations, curve_strength, etc.).

    Returns
    -------
    Any
        Backend-specific return type (e.g., matplotlib Axes).

    Raises
    ------
    ValueError
        If layout algorithm is not supported or backend is not registered.

    Examples
    --------
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> from phylozoo.visualize import plot_network
    >>>
    >>> net = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> ax = plot_network(net)
    """
    # Validate layout
    if layout != 'dag':
        raise ValueError(
            f"Layout '{layout}' is not supported. "
            "Currently only 'dag' layout is available."
        )

    # Get backend class
    BackendClass = get_backend(backend)

    # Create backend instance
    backend_instance = BackendClass()

    # Create or use existing axes
    if ax is not None:
        backend_instance._axes = ax
        backend_instance._figure = ax.figure
    else:
        backend_instance.create_figure()
        backend_instance.create_axes()

    # Get or create style
    if style is None:
        style = default_style()

    # Compute layout
    computed_layout = compute_dag_layout(network, **layout_kwargs)

    # Render edges
    network_obj = computed_layout.network
    positions = computed_layout.positions

    # Group parallel edges for offset calculation
    base_parallel_offset = 0.15
    parallel_offset_step = 0.1
    parallel_groups: dict[tuple[Any, Any], list[int]] = {}
    for (u, v, key), route in computed_layout.edge_routes.items():
        if route.edge_type.is_parallel:
            if (u, v) not in parallel_groups:
                parallel_groups[(u, v)] = []
            parallel_groups[(u, v)].append(key)

    # Render all edges
    for (u, v, key), route in computed_layout.edge_routes.items():
        # Calculate parallel edge offset if needed
        parallel_offset = 0.0
        if route.edge_type.is_parallel and (u, v) in parallel_groups:
            parallel_keys = sorted(parallel_groups[(u, v)])
            key_index = parallel_keys.index(key)
            direction = 1 if key_index % 2 == 0 else -1
            magnitude = base_parallel_offset + (key_index // 2) * parallel_offset_step
            parallel_offset = direction * magnitude

        backend_instance.render_edge(route, style, computed_layout, parallel_offset)

    # Render nodes
    leaves = network_obj.leaves
    hybrid_nodes = network_obj.hybrid_nodes
    root = network_obj.root_node

    for node, position in positions.items():
        if node == root:
            node_type = 'root'
        elif node in leaves:
            node_type = 'leaf'
        elif node in hybrid_nodes:
            node_type = 'hybrid'
        else:
            node_type = 'tree'

        backend_instance.render_node(node, position, node_type, style)

        # Add labels if enabled
        if style.with_labels:
            label = network_obj.get_label(node)
            if label:
                backend_instance.add_label(position, label, style, node_type)

    # Set axis properties (matplotlib-specific)
    if hasattr(backend_instance._axes, 'set_aspect'):
        backend_instance._axes.set_aspect('equal')
        backend_instance._axes.axis('off')

    # Return appropriate object based on backend
    if backend == 'pyqtgraph':
        return backend_instance._window
    else:
        return backend_instance._axes


def plot_tree(
    tree: 'DirectedPhyNetwork',
    layout: str = 'dag',
    style: NetworkStyle | None = None,
    backend: str = 'matplotlib',
    ax: Any | None = None,
    **layout_kwargs: Any,
) -> Any:
    """
    Plot a phylogenetic tree (DirectedPhyNetwork that is a tree).

    This is a convenience function that calls plot_network with tree-specific defaults.
    It's functionally equivalent to plot_network but provides a clearer API for trees.

    Parameters
    ----------
    tree : DirectedPhyNetwork
        The phylogenetic tree to plot (must be a tree, not a network).
    layout : str, optional
        Layout algorithm. Currently only 'dag' is available.
        By default 'dag'.
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
        Additional parameters for layout computation.

    Returns
    -------
    Any
        Backend-specific return type (e.g., matplotlib Axes).

    Examples
    --------
    >>> from phylozoo.core.network.dnetwork import DirectedPhyNetwork
    >>> from phylozoo.visualize import plot_tree
    >>>
    >>> tree = DirectedPhyNetwork(
    ...     edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> ax = plot_tree(tree)
    """
    return plot_network(tree, layout=layout, style=style, backend=backend, ax=ax, **layout_kwargs)
