"""
Public API for SemiDirectedPhyNetwork plotting.

This module provides the main plotting function for users.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Any

from phylozoo.core.network.sdnetwork.classifications import is_tree
from phylozoo.utils.exceptions import PhyloZooBackendError, PhyloZooLayoutError

from ..dnetwork.backends import get_backend
from ..dnetwork.styling import NetworkStyle, default_style
from ..utils.types import EdgeRoute
from .layout import compute_radial_layout

if TYPE_CHECKING:
    from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork


def plot_network(
    network: 'SemiDirectedPhyNetwork',
    layout: str = 'radial',
    style: NetworkStyle | None = None,
    backend: str = 'matplotlib',
    ax: Any | None = None,
    show: bool = False,
    **layout_kwargs: Any,
) -> Any:
    """
    Plot a SemiDirectedPhyNetwork tree using radial layout.

    This function plots semi-directed phylogenetic networks that are trees
    using a circular/radial layout. The root is positioned at the center
    and leaves are arranged on the outer circle.

    Parameters
    ----------
    network : SemiDirectedPhyNetwork
        The network to plot. Must be a tree (use is_tree() to check).
    layout : str, optional
        Layout algorithm. Currently only 'radial' is available.
        By default 'radial'.
    style : NetworkStyle, optional
        Styling configuration. If None, uses default style.
        By default None.
    backend : str, optional
        Backend to use ('matplotlib' or 'pyqtgraph').
        By default 'matplotlib'.
    ax : Any, optional
        Existing axes/figure to plot on (backend-specific).
        By default None.
    show : bool, optional
        If True, automatically display the plot using the backend's show() method.
        By default False.
    **layout_kwargs
        Additional parameters for layout computation (e.g., radius,
        start_angle, angle_direction).

    Returns
    -------
    Any
        Backend-specific return type (e.g., matplotlib Axes).

    Raises
    ------
    PhyloZooLayoutError
        If layout algorithm is not supported or network is not a tree.
    PhyloZooBackendError
        If backend is not registered.

    Examples
    --------
    >>> from phylozoo.core.network.sdnetwork import SemiDirectedPhyNetwork
    >>> from phylozoo.viz.sdnetwork import plot_network
    >>>
    >>> net = SemiDirectedPhyNetwork(
    ...     undirected_edges=[(3, 1), (3, 2)],
    ...     nodes=[(1, {'label': 'A'}), (2, {'label': 'B'})]
    ... )
    >>> ax = plot_network(net)
    """
    # Validate layout
    if layout != 'radial':
        raise PhyloZooLayoutError(
            f"Layout '{layout}' is not supported. "
            "Currently only 'radial' layout is available for SemiDirectedPhyNetwork."
        )

    # Check if network is a tree
    if not is_tree(network):
        raise PhyloZooLayoutError(
            "Radial layout is only supported for trees. "
            "Use is_tree() to check if the network is a tree."
        )

    # Get backend class (reuse from dnetwork)
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
    computed_layout = compute_radial_layout(network, **layout_kwargs)

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

        # Backend expects DAGLayout but RadialLayout has same interface (network, positions)
        # The backend only uses layout.network and layout.positions, so RadialLayout works
        backend_instance.render_edge(route, style, computed_layout, parallel_offset)  # type: ignore[arg-type]

    # Render nodes
    # Convert to directed network to get node types
    # Suppress warnings from intermediate networks during conversion
    import warnings
    from phylozoo.utils.exceptions import (
        PhyloZooEmptyNetworkWarning,
        PhyloZooSingleNodeNetworkWarning,
    )
    from phylozoo.core.network.sdnetwork.derivations import to_d_network
    
    # Only render nodes that are in the original network
    # (to_d_network may create subdivision nodes that we don't want to render)
    original_nodes = set(network._graph.nodes)
    
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=PhyloZooEmptyNetworkWarning)
        warnings.filterwarnings('ignore', category=PhyloZooSingleNodeNetworkWarning)
        d_network = to_d_network(network)
    leaves = d_network.leaves
    hybrid_nodes = d_network.hybrid_nodes
    root = d_network.root_node

    # Only render nodes that exist in the original network
    for node, position in positions.items():
        # Skip nodes that aren't in the original network (e.g., subdivision nodes)
        if node not in original_nodes:
            continue
            
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
            label = network.get_label(node)
            if label:
                # For radial layout, position labels radially outward from center
                x, y = position
                if node_type == 'leaf':
                    # Leaves: position label further out from center
                    angle = math.atan2(y, x)
                    label_x = x + style.label_offset * 1.5 * math.cos(angle)
                    label_y = y + style.label_offset * 1.5 * math.sin(angle)
                    # Use backend's add_label but with adjusted position
                    backend_instance.add_label((label_x, label_y), label, style, node_type)
                else:
                    # Internal nodes: use default positioning
                    backend_instance.add_label(position, label, style, node_type)

    # Set axis properties (matplotlib-specific)
    if hasattr(backend_instance._axes, 'set_aspect'):
        backend_instance._axes.set_aspect('equal')
        backend_instance._axes.axis('off')

    # Show plot if requested
    if show:
        backend_instance.show()

    # Return appropriate object based on backend
    if backend == 'pyqtgraph':
        return backend_instance._window
    else:
        return backend_instance._axes
